
-- import qualified Control.Monad.State as State
import Data.List (sortBy, transpose)
import Data.Function (on)
import Data.Maybe
import Data.Vector (freeze, thaw, MVector, Vector, (!))
import qualified Data.Vector as V
import qualified Data.Vector.Algorithms.Intro as VI
import Control.Applicative
import Control.Arrow
import Control.Monad
import Control.Monad.ST
import Control.Monad.Random
import System.Environment (getArgs)
import System.Random
import System.Random.Shuffle
import Text.CSV

import Statistics.Sample

import JohnFuns
import Random

import Stats
import Rational
import Anderson

-- Convenience functions


-- Science etc

medinSchafferTask :: [Double] -> (Stims, [PDFFromSample])
medinSchafferTask binomAlphas = (medinSchafferStims, andersondists)
  where 
    andersondists = replicate 5 binom_prior
    binom_prior = binomialPosterior binomAlphas
    medinSchafferStims = V.fromList $ map (V.fromList . (map Just)) medinSchafferItems
    medinSchafferItems = [[1,1,1,1,1], 
                          [1,0,1,0,1], 
                          [1,0,1,1,0], 
                          [0,0,0,0,0], 
                          [0,1,0,1,1], 
                          [0,1,0,0,0]]
    

testMedinSchaffer = do
    let (task, dists) = medinSchafferTask [1,1]
    let couplingParam = dirichletProcess 1.0
    andersonSample EncodeActual (couplingParam, dists) task

onedtask :: [(Double, Double)] -> Int -> IO (Stims, [PDFFromSample])
onedtask params n = do
    samples <- forM params (\(mu, sigma) -> evalRandIO $ (take n) . (map Just) <$> (normalsM mu sigma))
    let stims = map (\x -> [x]) $ concat samples
    items <- evalRandIO $ shuffleM stims
    let tpriors = [tPosterior (mean itemvec, stdDev itemvec, 1, 1) | itemset <- (transpose . (map catMaybes)) items, let itemvec = V.fromList itemset  ]
    return (V.fromList $ map V.fromList items, tpriors)


twodtask :: [(Double, Double)] -> Int -> IO (Stims, [PDFFromSample])
twodtask params n = do
    let mergeDims = (\(x,y) -> zipWith (\x y -> [x,y]) x y)
    stims <- forM params (\(mu, sigma) -> evalRandIO $ (mergeDims . (splitAt n) . (take $ n*2) . map Just) <$> (normalsM mu sigma))
    items <- evalRandIO $ shuffleM (concat stims)
    let tpriors = [tPosterior (mean itemvec, stdDev itemvec, 1, 1) | itemset <- (transpose . map catMaybes) items, let itemvec = V.fromList itemset  ]
    return (V.fromList $ map V.fromList items, tpriors)

zeithamovaMaddox :: (Double, Double) -> Int -> IO (Stims, [PDFFromSample])
zeithamovaMaddox (contalpha, contlambda) n = do
    -- length bimodal
    let bimodmean1 = 187.5
    let bimodmean2 = 412.5
    let bimodsd = 12.5
    let unimodmean = 45
    let unimodsd = 15
    astims <-  evalRandIO $ map (map Just) <$> binormals bimodmean1 unimodmean bimodsd unimodsd n
    bstims <- evalRandIO $ map (map Just) <$> binormals bimodmean2 unimodmean bimodsd unimodsd n
    items <- evalRandIO $ shuffleM (astims ++ bstims)
    let itemswithlabels = V.fromList $ map (\x -> V.snoc (V.fromList x) Nothing) items
    let tpriors = [tPosterior (mean itemvec, stdDev itemvec, contalpha, contlambda) | itemset <- (transpose . (map catMaybes)) items , let itemvec = V.fromList itemset ]
    let binomprior =  binomialPosterior [1, 1]
    return (itemswithlabels, tpriors ++ [binomprior])

randomInSquare :: (RandomGen g, Fractional a, Random a) => (a, a) -> (a, a) -> Rand g [a]
randomInSquare xbounds ybounds = do
    x <- getRandomR xbounds
    y <- getRandomR ybounds
    return [x, y]

randomsInSquare :: (RandomGen g, Fractional a, Random a) => (a, a) -> (a, a) -> Int -> Rand g [[a]]
randomsInSquare xbounds ybounds n = sequence $ replicate n (randomInSquare xbounds ybounds)

mcdonnellTask :: (Double, Double) -> Int -> Int -> IO (Stims, [PDFFromSample])
mcdonnellTask (contalpha, contlambda) n nlab = do
    let (boxmultiplier, nrem) = quotRem n 28
    let (perCat, nlabrem) = quotRem nlab 2
    let unlab = n - nlab
    when (nrem /=0 ) $ error $ "McDonnell Task n must be divisible by 28. Instead it was " ++ show n 
    when (nlabrem /= 0) $ error $ "McDonnell Task nlab must be divisible by 2. Instead it was " ++ show nlab
    let bimodBounds = map (,,) [(0, 0.2), (0.8, 1)]
    let boxBorders = [0,0.2..1]
    let boxes = zip (init boxBorders) (tail boxBorders)
    let boxcounts = (map (*boxmultiplier) [2,3,4,3,2])
    let stimArgs = map uncurry bimodBounds <*> (zip boxes boxcounts)
    stimlocsGrouped <- evalRandIO $ sequence $ map (\(xs, ys, n) -> randomsInSquare xs ys n) stimArgs
    let stimlocs = concat stimlocsGrouped
    let labs = (replicate perCat $ Just 0) ++ (replicate unlab Nothing) ++ (replicate perCat $ Just 1)
    let stims =  zipWith (\loc maybelab -> (map Just loc) ++ [maybelab])  stimlocs labs
    shuffledstims <- evalRandIO $ shuffleM stims
    let tpriors = [tPosterior (mean itemvec, stdDev itemvec, contalpha, contlambda) | itemset <- transpose stimlocs, let itemvec = V.fromList itemset ]
    let binomprior =  binomialPosterior [1, 1]
    return (V.fromList $ map V.fromList shuffledstims, tpriors ++ [binomprior])


labfirstcompare :: Ord a => Maybe a -> Maybe a -> Ordering
labfirstcompare (Just _)  (Nothing) = LT
labfirstcompare (Nothing) (Just _)  = GT
labfirstcompare _         _         = EQ

lablastcompare :: Ord a => Maybe a -> Maybe a -> Ordering
lablastcompare (Nothing) (Just _)  = LT
lablastcompare (Just _)  (Nothing) = GT
lablastcompare _         _         = EQ

data SortOrder = Interspersed | LabeledFirst | LabeledLast deriving (Show)

mcdonnellTaskOrdered :: SortOrder -> (Double, Double) -> Int -> Int -> IO (Stims, [PDFFromSample])
mcdonnellTaskOrdered Interspersed priors n nlab = mcdonnellTask priors n nlab
mcdonnellTaskOrdered order (contalpha, contlambda) n nlab = do
    (task, priors) <- mcdonnellTask (contalpha, contlambda) n nlab
    thawed <- V.unsafeThaw task
    VI.sortBy (\x y -> comparison (V.last x) (V.last y)) thawed
    V.unsafeFreeze thawed
    return (task, priors)
  where
    comparison = case order of LabeledFirst -> labfirstcompare
                               LabeledLast -> lablastcompare

runGridTest :: ClusterPrior -> [PDFFromSample] -> Stims -> Partition -> Vector (Stim, Double)
runGridTest cprior distributions stimuli assignments = V.zip grid labels
  where 
    labels = V.map (head . (infer cprior distributions stimuli assignments)) grid
    grid = V.fromList $ (\x y -> V.fromList [Just x, Just y, Nothing]) <$> [ 0,(1/6)..1 ] <*> [ 0,(1/6)..1 ]


testContinuous = do
    let mu1 = 1
    let sigma1 = 1
    let mu2 = -1
    let sigma2 = 1
    let n = 3
    putStrLn "Testing on:" 
    putStrLn $ "mu1=" ++ (show mu1) ++ " sigma1=" ++ (show sigma1) ++ " n=" ++ (show n)
    putStrLn $ "mu2=" ++ (show mu2) ++ " sigma2=" ++ (show sigma2) ++ " n=" ++ (show n)
    (task, distpriors) <- twodtask [(mu1, sigma1), (mu2, sigma2)] n
    
    let prior  = (dirichletProcess 1.0, distpriors)
    partition <- evalRandIO $ andersonSample EncodeActual prior task
    V.forM_ (V.zip task (V.map (fromMaybe (-1)) partition)) print

testZeithamova = do
    (task, distpriors) <- zeithamovaMaddox (1, 1) 100
    
    let prior  = (dirichletProcess 1.0, distpriors)
    partition <- evalRandIO $ andersonSample EncodeActual prior task
    V.forM_ (V.zip task (V.map (fromMaybe (-1)) partition)) print

testTVTask = do
    args <- getArgs
    let cparam = if length args > 0 then read (args!!0) else 1
    let maxlab = 16
    let (nlab, nounlab) = if length args > 1 then (\x -> if x<0 then (maxlab, True) else (x, False) ) $ read (args!!1) 
                                             else (maxlab, False)
    let orderarg = if length args > 2 then (args!!2) else "interspersed"
    let encodearg = if length args > 3 then (args!!3) else "actual"
    
    let order = case orderarg of "interspersed" -> Interspersed
                                 "labfirst"     -> LabeledFirst
                                 "lablast"      -> LabeledLast
                                 otherwise      -> error $ "Inappropriate order: " ++ orderarg ++ "; order should be one of interspersed, labfirst, lablast."
    let encoding = case encodearg of "guess" -> EncodeGuess
                                     "softguess" -> EncodeGuessSoft
                                     otherwise -> EncodeActual
    
    -- Set up priors
    let filterfun = if nounlab then V.filter (isJust . V.last) else id
    (task, distpriors) <- first filterfun <$> mcdonnellTaskOrdered order (1, 1) (28*4) nlab
    -- print $ sortBy (compare `on` fst) $ map (\((Just bimod):(Just unimod):xs) -> (bimod, unimod)) task
    
    let prior  = (dirichletProcess cparam, distpriors)
    
    -- Now run the model
    partition <- evalRandIO $ andersonSample encoding prior task
    
    -- Dump all those mabies
    let demabify = V.map (fromMaybe (-9))
    let demabified = V.map demabify task
    let results = map (\(x,Just y) -> V.toList $ V.snoc x (fromIntegral y)) $ (filter (isJust . snd)) $ zip (V.toList demabified) (V.toList partition)
    let resultstrings = map (map show) results
    
    -- Print it out
    putStrLn $ printCSV $ map ("STIM":) resultstrings
    
    putStrLn $ printCSV $ map (("CLUST":) . (map show)) $ summarizeClusters distpriors task partition
    
    -- Get inference at each point in a grid
    putStrLn $ printCSV $ map (("INFER":) . (map show) . \(stim, resp) -> (take 2 . V.toList . demabify) stim ++ [resp] ) $ V.toList $ (uncurry runGridTest) prior task partition

main = do
    -- testMedinSchaffer
    -- testContinuous
    -- testZeithamova
    testTVTask

