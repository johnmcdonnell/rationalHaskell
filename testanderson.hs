
-- import qualified Control.Monad.State as State
import Data.List (sortBy, transpose)
import Data.Function (on)
import Data.Maybe
import Control.Applicative
import Control.Arrow
import Control.Monad
import Control.Monad.Random
import System.Environment (getArgs)
import System.Random
import System.Random.Shuffle
import Text.CSV

import Math.Statistics

import Graphics.Gnuplot.Simple

import JohnFuns
import Random

import Stats
import Rational
import Anderson

-- Convenience functions

-- plotpoints :: [(Double, Double)] -> IO ()
-- plotpoints = plotListStyle [] defaultStyle{ plotType=Dots }


-- Science etc

medinSchafferTask :: [Double] -> ([Stim], [PDFFromSample])
medinSchafferTask binomAlphas = (medinSchafferStims, andersondists)
  where 
    andersondists = replicate 5 binom_prior
    binom_prior = binomialPosterior binomAlphas
    medinSchafferStims = map (map Just) medinSchafferItems
    medinSchafferItems = [[1,1,1,1,1], 
                          [1,0,1,0,1], 
                          [1,0,1,1,0], 
                          [0,0,0,0,0], 
                          [0,1,0,1,1], 
                          [0,1,0,0,0]]
    

testMedinSchaffer = do
    let (task, dists) = medinSchafferTask [1,1]
    let couplingParam = dirichletProcess 1.0
    andersonSample (couplingParam, dists) task

onedtask :: [(Double, Double)] -> Int -> IO ([Stim], [PDFFromSample])
onedtask params n = do
    samples <- forM params (\(mu, sigma) -> evalRandIO $ (take n) . (map Just) <$> (normalsM mu sigma))
    let stims = map (\x -> [x]) $ concat samples
    items <- evalRandIO $ shuffleM stims
    let tpriors = [tPosterior (mean itemset, stddev itemset, 1, 1) | itemset <- (transpose . (map catMaybes)) items ]
    return (items, tpriors)


twodtask :: [(Double, Double)] -> Int -> IO ([Stim], [PDFFromSample])
twodtask params n = do
    let mergeDims = (\(x,y) -> zipWith (\x y -> [x,y]) x y)
    stims <- forM params (\(mu, sigma) -> evalRandIO $ (mergeDims . (splitAt n) . (take $ n*2) . map Just) <$> (normalsM mu sigma))
    items <- evalRandIO $ shuffleM (concat stims)
    let tpriors = [tPosterior (mean itemset, stddev itemset, 1, 1) | itemset <- (transpose . map catMaybes) items ]
    return (items, tpriors)

tvLabelFun :: Stim -> Stim
tvLabelFun [ Just bimod, Just unimod ] 
  | bimod < 300 && unimod > 55 = [ Just bimod, Just unimod, Just 0  ]
  | bimod > 300 && unimod < 35 = [ Just bimod, Just unimod, Just 1  ]
  | otherwise                  = [ Just bimod, Just unimod, Nothing ]

tvTask :: (Double, Double) -> Int -> IO ([Stim], [PDFFromSample])
tvTask (contalpha, contlambda) n = do
    -- length bimodal
    let bimodmean1 = 187.5
    let bimodmean2 = 412.5
    let bimodsd = 12.5
    let unimodmean = 45
    let unimodsd = 15
    astims <- evalRandIO $  (map (map Just)) <$> binormals bimodmean1 unimodmean bimodsd unimodsd n
    bstims <- evalRandIO $ (map (map Just)) <$> binormals bimodmean2 unimodmean bimodsd unimodsd n
    items <- evalRandIO $ shuffleM (astims ++ bstims)
    let itemswithlabels = map tvLabelFun items
    let tpriors = [tPosterior (0.5, 0.2, contalpha, contlambda) | itemset <- (transpose . (map catMaybes)) items ]
    let binomprior =  binomialPosterior [1, 1]
    return (itemswithlabels, tpriors ++ [binomprior])
    

zeithamovaMaddox :: (Double, Double) -> Int -> IO ([Stim], [PDFFromSample])
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
    let itemswithlabels = map (\x -> x ++ [Nothing]) items
    let tpriors = [tPosterior (mean itemset, stddev itemset, contalpha, contlambda) | itemset <- (transpose . (map catMaybes)) items ]
    let binomprior =  binomialPosterior [1, 1]
    return (itemswithlabels, tpriors ++ [binomprior])

randomInSquare :: (RandomGen g, Fractional a, Random a) => (a, a) -> (a, a) -> Rand g [a]
randomInSquare xbounds ybounds = do
    x <- getRandomR xbounds
    y <- getRandomR ybounds
    return [x, y]

randomsInSquare :: (RandomGen g, Fractional a, Random a) => (a, a) -> (a, a) -> Int -> Rand g [[a]]
randomsInSquare xbounds ybounds n = sequence $ replicate n (randomInSquare xbounds ybounds)

mcdonnellTask :: (Double, Double) -> Int -> Int -> IO ([Stim], [PDFFromSample])
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
    let tpriors = [tPosterior (mean itemset, stddev itemset, contalpha, contlambda) | itemset <- transpose stimlocs ]
    let binomprior =  binomialPosterior [1, 1]
    return (shuffledstims, tpriors ++ [binomprior])


labfirstcompare :: Ord a => Maybe a -> Maybe a -> Ordering
labfirstcompare (Just _)  (Nothing) = LT
labfirstcompare (Nothing) (Just _)  = GT
labfirstcompare _         _         = EQ

lablastcompare :: Ord a => Maybe a -> Maybe a -> Ordering
lablastcompare (Nothing) (Just _)  = LT
lablastcompare (Just _)  (Nothing) = GT
lablastcompare _         _         = EQ

data SortOrder = Interspersed | LabeledFirst | LabeledLast deriving (Show)

mcdonnellTaskOrdered :: SortOrder -> (Double, Double) -> Int -> Int -> IO ([Stim], [PDFFromSample])
mcdonnellTaskOrdered order (contalpha, contlambda) n nlab = (first orderfun) <$> (mcdonnellTask (contalpha, contlambda) n nlab)
-- mcdonnellTaskOrdered order = (first (orderfun . last)) <$> mcdonnellTask
  where
    orderfun = case order of 
                      Interspersed -> id
                      LabeledFirst -> sortBy (\x y ->  labfirstcompare (last x) (last y))
                      LabeledLast  -> sortBy (\x y ->  lablastcompare (last x) (last y))

runGridTest :: ClusterPrior -> [PDFFromSample] -> [Stim] -> Partition -> [(Stim, Double)]
runGridTest cprior distributions stimuli assignments = zip grid labels
  where 
    labels = map head $ map (infer cprior distributions stimuli assignments) grid
    grid = (\x y -> [Just x, Just y, Nothing]) <$> [ 0,(1/6)..1 ] <*> [ 0,(1/6)..1 ]


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
    partition <- andersonSample prior task
    forM_ (zip task (catMaybes partition)) print

testZeithamova = do
    (task, distpriors) <- zeithamovaMaddox (1, 1) 100
    
    let prior  = (dirichletProcess 1.0, distpriors)
    partition <- andersonSample prior task
    forM_ (zip task (catMaybes partition)) print

testTVTask = do
    args <- getArgs
    let cparam = if length args > 0 then read (args!!0) else 1
    let nlab = if length args > 1 then read (args!!1) else 16
    let orderarg = if length args > 2 then (args!!2) else "interspersed"
    
    let order = case orderarg of "interspersed" -> Interspersed
                                 "labfirst"     -> LabeledFirst
                                 "lablast"      -> LabeledLast
                                 otherwise      -> error $ "Inappropriate order: " ++ orderarg ++ "; order should be one of interspersed, labfirst, lablast."
    print $ "Order was " ++ show order;
    
    -- Set up priors
    (task, distpriors) <- mcdonnellTaskOrdered order (1, 1) (28*4) nlab
    -- print $ sortBy (compare `on` fst) $ map (\((Just bimod):(Just unimod):xs) -> (bimod, unimod)) task
    
    let prior  = (dirichletProcess cparam, distpriors)
    
    -- Now run the model
    partition <- andersonSample prior task
    
    -- Dump all those mabies
    let demabify [bimod,unimod,lab] = [fromJust bimod, fromJust unimod, fromMaybe (-9) lab]
    let demabified = map demabify task
    let results = map (\(x,Just y) -> x ++ [fromIntegral y]) $ (filter (isJust . snd)) $ zip demabified partition
    let resultstrings = map (map show) results
    
    -- Print it out
    putStrLn $ printCSV $ map ("STIM":) resultstrings
    
    -- Plot it
    --let grouped = map (map (\([x,y,_],_) -> (x,y))) $ groupsortBy snd results
    -- plotListsStyle [] $ map (\x -> (defaultStyle {plotType = Dots, lineSpec = CustomStyle [PointType 5]}, x)) grouped
    --plotListsStyle [LineStyle 1 [PointType 1, PointSize 0.8]] $ map (\x -> (defaultStyle {plotType = Dots}, x)) grouped
    
    -- Get inference at each point in a grid
    putStrLn $ printCSV $ map (("INFER":) . (map show) . \(stim, resp) -> (take 2 . demabify) stim ++ [resp] ) $ (uncurry runGridTest) prior task partition 

main = do
    -- testMedinSchaffer
    -- testContinuous
    -- testZeithamova
    testTVTask

