module Tasks (  Task
              , medinSchafferTask
              , onedtask
              , twodtask
              , zeithamovaMaddox
              , SortOrder (Interspersed, LabeledFirst, LabeledLast)
              , mcdonnellTaskOrdered
              , gridTest
             ) where

import Data.List (sortBy, transpose)
import Data.Maybe (catMaybes, fromJust)
import Data.Vector (freeze, thaw, MVector, Vector, (!))
import qualified Data.Vector as V
import qualified Data.Vector.Algorithms.Intro as VI
import Control.Applicative
import Control.Monad
import Control.Monad.Random
import System.Random.Shuffle
import Statistics.Sample

import Utils
import Random

import Rational
import Stats


-- Encapsulates the task and its prior.
type Task = (Stims, [PDFFromSample])

medinSchafferTask :: [Double] -> Task
medinSchafferTask binomAlphas = (medinSchafferStims, andersondists)
  where 
    andersondists = replicate 5 binom_prior
    binom_prior = bernoulliPosterior binomAlphas
    medinSchafferStims = V.fromList $ map (V.fromList . (map Just)) medinSchafferItems
    medinSchafferItems = [[1,1,1,1,1], 
                          [1,0,1,0,1], 
                          [1,0,1,1,0], 
                          [0,0,0,0,0], 
                          [0,1,0,1,1], 
                          [0,1,0,0,0]]
    
onedtask :: RandomGen g => [(Double, Double)] -> Int -> Rand g Task
onedtask params n = do
    samples <- forM params (\(mu, sigma) -> (take n) . (map Just) <$> (normalsM mu sigma))
    let stims = map (\x -> [x]) $ concat samples
    items <- shuffleM stims
    let tpriors = [tPosterior (mean itemvec, stdDev itemvec, 1, 1) | itemset <- (transpose . (map catMaybes)) items, let itemvec = V.fromList itemset  ]
    return (V.fromList $ map V.fromList items, tpriors)


twodtask :: (RandomGen g) => [(Double, Double)] -> Int -> Rand g Task
twodtask params n = do
    let mergeDims = (\(x,y) -> zipWith (\x y -> [x,y]) x y)
    stims <- forM params (\(mu, sigma) -> (mergeDims . (splitAt n) . (take $ n*2) . map Just) <$> (normalsM mu sigma))
    items <- shuffleM (concat stims)
    let tpriors = [tPosterior (mean itemvec, stdDev itemvec, 1, 1) | itemset <- (transpose . map catMaybes) items, let itemvec = V.fromList itemset  ]
    return (V.fromList $ map V.fromList items, tpriors)
 
zeithamovaMaddox :: (RandomGen g) => (Double, Double) -> Int -> Rand g Task
zeithamovaMaddox (contalpha, contlambda) n = do
    -- length bimodal
    let bimodmean1 = 187.5
    let bimodmean2 = 412.5
    let bimodsd = 12.5
    let unimodmean = 45
    let unimodsd = 15
    astims <-  map (map Just) <$> binormals bimodmean1 unimodmean bimodsd unimodsd n
    bstims <- map (map Just) <$> binormals bimodmean2 unimodmean bimodsd unimodsd n
    items <- shuffleM (astims ++ bstims)
    let itemswithlabels = V.fromList $ map (\x -> V.snoc (V.fromList x) Nothing) items
    let tpriors = [tPosterior (mean itemvec, stdDev itemvec, contalpha, contlambda) | itemset <- (transpose . (map catMaybes)) items , let itemvec = V.fromList itemset ]
    let binomprior =  bernoulliPosterior [1, 1]
    return (itemswithlabels, tpriors ++ [binomprior])


randomInSquare :: (RandomGen g, Fractional a, Random a) => (a, a) -> (a, a) -> Rand g [a]
randomInSquare xbounds ybounds = do
    x <- getRandomR xbounds
    y <- getRandomR ybounds
    return [x, y]

randomsInSquare :: (RandomGen g, Fractional a, Random a) => (a, a) -> (a, a) -> Int -> Rand g [[a]]
randomsInSquare xbounds ybounds n = sequence $ replicate n (randomInSquare xbounds ybounds)

mcdonnellTask :: RandomGen g => (Double, Double) -> Int -> Int -> Rand g Task
mcdonnellTask (contalpha, contlambda) n nlab = do
    let (boxmultiplier, nrem) = quotRem n 28
    let (perCat, nlabrem) = quotRem nlab 2
    let unlab = n - nlab
    when (nrem /= 0) $ error $ "McDonnell Task n must be divisible by 28. Instead it was " ++ show n 
    when (nlabrem /= 0) $ error $ "McDonnell Task nlab must be divisible by 2. Instead it was " ++ show nlab
    let bimodBounds = map (,,) [(0, 0.2), (0.8, 1)]
    let boxBorders = [0,0.2..1]
    let boxes = zip (init boxBorders) (tail boxBorders)
    let boxcounts = (map (*boxmultiplier) [2,3,4,3,2])
    let stimArgs = map uncurry bimodBounds <*> (zip boxes boxcounts)
    stimlocsGrouped <- sequence $ map (\(xs, ys, n) -> randomsInSquare xs ys n) stimArgs
    let stimlocs = concat stimlocsGrouped
    let labs = (replicate perCat $ Just 0) ++ (replicate unlab Nothing) ++ (replicate perCat $ Just 1)
    let stims =  zipWith (\loc maybelab -> (map Just loc) ++ [maybelab])  stimlocs labs
    shuffledstims <- shuffleM stims
    let tpriors = [tPosterior (mean itemvec, stdDev itemvec, contalpha, contlambda) | itemset <- transpose stimlocs, let itemvec = V.fromList itemset ]
    let binomprior =  bernoulliPosterior [1, 1]
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

mcdonnellTaskOrdered :: RandomGen g => SortOrder -> (Double, Double) -> Int -> Int -> Rand g Task
mcdonnellTaskOrdered Interspersed priors n nlab = mcdonnellTask priors n nlab
mcdonnellTaskOrdered order (contalpha, contlambda) n nlab = do
    (task, priors) <- mcdonnellTask (contalpha, contlambda) n nlab
    let sortedtask = V.modify tasksort task
    return (sortedtask, priors)
  where
    tasksort v = VI.sortBy (\x y -> comparison (V.last x) (V.last y)) v
    comparison = case order of LabeledFirst -> labfirstcompare
                               LabeledLast -> lablastcompare

gridTest :: Stims
gridTest = V.fromList $ (\x y -> V.fromList [Just x, Just y, Nothing]) <$> [ 0,(1/6)..1 ] <*> [ 0,(1/6)..1 ]


vandistTask :: RandomGen g => (Double, Double) -> Int -> Double -> Rand g Task
vandistTask (contalpha, contlambda) n proplab = do
    let nlab = floor $ proplab * (fromIntegral n)
    let (nperCat, nrem) = quotRem n 2
    let (perCat, nlabrem) = quotRem nlab 2
    let unlab = n - nlab
    when ((fromIntegral nlab) /= (proplab * (fromIntegral n))) $ error "n did not divid into proplab."
    when (nrem /= 0) $ error $ "n must be divisible by 2. Instead it was " ++ show n ++ "."
    when (nlabrem /= 0) $ error $ "# of labeled items must be divisible by 2. Instead it was " ++ show nlab ++ "."
    let meanX1 = 40
        meanX2 = 60
        meanN1 = 60
        meanN2 = 40
        sd = 11.88
        rho = 0.99
    xs <- map ((++[Just 0]) . (map Just)) <$> binormalscov meanX1 meanX2 sd sd rho nperCat
    ns <- map ((++[Just 1]) . (map Just)) <$> binormalscov meanN1 meanN2 sd sd rho nperCat
    let newxs = (take perCat xs) ++ (map (\[d1,d2,lab] -> [d1,d2,Nothing]) (drop perCat xs))
    let newns = (take perCat ns) ++ (map (\[d1,d2,lab] -> [d1,d2,Nothing]) (drop perCat ns))
    shuffledstims <- shuffleM (newxs ++ newns)
    let dim1 = map (!!0) shuffledstims
        dim2 = map (!!1) shuffledstims
    let tpriors = [tPosterior (mean itemvec, stdDev itemvec, contalpha, contlambda) | itemset <- [dim1, dim2], let itemvec = (V.fromList . (map fromJust)) itemset ]
        binomprior =  bernoulliPosterior [1, 1]
    return (V.fromList $ map V.fromList shuffledstims, tpriors ++ [binomprior])

