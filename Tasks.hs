module Tasks (  Task
              , medinSchafferTask
              , onedtask
              , twodtask
              , zeithamovaMaddox
              , SortOrder (Interspersed, LabeledFirst, LabeledLast)
              , mcdonnellTaskOrdered
              , gridTest
              , vandistTask
             ) where

import Data.List (sortBy, transpose)
import Data.Maybe (catMaybes, fromJust)
import Data.Vector (freeze, thaw, MVector, Vector, (!))
import qualified Data.Vector as V
import qualified Data.Vector.Algorithms.Intro as VI
import qualified Numeric.LinearAlgebra as LA
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

scal2mat :: LA.Element a => a -> (LA.Matrix a)
scal2mat x = LA.fromLists [[x]]

zeithamovaMaddox :: (RandomGen g) => (Double, Double) -> Int -> Rand g Task
zeithamovaMaddox (contalpha, contlambda) n = do
    -- length bimodal
    let meansA = LA.fromList [187.5, 45]
        meansB = LA.fromList [412.5, 45]
        sigmas = (12.5, 12.5)
        rho = 0
    astims <- binormals meansA sigmas rho n
    bstims <- binormals meansB sigmas rho n
    let withlabels = LA.fromBlocks [[astims, scal2mat $ 0/0], [bstims, scal2mat $ 0/0]]
    shuffled <- shuffleM $ LA.toLists withlabels
    let mabify = (\x -> if isNaN x then Nothing else Just x)
        stims = V.fromList $ map (V.fromList . (map mabify)) shuffled
    -- let stims = map ((V.snoc Nothing) . V.fromList . (map Just)) shuffled
    let dims = init $ LA.toColumns withlabels
        tpriors = [tPosterior (mean dim, stdDev dim, contalpha, contlambda) | dim <- dims]
    let binomprior =  bernoulliPosterior [1, 1]
    return (stims, tpriors ++ [binomprior])


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
    let meansX = LA.fromList [40, 60]
        meansN = LA.fromList [60, 40]
        sigmas = (11.88, 11.88)
        rho = 0.99
    xstims <- binormals meansX sigmas rho nperCat
    nstims <- binormals meansN sigmas rho nperCat
    let withlabels = LA.fromBlocks [[debug $ LA.takeRows perCat xstims, scal2mat 0], 
                                    [LA.dropRows perCat xstims, scal2mat $ 0/0],
                                    [LA.takeRows perCat nstims, scal2mat 1],
                                    [LA.dropRows perCat nstims, scal2mat $ 0/0]]
    shuffled <- shuffleM $ LA.toLists withlabels
    let mabify = (\x -> if isNaN x then Nothing else Just x)
        stims = V.fromList $ map (V.fromList . (map mabify)) shuffled
    let dims = init $ LA.toColumns withlabels
        tpriors = [tPosterior (mean dim, stdDev dim, contalpha, contlambda) | dim <- dims]
        binomprior =  bernoulliPosterior [1, 1]
    return (stims, tpriors ++ [binomprior])


