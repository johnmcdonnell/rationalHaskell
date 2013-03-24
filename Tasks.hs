{-# LANGUAGE DeriveDataTypeable #-}

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
import Data.Maybe (catMaybes, fromJust, fromMaybe)
import Data.Vector (freeze, thaw, MVector, Vector, (!))
import qualified Data.Vector as V
import qualified Data.Vector.Algorithms.Intro as VI
import qualified Numeric.LinearAlgebra as LA
import Control.Applicative
import Control.Monad
import Control.Monad.Random
import System.Random.Shuffle
import Statistics.Sample
import System.Console.CmdArgs

-- Internal libraries
import Utils
import Random
import Types
import Rational


-- Encapsulates the task and its prior.
type Task = (Stims, [PDFFromSample])

-- Generates priors given a task, the confidences, and a noise factor.
getPriorsRandom :: RandomGen g => (Double, Double, Double) -> [[Maybe Double]] -> Rand g [PDFFromSample]
getPriorsRandom (contalpha, contlambda, noise) stims  = do 
    sds <- sequence $ replicate n (lognormalM (log  pooledsd) noise)
    let tpriors = map (\(mu, sd) -> tPosterior (mu, sd, contalpha, contlambda)) (zip means sds)
    return $ tpriors ++ [binomprior]
  where
    means = map mean $ LA.toColumns stimmat
    pooledsd = stdDev $ LA.flatten stimmat
    stimmat = LA.fromLists $ map ((map fromJust) . init) stims
    binomprior =  bernoulliPosterior [1, 1]
    n = (length $ head stims) - 1

-- Generates priors given a task, the confidences, and a dimension bias.
-- Only works in the case of two non-label dimensions, unfortunately
-- var = sigma1 * sigma2 (this is fixed: pooledsd^2)
-- bias = sigma1 / sigma2 log(bias)
getPriorsBias :: (Double, Double, Maybe Double, Double) -> [[Maybe Double]] -> [PDFFromSample]
getPriorsBias (a0, lambda0, sigma0, logbias) stims  = tpriors ++ [binomprior]
  where
    sqrtbias = sqrt $ exp logbias
    tpriors = [tPosterior (mean dim, sigma, a0, lambda0) | (dim, sigma) <- zip (LA.toColumns stimmat) sigmas]
    sigmas = [pooledsd*sqrtbias, pooledsd/sqrtbias]
    pooledsd = fromMaybe ((stdDev . LA.flatten) stimmat) sigma0
    stimmat = LA.fromLists $ map ((map fromJust) . init) stims
    binomprior =  bernoulliPosterior [1, 1]

medinSchafferTask :: [Double] -> Task
medinSchafferTask bernoulliAlphas = (medinSchafferStims, andersondists)
  where 
    andersondists = replicate 5 bernoulli_prior
    bernoulli_prior = bernoulliPosterior bernoulliAlphas
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

zeithamovaMaddoxStims :: (RandomGen g) => Int -> Rand g [[Maybe Double]]
zeithamovaMaddoxStims n = do
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
    return $ map (map mabify) shuffled

zeithamovaMaddox :: (RandomGen g) => (Double, Double, Maybe Double, Double) -> Int -> Rand g Task
zeithamovaMaddox priorparams n = do
    stims <- zeithamovaMaddoxStims n
    priors <- return $ getPriorsBias priorparams stims
    let vectorizedstims = V.fromList $ (map V.fromList) $ stims
    return (vectorizedstims, priors)

randomInSquare :: (RandomGen g, Fractional a, Random a) => (a, a) -> (a, a) -> Rand g [a]
randomInSquare xbounds ybounds = do
    x <- getRandomR xbounds
    y <- getRandomR ybounds
    return [x, y]

randomsInSquare :: (RandomGen g, Fractional a, Random a) => (a, a) -> (a, a) -> Int -> Rand g [[a]]
randomsInSquare xbounds ybounds n = sequence $ replicate n (randomInSquare xbounds ybounds)

mcdonnellTaskStims :: RandomGen g => Int -> Int -> Rand g [[Maybe Double]]
mcdonnellTaskStims n nlab = do
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
    return shuffledstims

mcdonnellTask :: RandomGen g => (Double, Double, Maybe Double, Double) -> Int -> Int -> Rand g Task
mcdonnellTask priorparams n nlab = do
    stims <- mcdonnellTaskStims n nlab
    priors <- return $ getPriorsBias priorparams stims
    return (V.fromList $ (map V.fromList) $ stims, priors)

labfirstcompare :: Maybe a -> Maybe a -> Ordering
labfirstcompare (Just _)  (Nothing) = LT
labfirstcompare (Nothing) (Just _)  = GT
labfirstcompare _         _         = EQ

lablastcompare :: Maybe a -> Maybe a -> Ordering
lablastcompare (Nothing) (Just _)  = LT
lablastcompare (Just _)  (Nothing) = GT
lablastcompare _         _         = EQ

data SortOrder = Interspersed | LabeledFirst | LabeledLast deriving (Show, Data, Typeable)

mcdonnellTaskOrdered :: RandomGen g => SortOrder -> (Double, Double, Maybe Double, Double) -> Int -> Int -> Rand g Task
mcdonnellTaskOrdered Interspersed priors n nlab = mcdonnellTask priors n nlab
mcdonnellTaskOrdered order priorparams n nlab = do
    (task, priors) <- mcdonnellTask priorparams n nlab
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
    let withlabels = LA.fromBlocks [[LA.takeRows perCat xstims, scal2mat 0], 
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


