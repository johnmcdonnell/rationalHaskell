{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Types (Stim,
              Stims,
              Partition,
              validatePartition,
              clusterItems,
              PDFFromSample,
              ClusterPrior,
              bernoulliPosterior,
              tPosterior,
              dirichletProcess,
              Params,
              Model (..)
              ) where

import Data.Maybe (Maybe, fromJust, isJust)
import qualified Data.Vector as V
import qualified Data.Map as Map
import Data.List (group, sort, nub)
import Control.Monad.Random
import Control.Monad.Reader
import Control.Monad.State

import qualified Statistics.Distribution as StatDist
import Statistics.Distribution.StudentT as StatDist.T
import Statistics.Sample

import Utils

-- * Stim datatype and associated functions
type Stim = V.Vector (Maybe Double)

-- * Stim datatype and associated functions
type Stims = V.Vector Stim

-- * Partition datatype and associated functions
type Partition = V.Vector (Maybe Int)

-- | Partitions must have all members from [0..max]
validatePartition :: Partition -> Bool
validatePartition part 
  | V.null justs = True
  | otherwise    = minzero && correctn
  where
    correctn = length uniques == (max+1)
    minzero = min==0
    uniques = (nub . V.toList) justs
    max = V.maximum justs
    min = V.minimum justs
    justs = V.map fromJust $ V.filter isJust part

-- | Extract the items that have been assigned to the cluster.
clusterItems :: Partition -> Stims -> [Stims]
clusterItems assignments stims = takeWhile (not . V.null) clusters
  where
    clusters = map (\i-> V.map (stims V.!) $ V.elemIndices (Just i) assignments) [0..]

-- * Type of function to derive PDFs from a sample of points

type PDFFromSample = V.Vector Double -> (Maybe Double -> Double, Double)

fillerDistribtution :: PDFFromSample
fillerDistribtution sample = ((\x -> 0.5), 0.5)

-- | t posterior, used for continuous-valued dimensions
tPosterior :: (Double, Double, Double, Double) -> PDFFromSample
tPosterior (mu0, sigma0, alpha0, lambda0) sample  = (pdfFun, mui)
  where
    pdfFun (Just query) = StatDist.density (studentTUnstandardized alphai mui tsigma) query
    pdfFun Nothing = 1
    tsigma = sigmai * (sqrt (1 + (1/lambdai)))
    sigmai = sqrt $ (alpha0*(sq sigma0) + 
                      (n-1)*varianceunb +
                      (lambda0*n/lambdai)*(sq (mu0-xbar))
                    ) / alphai
    mui = (lambda0*mu0 + n * xbar) / lambdai
    (xbar, varianceunb) = meanVarianceUnb sample
    lambdai = lambda0 + n
    alphai = alpha0 + n
    n = fromIntegral $ V.length sample

-- | Bernoulli posterior, used for discrete-valued dimensions
bernoulliPosterior :: [Double] -> PDFFromSample
bernoulliPosterior alphas sample = (pdfFun, pdfFun (Just 1))
  where
    pdfFun Nothing = 1
    pdfFun (Just query) = (cj + alphas!!(floor query)) / (n + alpha0)
      where
        cj = (fromIntegral . V.length . (V.filter (==query))) sample
    ns = map (fromIntegral . snd) $ Map.toAscList $ countUnique $ V.toList $ sample
    n = fromIntegral $ V.length sample
    alpha0 = sum alphas

-- * Prior on partitions

type ClusterPrior = [Int] -> [Double]

-- Takes a list of assignments (assumed to be consecutive ints) and returns the
-- likelihood for each cluster, including a new cluster.
dirichletProcess :: Double -> ClusterPrior
dirichletProcess _ [] = [1]
dirichletProcess alpha assignments = pPartitions ++ [pNew]
  where
    pNew = if isInfinite alpha then 1 else alpha / (n+alpha)
    pPartitions = map (\x -> (fromIntegral x) / (n+alpha)) counts
    counts = map length $ (group . sort) assignments
    n = fromIntegral $ length assignments

type Params = (ClusterPrior, [PDFFromSample])

-- ^ Monad for model tracking state
newtype Model a = Model {
      runM :: ReaderT (Params, Stims) (StateT Partition (Rand StdGen)) a
          } deriving (Monad, MonadReader (Params, Stims), MonadState Partition)

