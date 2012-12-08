{-# LANGUAGE Rank2Types #-}

module Rational (clusterPosterior,
                 clusterItems,
                 infer,
                 Stim,
                 Stims,
                 Partition,
                 validatePartition,
                 clusterPredictions,
                 summarizeClusters
                ) where

import Data.Function (on)
import Data.Maybe
import Data.List (nub, group, groupBy, transpose, sortBy, splitAt)
import qualified Data.Set as Set
import Data.Vector (freeze, thaw, MVector, Vector, (!))
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM

import Control.Exception
import Control.Monad
import Control.Monad.ST

import Statistics.Sample

import Utils
import Stats

-- * Convenience functions
-- * Stim dattype and associated functions
type Stim = Vector (Maybe Double)

-- * Stim dattype and associated functions
type Stims = Vector Stim

-- * Partition dattype and associated functions
type Partition = Vector (Maybe Int)

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


clusterItems :: Partition -> Stims -> [Stims]
clusterItems assignments stims = takeWhile (not . V.null) clusters
  where
    clusters = map (\i-> V.map (stims!) $ V.elemIndices (Just i) assignments) [0..]

-- | Predictions for each cluster
clusterPredictions :: [PDFFromSample] -> Stims -> [Double]
clusterPredictions pdfs stims = zipWith distSummary pdfs dims
  where
    distSummary dist = snd . dist . catMaybeV
    dims = V.toList $ vectranspose ndims stims 
    ndims = V.length $ V.head stims

-- | Likelihood that a stimulus belongs to each cluster
clusterPosterior :: ClusterPrior -> [PDFFromSample] -> Stims -> Partition -> Stim -> [Double]
clusterPosterior cprior distributions stimuli assignments newstim = map (/norm) posterior
  where
    norm = sum posterior
    posterior = zipWith (*) clustPriors clustLikelihoods
    -- NOTE: watch out that the stim to be assigned isn't assigned in assignments.
    clustPriors = cprior $ (catMaybes . V.toList) assignments
    clustLikelihoods = map clustLik clusts ++ [emptyClustLik]
    clustLik clust = V.product $ V.zipWith3 (\dist sample query -> (fst . dist) (V.map fromJust $ (V.filter isJust) sample) query ) distv (vectranspose ndim clust) newstim
    emptyClustLik = V.product $ V.zipWith (\dist query -> (fst . dist) V.empty query) distv newstim
    clusts = clusterItems assignments stimuli
    stimlength = (V.length . V.head) stimuli
    distv = V.fromList distributions
    ndim = length distributions

-- | For all missing items in the querystim, infers the most likely possible
-- | outcome, taking all possible cluster assignments into account.
infer :: (ClusterPrior, [PDFFromSample]) -> Stims -> Partition -> Stim -> ([Double], [Double])
infer (cprior, distributions) stimuli assignments querystim = (prediction, post)
  where
    prediction = map inferDim queryDims
    inferDim i = sum $ (zipWith (*) post) (clustPreds!!i)
    clustPreds = transpose $ map (clusterPredictions distributions) clusters
    clusters = clusterItems assignments stimuli
    queryDims = V.toList $ V.map snd $ V.filter (isNothing . fst) $ V.zip querystim (V.enumFromN 0 (V.length querystim))
    post = clusterPosterior cprior distributions stimuli assignments querystim


-- | Cluster centroids for all clusters in all dimensions
summarizeClusters :: [PDFFromSample] -> Stims -> Partition -> [[Double]]
summarizeClusters pdfs stims partition = map summfun clusts
  where
    summfun = clusterPredictions pdfs
    clusts = clusterItems partition stims

