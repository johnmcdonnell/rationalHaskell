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
import Types

-- | Predictions for each cluster
clusterPredictions :: [PDFFromSample] -> Stims -> [Double]
clusterPredictions pdfs stims 
  | V.null stims = zipWith distSummary pdfs (repeat V.empty)
  | otherwise    = zipWith distSummary pdfs dims
  where
    distSummary dist = snd . dist . catMaybeV
    dims = V.toList $ vectranspose ndims stims 
    ndims = V.length $ V.head stims

-- | Likelihood that a stimulus belongs to each cluster
clusterPosterior :: ClusterPrior -> [PDFFromSample] -> Stims -> Partition -> Stim -> [Double]
clusterPosterior cprior distributions stimuli assignments newstim = map (/norm) posterior
  where
    norm = sum posterior
    -- NOTE: stick a debug in here to get the Anderson predictions.
    posterior = zipWith (*) clustPriors clustLikelihoods
    clustPriors = cprior $ (catMaybes . V.toList) assignments
    clustLikelihoods = map clustLik clusts ++ [emptyClustLik]
    clustLik clust = V.product $ V.zipWith3 featureLik distv (vectranspose ndim clust) newstim
    featureLik dist sample query = (fst . dist) (catMaybeV sample) $ query
    emptyClustLik = V.product $ V.zipWith (\dist query -> (fst . dist) V.empty query) distv newstim
    clusts = clusterItems assignments stimuli
    stimlength = (V.length . V.head) stimuli
    distv = V.fromList distributions
    ndim = length distributions

-- | For all missing items in the querystim, infers the most likely possible
-- | outcome, taking all possible cluster assignments into account.
infer :: (ClusterPrior, [PDFFromSample]) -> Stims -> Partition
         -> Stim       -- ^ Stimulus being queried.
         -> ([Double], [Double]) -- ^ (List of predictions for the dimension features, and List of posteriors for the clusters, INCLUDING EMPTY)
infer (cprior, distributions) stimuli assignments querystim = (prediction, post)
  where
    prediction = map inferDim queryDims
    inferDim i = sum $ (zipWith (*) post) (clustPreds!!i)
    clustPreds = transpose $ map (clusterPredictions distributions) (clusters++[V.empty])
    clusters = clusterItems assignments stimuli
    queryDims = V.toList $ V.map snd $ V.filter (isNothing . fst) $ V.zip querystim (V.enumFromN 0 (V.length querystim))
    post = clusterPosterior cprior distributions stimuli assignments querystim

-- | Cluster centroids for all clusters in all dimensions
summarizeClusters :: [PDFFromSample] -> Stims -> Partition -> [[Double]]
summarizeClusters pdfs stims partition = map summfun clusts
  where
    summfun = clusterPredictions pdfs
    clusts = clusterItems partition stims

