{-# LANGUAGE Rank2Types #-}

module Rational (
                 clusterPosterior,
                 clusterItems,
                 infer,
                 Stim,
                 Stims,
                 Partition,
                 clusterPredictions,
                 summarizeClusters
                ) where

import Data.List (nub, group, groupBy, transpose, sortBy, splitAt)
import Data.Function (on)
import Data.Maybe

import Control.Monad
import Control.Monad.ST
import Data.Vector (freeze, thaw, MVector, Vector, (!))
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM

import Statistics.Sample

import JohnFuns
import Stats

-- * Convenience functions
-- * Stim dattype and associated functions
type Stim = Vector (Maybe Double)

-- * Stim dattype and associated functions
type Stims = Vector Stim

-- * Partition dattype and associated functions
type Partition = Vector (Maybe Int)

-- Safely remove item from Partition 
dropAssignment :: Partition -> Int -> Partition
dropAssignment assignments i = runST $ do
    let clustm = assignments!i
    let clust = fromJust clustm
    let nInClust = V.length $ V.filter (==clustm) assignments
    let iterfun x = if maybe True (<clust) x then x else fmap (\x->x-1) x :: Maybe Int
    let n = V.length assignments
    thawedassign <- V.unsafeThaw assignments
    VM.unsafeWrite thawedassign i Nothing
    if (nInClust > 1) then return ()
                      else forM_ [0..n-1] (\j ->
                          (VM.unsafeRead thawedassign j )>>=
                          (VM.unsafeWrite thawedassign j) . iterfun)
    V.unsafeFreeze thawedassign 

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
infer :: ClusterPrior -> [PDFFromSample] -> Stims -> Partition -> Stim -> [Double]
infer cprior distributions stimuli assignments querystim = map inferDim querydims
  where
    inferDim i = sum $ (zipWith (*) post) (clustPreds!!i)
    clustPreds = transpose $ map (clusterPredictions distributions) clusters
    clusters = clusterItems assignments stimuli
    querydims = V.toList $ V.map snd $ V.filter (isNothing . fst) $ V.zip querystim (V.enumFromN 0 (V.length querystim))
    post = clusterPosterior cprior distributions stimuli assignments querystim


-- | Cluster centroids for all clusters in all dimensions
summarizeClusters :: [PDFFromSample] -> Stims -> Partition -> [[Double]]
summarizeClusters pdfs stims partition = map summfun clusts
  where
    summfun = clusterPredictions pdfs
    clusts = clusterItems partition stims

