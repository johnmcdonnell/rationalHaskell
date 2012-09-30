{-# LANGUAGE Rank2Types #-}

module Rational (
                 clusterPosterior,
                 clusterItems,
                 infer,
                 Stim,
                 Stims,
                 Partition
                ) where

import Data.List (nub, group, groupBy, transpose, sortBy, splitAt)
import Data.Function (on)
import Data.Maybe

import Control.Monad
import Control.Monad.ST
import Data.Vector (freeze, thaw, MVector, Vector, (!))
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM


import JohnFuns
import Stats

-- * Convenience functions
safetail :: [a] -> [a]
safetail [] = []
safetail x = tail x

gatherBy :: Ord b => (a -> b) -> [a] -> [[a]]
gatherBy f = (groupBy ((==) `on` f)) . (sortBy (compare `on` f))

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
    
vectranspose :: Int -> Vector (Vector a) -> Vector (Vector a)
vectranspose mindim vec 
  | V.null vec       = V.replicate mindim V.empty
  | V.length vec == 1 = V.map V.singleton (vec!0)
  | otherwise        = V.map (\c -> V.map (\r -> vec!r!c ) (V.enumFromN 0 rows)) (V.enumFromN 0 cols)
  where
    cols = (V.length . V.head) vec
    rows = V.length vec


-- Likelihood that a stimulus belongs to each cluster
clusterPosterior :: ClusterPrior -> [PDFFromSample] -> Stims -> Partition -> Stim -> [Double]
clusterPosterior cprior distributions stimuli assignments newstim = map (/norm) posterior
  where
    norm = sum posterior
    posterior = zipWith (*) clustPriors clustLikelihoods
    -- NOTE: watch out that the stim to be assigned isn't assigned in assignments.
    clustPriors = cprior $ (catMaybes . V.toList) assignments
    clustLikelihoods = map clustLik clusts ++ [emptyClustLik]
    clustLik clust = V.product $ V.zipWith3 (\dist sample query -> (fst . dist) ((catMaybes . V.toList) sample) query ) distv (vectranspose ndim clust) newstim
    emptyClustLik = V.product $ V.zipWith (\dist query -> (fst . dist) [] query) distv newstim
    clusts = clusterItems assignments stimuli
    stimlength = (V.length . V.head) stimuli
    distv = V.fromList distributions
    ndim = length distributions

-- | For all missing items in the querystim, infers the most likely possible
-- | outcome, taking all possible cluster assignments into account.
infer :: ClusterPrior -> [PDFFromSample] -> Stims -> Partition -> Stim -> [Double]
infer cprior distributions stimuli assignments querystim = map inferDim querydims
  where
    inferDim i = V.sum $ V.zipWith (*) post (clustPreds!i)
    clustPreds = vectranspose ndim $ V.map (V.fromList . (zipWith (\d c -> (snd . d . catMaybes . V.toList) c) distributions) . V.toList . (vectranspose ndim)) clusters
    clusters = V.fromList $ clusterItems assignments stimuli
    querydims = V.toList $ V.map snd $ V.filter (isNothing . fst) $ V.zip querystim (V.enumFromN 0 (V.length querystim))
    post = V.fromList $ clusterPosterior cprior distributions stimuli assignments querystim
    ndim = length distributions

