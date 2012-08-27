
module Rational (
                 dropByIndex,
                 replaceAtIndex,
                 clusterPosterior,
                 Stim,
                 Partition
                ) where

import Data.List (nub, group, groupBy, transpose, sortBy, splitAt)
import Data.Function (on)
import Data.Maybe (fromJust, catMaybes, maybe)

import JohnFuns (debug)
import Stats

-- * Convenience functions
safetail :: [a] -> [a]
safetail [] = []
safetail x = tail x

dropByIndex l i newval = take i l ++ drop (i+1) l
replaceAtIndex i newval l = (\(x, y) -> x ++ (newval:safetail y) ) $ splitAt i l

gatherBy :: Ord b => (a -> b) -> [a] -> [[a]]
gatherBy f = (groupBy ((==) `on` f)) . (sortBy (compare `on` f))

-- * Stim dattype and associated functions
type Stim = [Double]

-- * Partition dattype and associated functions
type Partition = [Maybe Int]

-- Safely remove item from Partition 
dropAssignment :: Partition -> Int -> Partition
dropAssignment assignments i
  | clustmembers /= [] = assignDropped
  | nclusts == clust   = assignDropped
  | otherwise          = map iterfun assignDropped
  where
    iterfun x     = if maybe False (<clust) x then x else fmap (+1) x
    nclusts       = (length . nub) assignRest
    clustmembers  = filter (==clust) assignRest
    assignRest    = catMaybes assignDropped
    assignDropped = replaceAtIndex i Nothing assignments
    clust         = fromJust $ assignments!!i

-- Likelihood that a stimulus belongs to each cluster
clusterPosterior :: ClusterPrior -> [PDFFromSample] -> [Stim] -> Partition -> Stim -> [Double]
clusterPosterior cprior distributions stimuli assignments newstim = map (/norm) posterior
  where
    norm = sum posterior
    posterior = zipWith (*) clustPriors clustLikelihoods
    -- NOTE: watch out that the stim to be assigned isn't assigned in assignments.
    clustPriors = cprior (catMaybes assignments)
    clustLikelihoods = debug $ map clustLik clusts ++ [emptyClustLik]
    clustLik clust = product $ zipWith3 (\dist sample query -> dist sample query ) distributions (transpose clust) newstim
    emptyClustLik = product $ zipWith3 (\dist sample query -> dist sample query ) distributions (replicate stimlength []) newstim
    clusts = (map (map snd) $ gatherBy fst $ (zip assignments stimuli))
    stimlength = (length . head) stimuli


