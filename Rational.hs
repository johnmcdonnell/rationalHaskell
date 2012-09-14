
module Rational (
                 dropByIndex,
                 replaceAtIndex,
                 clusterPosterior,
                 clusterItems,
                 infer,
                 Stim,
                 Partition
                ) where

import Data.List (nub, group, groupBy, transpose, sortBy, splitAt)
import Data.Function (on)
import Data.Maybe

import JohnFuns
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
type Stim = [Maybe Double]

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

clusterItems :: Partition -> [Stim] -> [[Stim]]
clusterItems = curry (takeStims . groupByClust . takeIncluded . (uncurry zip))
  where
    takeStims = map (map snd)
    groupByClust = groupsortBy fst
    takeIncluded = filter (isJust . fst)

-- Likelihood that a stimulus belongs to each cluster
clusterPosterior :: ClusterPrior -> [PDFFromSample] -> [Stim] -> Partition -> Stim -> [Double]
clusterPosterior cprior distributions stimuli assignments newstim = map (/norm) posterior
  where
    norm = sum posterior
    posterior = zipWith (*) clustPriors clustLikelihoods
    -- NOTE: watch out that the stim to be assigned isn't assigned in assignments.
    clustPriors = cprior (catMaybes assignments)
    clustLikelihoods = map clustLik clusts ++ [emptyClustLik]
    clustLik clust = product $ zipWith3 (\dist sample query -> (fst . dist) (catMaybes sample) query ) distributions (transpose clust) newstim
    emptyClustLik = product $ zipWith (\dist query -> (fst . dist) [] query) distributions newstim
    clusts = clusterItems assignments stimuli
    stimlength = (length . head) stimuli

-- | For all missing items in the querystim, infers the most likely possible
-- | outcome, taking all possible cluster assignments into account.
infer :: ClusterPrior -> [PDFFromSample] -> [Stim] -> Partition -> Stim -> [Double]
infer cprior distributions stimuli assignments querystim = map inferDim querydims
  where
    inferDim i = sum $ zipWith (*) post (clustPreds!!i)
    clustPreds = debug $ transpose $ map ((zipWith (\d c -> (snd . d . catMaybes) c) distributions) . transpose) clusters
    clusters = clusterItems assignments stimuli
    querydims = map snd $ filter (isNothing . fst) $ zip querystim [0..]
    post = clusterPosterior cprior distributions stimuli assignments querystim

