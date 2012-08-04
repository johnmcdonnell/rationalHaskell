
import Data.List (nub, group, groupBy, transpose, sort, splitAt)
import Data.Function (on)
import Data.Maybe (fromJust, catMaybes, maybe)
import Stats

-- * Convenience functions
changeByIndex l i newval = take i l ++ drop (i+1) l
replaceAtIndex l i newval = (\(x, y) -> x ++ (newval:tail y) ) $ splitAt i l

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
    assignDropped = replaceAtIndex assignments i Nothing
    clust         = fromJust $ assignments!!i


-- Likelihood that a stimulus belongs to each cluster
assignDist :: [PDFFromSample] -> [Stim] -> Partition -> Stim -> [Double]
assignDist distributions stimuli assignments newstim = clustLikelihoods
  where
    clustLikelihoods = map clustLik (clusts ++ [[]] )
    clustLik clust = product $ zipWith3 (\dist sample query -> dist sample query ) distributions (transpose clust) newstim
    clusts = map (map snd) $ groupBy ((==) `on` fst) $ sortBy (compare `on` fst) $ (zip assignments stimuli)


