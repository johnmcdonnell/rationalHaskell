{-# LANGUAGE TupleSections #-}

module JohnFuns 
    (
    debug,
    debugMessage,
    countUnique,
    safetail,
    gatherBy,
    gather,
    vectranspose
    ) where

import qualified Data.List as List
import qualified Data.Map as Map
import qualified Data.Vector as V
import Data.Vector (Vector, (!))
import Debug.Trace
import Data.Function (on)
import Data.Typeable (Typeable, cast)
import Data.Maybe (fromMaybe)

-- * Unsafe IO for 

-- | Unsafe IO, just spits out the argument.
debug v = trace (show v) v

-- | Unsafe IO, with an added message.
debugMessage :: (Show a) => String -> a -> a
debugMessage message v = trace (message ++ show v) v

    
-- * Things that should be in Data.List

-- | Tail which returns an empty list if empty
safetail :: [a] -> [a]
safetail [] = []
safetail x = tail x

-- | Returns a list of items, grouped by the given function.
gatherBy :: Ord b => (a -> b) -> [a] -> [[a]]
gatherBy f = (List.groupBy ((==) `on` f)) . (List.sortBy (compare `on` f))

-- | Returns a list of items, grouped by identity.
gather :: (Ord a, Ord b) => (a -> b) -> [a] -> [[a]]
gather f = gatherBy id

-- | Returns map from each value to how often it occured.
countUnique :: Ord a => [a] -> Map.Map a Int
countUnique = List.foldl' additem Map.empty
  where
    additem = flip $ flip (Map.insertWith (+)) 1

-- * Things that should be in Data.Vector
vectranspose :: Int -> Vector (Vector a) -> Vector (Vector a)
vectranspose mindim vec 
  | V.null vec       = V.replicate mindim V.empty
  | V.length vec == 1 = V.map V.singleton (vec!0)
  | otherwise        = V.map (\c -> V.map (\r -> vec!r!c ) (V.enumFromN 0 rows)) (V.enumFromN 0 cols)
  where
    cols = (V.length . V.head) vec
    rows = V.length vec


