{-# LANGUAGE TupleSections #-}

module JohnFuns 
    (
    groupsortBy,
    debug,
    debugMessage,
    countUnique
    ) where

import qualified Data.List as List
import qualified Data.Map as Map
import Debug.Trace
import Data.Function (on)
import Data.Typeable (Typeable, cast)
import Data.Maybe (fromMaybe)

-- | Unsafe IO, just spits out the argument.
debug v = trace (show v) v

-- | Unsafe IO, with an added message.
debugMessage :: (Show a) => String -> a -> a
debugMessage message v = trace (message ++ show v) v

-- |first sort by the thing, then group by it. (this is maybe preferred behavior for group)
groupsortBy :: Ord b => (a -> b) -> [a] -> [[a]]
groupsortBy func = List.groupBy ((==) `on` func) . (List.sortBy (compare `on` func))


-- | Returns map from each value to a 
countUnique :: Ord a => [a] -> Map.Map a Int
countUnique = List.foldl' additem Map.empty
  where
    additem = flip $ flip (Map.insertWith (+)) 1
    


