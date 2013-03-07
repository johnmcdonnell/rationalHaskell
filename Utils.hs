{-# LANGUAGE TupleSections #-}

module Utils
    (
    debug,
    debugMessage,
    both,
    countUnique,
    maxWithArg,
    argmax,
    safetail,
    gatherBy,
    gather,
    vectranspose,
    catMaybeV,
    inplacewrite,
    Down (Down),
    sq
    ) where

import qualified Data.List as List
import qualified Data.Map as Map
import qualified Data.Vector as V
import Data.Vector (Vector, (!))
import qualified Data.Foldable as Fold
import qualified Data.Vector.Mutable as VM
import Debug.Trace
import Data.Function (on)
import Data.Typeable (Typeable, cast)
import Data.Maybe (fromMaybe, fromJust, isJust)
import Control.Arrow
-- import Control.Monad.ST

-- * Unsafe IO for 

-- | Unsafe IO, just spits out the argument.
debug v = trace (show v) v

-- | Unsafe IO, with an added message.
debugMessage :: (Show a) => String -> a -> a
debugMessage message v = trace (message ++ show v) v

-- * Fun with arrows

both :: Arrow a => a b c -> a (b, b) (c, c)
both f = f *** f

    
-- * Things that should be in Data.List

-- | Returns (argmax, max)
maxWithArg :: (Ord a) => [a] -> (Int, a)
maxWithArg items = List.foldl1' compare loopvals
  where
    compare opt next = if (snd next) > (snd opt) then next else opt
    loopvals = zip [0..] items

argmax :: Ord a => [a] -> Int
argmax = fst . maxWithArg 

-- | Tail which returns an empty list if empty
safetail :: [a] -> [a]
safetail [] = []
safetail x = tail x

-- | Returns a list of items, grouped by the given function.
gatherBy :: Ord b => (a -> b) -> [a] -> [[a]]
gatherBy f = (List.groupBy ((==) `on` f)) . (List.sortBy (compare `on` f))

-- | Returns a list of items, grouped by identity.
gather :: Ord a => [a] -> [[a]]
gather = gatherBy id

-- | Returns map from each value to how often it occured.
countUnique :: (Fold.Foldable f, Ord a) => f a -> Map.Map a Int
countUnique = Fold.foldl' additem Map.empty
  where
    additem sofar newkey = Map.insertWith (+) newkey 1 sofar

-- * Things that should be in Data.Vector
vectranspose :: Int -> Vector (Vector a) -> Vector (Vector a)
vectranspose mindim vec 
  | V.null vec       = V.replicate mindim V.empty
  | V.length vec == 1 = V.map V.singleton (vec!0)
  | otherwise        = V.map (\c -> V.map (\r -> vec!r!c ) (V.enumFromN 0 rows)) (V.enumFromN 0 cols)
  where
    cols = (V.length . V.head) vec
    rows = V.length vec

-- |Vector of the Just values in a vector
catMaybeV :: Vector (Maybe a) -> Vector a
catMaybeV = (V.map fromJust) . (V.filter isJust)


inplacewrite index newval = V.modify modfun
  where
    modfun v = VM.write v index newval

-- | For backwards compatibility 
newtype Down a = Down a deriving (Eq)

instance Ord a => Ord (Down a) where
    compare (Down x) (Down y) = y `compare` x

sq x = x*x
