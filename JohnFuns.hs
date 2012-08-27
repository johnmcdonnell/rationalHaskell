
module JohnFuns 
    (
    (=.<),
    groupsortBy,
    debug,
    debugMessage
    ) where

import qualified Data.List as L
import Debug.Trace
import Data.Function (on)

-- | Unsafe IO, just spits out the argument.
debug v = traceShow v v

-- | Unsafe IO, with an added message.
debugMessage :: (Show a) => String -> a -> a
debugMessage message v = trace (message ++ show v) v

-- |first sort by the thing, then group by it. (this is maybe preferred behavior for group)
groupsortBy :: Ord b => (a -> b) -> [a] -> [[a]]
groupsortBy func = L.groupBy ((==) `on` func) . (L.sortBy (compare `on` func))

-- |infix version of fmap.
infixr 9 =.<
a =.< f = fmap f a

