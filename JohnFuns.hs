
module JohnFuns 
    (
    groupsortBy,
    debug,
    debugMessage
    ) where

import qualified Data.List as L
import Debug.HTrace
import Data.Function (on)
import Data.Typeable (Typeable, cast)
import Data.Maybe (fromMaybe)

-- | Unsafe IO, just spits out the argument.
debug v = htrace (show v) v

-- | Unsafe IO, with an added message.
debugMessage :: (Show a) => String -> a -> a
debugMessage message v = htrace (message ++ show v) v

-- |first sort by the thing, then group by it. (this is maybe preferred behavior for group)
groupsortBy :: Ord b => (a -> b) -> [a] -> [[a]]
groupsortBy func = L.groupBy ((==) `on` func) . (L.sortBy (compare `on` func))

-- | Casting 

-- toString :: Show a => a -> String
-- toString
