
module Random (
               normalsM,
               binormals
               ) where

import Control.Monad
import Control.Monad.Random
import System.Random (Random, StdGen, getStdGen, mkStdGen, random)
import System.Random.Shuffle
import Data.Random.Normal
import Control.Applicative ((<$>))

import JohnFuns

-- Thanks to: http://stackoverflow.com/questions/2110535/sampling-sequences-of-random-numbers-in-haskell/2125329#2125329
-- for a start on this.

normalsM :: RandomGen g => Double -> Double -> Rand g [Double]
normalsM mu sd = do
    s <- getSplit
    return $ normals' (mu, sd) s

binormals :: (RandomGen g) => Double -> Double -> Double -> Double -> Int -> Rand g [[Double]]
binormals mean1 mean2 sd1 sd2 n = do
    first <- (take n) <$> (normalsM mean1 sd1)
    second <- (take n) <$> (normalsM mean2 sd2)
    return $ zipWith (\x y -> [x,y]) first second


