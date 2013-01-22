module Random (
               normalsM,
               binormals,
               ) where

import Control.Monad
import Control.Monad.Random
import System.Random (Random, StdGen, getStdGen, mkStdGen, random)
import System.Random.Shuffle
import Data.Random.Normal
import Control.Applicative ((<$>))
import qualified Numeric.LinearAlgebra as LA

import Utils

-- Thanks to: http://stackoverflow.com/questions/2110535/sampling-sequences-of-random-numbers-in-haskell/2125329#2125329
-- for a start on this.

normalsM :: RandomGen g => Double -> Double -> Rand g [Double]
normalsM mu sd = do
    s <- getSplit
    return $ normals' (mu, sd) s

-- binormalscov :: (RandomGen g) => Double -> Double -> Double -> Double -> Double -> Int -> Rand g [[Double]]
-- binormalscov mean1 mean2 sd1 sd2 cov n = do
--     s <- getSplit
--     let samples = (take (n*2)) $ (normals s)
--     let sampleMat = (n LA.>< 2) samples
--     let covmat = LA.fromLists [[sd1, cov], [cov, sd2]]
--     let rotated = sampleMat LA.<> covmat
--     return $ map (zipWith (+) [mean1, mean2]) (LA.toLists rotated)

bivariatecovmat :: (Double, Double) -> Double -> LA.Matrix Double
bivariatecovmat (sigma1, sigma2) rho = LA.fromLists [[sigma1*sigma1, rho*sigma1*sigma2], [rho*sigma1*sigma2, sigma2*sigma2]]

binormals :: (RandomGen g) => LA.Vector Double -> (Double, Double) -> Double -> Int -> Rand g (LA.Matrix Double)
binormals means sigmas rho n = do
    seed <- getRandom
    return $ LA.gaussianSample seed n means (bivariatecovmat sigmas rho)

