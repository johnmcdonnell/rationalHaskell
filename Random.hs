module Random (
               normalsM,
               binormals,
               ) where

import Control.Monad
import Control.Monad.Random
import System.Random (Random, StdGen, getStdGen, mkStdGen, random)
import System.Random.Shuffle
import Control.Applicative ((<$>))
import Statistics.Distribution
import Statistics.Distribution.Normal
import qualified Numeric.LinearAlgebra as LA

import Utils

normalsM :: RandomGen g => Double -> Double -> Rand g [Double]
normalsM mu sd = do
    let dist = normalDistr mu sd
    randoms <- getRandomRs (0, 1)
    return $ map (quantile dist) randoms

bivariatecovmat :: (Double, Double) -> Double -> LA.Matrix Double
bivariatecovmat (sigma1, sigma2) rho = LA.fromLists [[sigma1*sigma1, rho*sigma1*sigma2], [rho*sigma1*sigma2, sigma2*sigma2]]

binormals :: (RandomGen g) => LA.Vector Double -> (Double, Double) -> Double -> Int -> Rand g (LA.Matrix Double)
binormals means sigmas rho n = do
    seed <- getRandom
    return $ LA.gaussianSample seed n means (bivariatecovmat sigmas rho)

