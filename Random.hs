module Random (lognormalM,
               normalM,
               normalsM,
               binormals
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

lognormalM :: RandomGen g => Double -> Double -> Rand g Double
lognormalM mu sd = liftM exp (normalM mu sd)

normalM :: RandomGen g => Double -> Double -> Rand g Double
normalM mu sd = do
    let dist = normalDistr mu sd
    random <- getRandomR (0, 1)
    return $ quantile dist random

normalsM :: RandomGen g => Double -> Double -> Rand g [Double]
normalsM mu sd = (sequence . repeat) $ normalM mu sd

bivariatecovmat :: (Double, Double) -> Double -> LA.Matrix Double
bivariatecovmat (sigma1, sigma2) rho = LA.fromLists [[sigma1*sigma1, rho*sigma1*sigma2], [rho*sigma1*sigma2, sigma2*sigma2]]

binormals :: (RandomGen g) => LA.Vector Double -> (Double, Double) -> Double -> Int -> Rand g (LA.Matrix Double)
binormals means sigmas rho n = do
    seed <- getRandom
    return $ LA.gaussianSample seed n means (bivariatecovmat sigmas rho)

