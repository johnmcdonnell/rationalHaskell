

module Stats (PDFFromSample,
              ClusterPrior,
              multinomialPosterior,
              tPosterior,
              dirichletProcess) where

import Data.List (group, sort)
import qualified Statistics.Distribution as StatDist
import Statistics.Distribution.StudentT as StatDist.T
import Math.Statistics

import JohnFuns (debug)

-- helper funs
tatval :: Double -- ^ df
       -> Double -- ^ x
       -> Double -- ^ PDF evaluated at x
tatval = StatDist.density . studentT

-- * Generating PDFs from a sample of points

type PDFFromSample = [Double] -> Double -> Double
fillerDistribtution :: PDFFromSample
fillerDistribtution sample query = 0.5

type ClusterPrior = [Int] -> [Double]

tPosterior :: (Double, Double, Double, Double) -> PDFFromSample
tPosterior prior sample query = debug $ tatval alphai x
  where
    x = (query - mui) / tstdev
    tstdev = sigmai * (sqrt (1 + (1/lambdai)))
    sigmai = sqrt $ (debug $ alpha0*(sq sigma0) + 
                      (n-1)*var +
                      (lambda0*n/(lambda0+n))*(sq (mu0-xbar))
                    ) / (alpha0 + n)
    mui = (lambda0*mu0 + n * xbar) / (lambda0 + n)
    var = pvar sample
    xbar = mean sample
    lambdai = lambda0 + n
    alphai = alpha0 + n
    n = fromIntegral $ length sample
    (mu0, sigma0, alpha0, lambda0) = prior

multinomialPosterior :: [Double] -> PDFFromSample
multinomialPosterior alphas sample query = (cj + (debug alphas)!!(debug $ floor query)) / (n + alpha0)
  where
    cj = fromIntegral $ length $ filter (==query) sample
    n = fromIntegral $ length sample
    alpha0 = sum alphas

-- Takes a list of assignments (assumed to be consecutive ints) and returns the
-- likelihood for each cluster, including a new cluster.
dirichletProcess :: Double -> ClusterPrior
dirichletProcess _ [] = [1]
dirichletProcess alpha assignments = pPartitions ++ [pNew]
  where
    pNew = if isInfinite alpha then 1 else alpha / (n+alpha)
    pPartitions = map (\x -> (fromIntegral x) / (n+alpha)) counts
    counts = map length $ (group . sort) assignments
    n = fromIntegral $ length assignments


