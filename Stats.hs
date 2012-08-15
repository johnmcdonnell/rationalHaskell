

module Stats (PDFFromSample,
              ClusterPrior,
              multinomialDistribution,
              dirichletProcess) where

import Data.List (group, sort)

import JohnFuns (debug)

type PDFFromSample = [Double] -> Double -> Double
fillerDistribtution :: PDFFromSample
fillerDistribtution sample query = 0.5

type ClusterPrior = [Int] -> [Double]

multinomialDistribution :: [Double] -> PDFFromSample
multinomialDistribution alphas sample query = (cj + alphas!!(floor query)) / (n + alpha0)
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


