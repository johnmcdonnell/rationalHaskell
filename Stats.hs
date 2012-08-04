

module Stats (PDFFromSample,
              multinomialDistribtution,
              dirichletProcess) where

import Data.List (group, sort)

type PDFFromSample = [Double] -> Double -> Double
fillerDistribtution :: PDFFromSample
fillerDistribtution sample query = 0.5

multinomialDistribtution :: Double -> [Double] -> PDFFromSample
multinomialDistribtution alpha0 alphas sample query = (cj +alphas!!(floor query)) / (n + alpha0)
  where
    cj = fromIntegral $ length $ filter (==query) sample
    n = fromIntegral $ length sample

-- Takes a list of assignments (assumed to be consecutive ints) and returns the
-- likelihood for each cluster, including a new cluster.
dirichletProcess :: [Int] -> Double -> [Double]
dirichletProcess [] _ = [1]
dirichletProcess assignments alpha = pPartitions ++ [pNew]
  where
    pNew = if isInfinite alpha then 1 else alpha / (n+alpha)
    pPartitions = map (\x -> (fromIntegral x) / (n+alpha)) counts
    counts = map length $ (group . sort) assignments
    n = fromIntegral $ length assignments
