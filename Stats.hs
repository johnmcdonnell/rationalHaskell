

module Stats (PDFFromSample,
              ClusterPrior,
              binomialPosterior,
              tPosterior,
              dirichletProcess) where

import Data.List (group, sort)
import qualified Statistics.Distribution as StatDist
import Statistics.Distribution.StudentT as StatDist.T
import Math.Statistics

import JohnFuns

-- helper funs
tatval :: Double -- ^ df
       -> Double -- ^ x
       -> Double -- ^ PDF evaluated at x
tatval = StatDist.density . studentT

-- * Generating PDFs from a sample of points

type PDFFromSample = [Double] -> (Maybe Double -> Double, Double)

fillerDistribtution :: PDFFromSample
fillerDistribtution sample = ((\x -> 0.5), 0.5)

type ClusterPrior = [Int] -> [Double]

tPosterior :: (Double, Double, Double, Double) -> PDFFromSample
tPosterior prior sample  = (pdfFun, mui)
  where
    pdfFun (Just query) = tatval alphai ((query - mui) / tstdev)
    pdfFun Nothing = 1
    tstdev = sigmai * (sqrt (1 + (1/lambdai)))
    sigmai = sqrt $ (alpha0*(sq sigma0) + 
                      (n-1)*variance +
                      (lambda0*n/(lambda0+n))*(sq (mu0-xbar))
                    ) / (alpha0 + n)
    mui = (lambda0*mu0 + n * xbar) / (lambda0 + n)
    variance = if n>0 then pvar sample else 0 -- TODO: is pvar correct here? Also is this the right way to deal with a small sample?
    xbar = if n>0 then mean sample else mu0
    lambdai = lambda0 + n
    alphai = alpha0 + n
    n = fromIntegral $ length sample
    (mu0, sigma0, alpha0, lambda0) = prior

-- Binomial values must be in the set [0..n-1] where n is the number of
-- possible values.
binomialPosterior :: [Double] -> PDFFromSample
binomialPosterior alphas sample = (pdfFun, pdfFun (Just 1))
  where
    pdfFun Nothing = 1
    pdfFun (Just query) = (cj + alphas!!(floor query)) / (n + alpha0)
      where
        cj = (fromIntegral . length . (filter (==query))) sample
    n = sum ns
    ns = map (fromIntegral . length) $ groupsortBy id sample
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


