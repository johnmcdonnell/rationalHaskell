module Stats (PDFFromSample,
              ClusterPrior,
              bernoulliPosterior,
              tPosterior,
              dirichletProcess) where

import qualified Data.Vector as V
import qualified Data.Map as Map
import Data.List (group, sort)
import qualified Statistics.Distribution as StatDist
import Statistics.Distribution.StudentT as StatDist.T
import Statistics.Sample
-- import Math.Statistics

import Utils

-- helper funs
tatval :: Double -- ^ df
       -> Double -- ^ x
       -> Double -- ^ PDF evaluated at x
tatval = StatDist.density . studentT

sq x = x*x

-- * Generating PDFs from a sample of points

type PDFFromSample = V.Vector Double -> (Maybe Double -> Double, Double)

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
    variance = if n>0 then varianceUnbiased sample else 0 -- TODO: is varianceUnbiased correct here? Also is this the right way to deal with a small sample?
    xbar = if n>0 then mean sample else mu0
    lambdai = lambda0 + n
    alphai = alpha0 + n
    n = fromIntegral $ V.length sample
    (mu0, sigma0, alpha0, lambda0) = prior

bernoulliPosterior :: [Double] -> PDFFromSample
bernoulliPosterior alphas sample = (pdfFun, pdfFun (Just 1))
  where
    pdfFun Nothing = 1
    pdfFun (Just query) = (cj + alphas!!(floor query)) / (n + alpha0)
      where
        cj = (fromIntegral . V.length . (V.filter (==query))) sample
    ns = map (fromIntegral . snd) $ Map.toAscList $ countUnique $ V.toList $ sample
    n = fromIntegral $ V.length sample
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

