
-- import qualified Control.Monad.State as State
import Data.List (transpose)
import Data.Maybe
import System.Random
import Control.Monad
import Control.Applicative

import Math.Statistics

import JohnFuns (debug)
import Random

import Stats
import Rational
import Anderson


medinSchafferTask :: [Double] -> ([Stim], [PDFFromSample])
medinSchafferTask multinomAlphas = (medinSchafferStims, andersondists)
  where 
    andersondists = replicate 5 multinom_prior
    multinom_prior = multinomialPosterior multinomAlphas
    medinSchafferStims = map (map Just) medinSchafferItems
    medinSchafferItems = [[1,1,1,1,1], 
                          [1,0,1,0,1], 
                          [1,0,1,1,0], 
                          [0,0,0,0,0], 
                          [0,1,0,1,1], 
                          [0,1,0,0,0]]
    

testMedinSchaffer = do
    let (task, dists) = medinSchafferTask [1,1]
    let couplingParam = 1.0
    andersonSample (couplingParam, dists) task

onedtask :: [(Double, Double)] -> Int -> IO ([Stim], [PDFFromSample])
onedtask params n = do
    samples <- forM params (\(mu, sigma) -> (take n) . (map Just) <$> runRandomIO (normals mu sigma))
    let stims = map (\x -> [x]) $ concat samples
    items <- runRandomIO $ shuffleNR stims (n * (length params))
    let tpriors = [tPosterior (mean itemset, stddev itemset, 1, 1) | itemset <- (transpose . (map catMaybes)) items ]
    return (items, tpriors)


twodtask :: [(Double, Double)] -> Int -> IO ([Stim], [PDFFromSample])
twodtask params n = do
    let mergeDims = (\(x,y) -> zipWith (\x y -> [x,y]) x y)
    stims <- forM params (\(mu, sigma) -> (mergeDims . (splitAt n) . (take $ n*2) . map Just) <$> runRandomIO (normals mu sigma))
    items <- runRandomIO $ shuffleNR (concat stims) (n * (length params))
    let tpriors = [tPosterior (mean itemset, stddev itemset, 1, 1) | itemset <- (transpose . map catMaybes) items ]
    return (items, tpriors)

binormalSamples :: Double -> Double -> Double -> Double -> Int -> IO [Stim]
binormalSamples mean1 mean2 sd1 sd2 n = do
    first <- map Just <$> runRandomIO (normals mean1 sd1)
    second <- map Just <$> runRandomIO (normals mean2 sd2)
    return $ take n $ zipWith (\x y -> [x,y]) first second
    

zeithamovaMaddox :: (Double, Double) -> Int -> IO ([Stim], [PDFFromSample])
zeithamovaMaddox (contalpha, contlambda) n = do
    -- length bimodal
    let bimodmean1 = 187.5
    let bimodmean2 = 412.5
    let bimodsd = 12.5
    let unimodmean = 45
    let unimodsd = 15
    astims <- binormalSamples bimodmean1 unimodmean bimodsd unimodsd n
    bstims <- binormalSamples bimodmean2 unimodmean bimodsd unimodsd n
    items <- runRandomIO $ shuffleNR (astims ++ bstims) (n * 2)
    let tpriors = [tPosterior (mean itemset, stddev itemset, contalpha, contlambda) | itemset <- (transpose . (map catMaybes)) items ]
    return (items, tpriors)

testContinuous = do
    let mu1 = 1
    let sigma1 = 1
    let mu2 = -1
    let sigma2 = 1
    let n = 3
    putStrLn "Testing on:" 
    putStrLn $ "mu1=" ++ (show mu1) ++ " sigma1=" ++ (show sigma1) ++ " n=" ++ (show n)
    putStrLn $ "mu2=" ++ (show mu2) ++ " sigma2=" ++ (show sigma2) ++ " n=" ++ (show n)
    (task, distpriors) <- twodtask [(mu1, sigma1), (mu2, sigma2)] n
    
    let prior  = (1.0, distpriors)
    partition <- andersonSample prior task
    forM_ (zip task (catMaybes partition)) print

testZeithamova = do
    (task, distpriors) <- zeithamovaMaddox (1, 1) 100
    
    let prior  = (1.0, distpriors)
    partition <- andersonSample prior task
    forM_ (zip task (catMaybes partition)) print

main = do
    -- testMedinSchaffer
    -- testContinuous
    testZeithamova

