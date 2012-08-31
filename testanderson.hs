
-- import qualified Control.Monad.State as State
import Data.List (sortBy, transpose)
import Data.Function (on)
import Data.Maybe
import Control.Applicative
import Control.Monad
import Control.Monad.Random
import System.Random
import System.Random.Shuffle

import Math.Statistics

import Graphics.Gnuplot.Simple

import JohnFuns (debug, groupsortBy)
import Random

import Stats
import Rational
import Anderson

-- Convenience functions

plotpoints :: [(Double, Double)] -> IO ()
plotpoints = plotListStyle [] defaultStyle{ plotType=Dots }



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
    samples <- forM params (\(mu, sigma) -> evalRandIO $ (take n) . (map Just) <$> (normalsM mu sigma))
    let stims = map (\x -> [x]) $ concat samples
    items <- evalRandIO $ shuffleM stims
    let tpriors = [tPosterior (mean itemset, stddev itemset, 1, 1) | itemset <- (transpose . (map catMaybes)) items ]
    return (items, tpriors)


twodtask :: [(Double, Double)] -> Int -> IO ([Stim], [PDFFromSample])
twodtask params n = do
    let mergeDims = (\(x,y) -> zipWith (\x y -> [x,y]) x y)
    stims <- forM params (\(mu, sigma) -> evalRandIO $ (mergeDims . (splitAt n) . (take $ n*2) . map Just) <$> (normalsM mu sigma))
    items <- evalRandIO $ shuffleM (concat stims)
    let tpriors = [tPosterior (mean itemset, stddev itemset, 1, 1) | itemset <- (transpose . map catMaybes) items ]
    return (items, tpriors)

    

zeithamovaMaddox :: (Double, Double) -> Int -> IO ([Stim], [PDFFromSample])
zeithamovaMaddox (contalpha, contlambda) n = do
    -- length bimodal
    let bimodmean1 = 187.5
    let bimodmean2 = 412.5
    let bimodsd = 12.5
    let unimodmean = 45
    let unimodsd = 15
    astims <-  evalRandIO $ map (map Just) <$> binormals bimodmean1 unimodmean bimodsd unimodsd n
    bstims <- evalRandIO $ map (map Just) <$> binormals bimodmean2 unimodmean bimodsd unimodsd n
    items <- evalRandIO $ shuffleM (astims ++ bstims)
    let itemswithlabels = map (\x -> x ++ [Nothing]) items
    let tpriors = [tPosterior (mean itemset, stddev itemset, contalpha, contlambda) | itemset <- (transpose . (map catMaybes)) items ]
    let multinomprior =  multinomialPosterior [1, 1]
    return (itemswithlabels, tpriors ++ [multinomprior])

tvLabelFun :: Stim -> Stim
tvLabelFun [ Just bimod, Just unimod ] 
  | bimod < 300 && unimod > 55 = [ Just bimod, Just unimod, Just 0  ]
  | bimod > 300 && unimod < 35 = [ Just bimod, Just unimod, Just 1  ]
  | otherwise                  = [ Just bimod, Just unimod, Nothing ]

tvTask :: (Double, Double) -> Int -> IO ([Stim], [PDFFromSample])
tvTask (contalpha, contlambda) n = do
    -- length bimodal
    let bimodmean1 = 187.5
    let bimodmean2 = 412.5
    let bimodsd = 12.5
    let unimodmean = 45
    let unimodsd = 15
    astims <- evalRandIO $  (map (map Just)) <$> binormals bimodmean1 unimodmean bimodsd unimodsd n
    bstims <- evalRandIO $ (map (map Just)) <$> binormals bimodmean2 unimodmean bimodsd unimodsd n
    items <- evalRandIO $ shuffleM (astims ++ bstims)
    let itemswithlabels = map tvLabelFun items
    let tpriors = [tPosterior (mean itemset, stddev itemset, contalpha, contlambda) | itemset <- (transpose . (map catMaybes)) items ]
    let multinomprior =  multinomialPosterior [1, 1]
    return (itemswithlabels, tpriors ++ [multinomprior])

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

testTVTask = do
    -- Set up priors
    (task, distpriors) <- tvTask (1, 1) 100
    print task
    print $ sortBy (compare `on` fst) $ map (\((Just bimod):(Just unimod):xs) -> (bimod, unimod)) task
    -- plotpoints $ map (\((Just bimod):(Just unimod):xs) -> (bimod, unimod)) task
    let prior  = (1.0, distpriors)
    
    -- Now run the model
    partition <- andersonSample prior task

    -- Dump all those mabies
    let demabify [bimod,unimod,lab] = [fromJust bimod, fromJust unimod, fromMaybe 1 lab]
    let demabified = map demabify task
    let results = map (\(x,Just y) -> (x,y)) $ (filter (isJust . snd)) $ zip demabified partition
    
    -- Print it out
    --forM_ results print

    -- Plot it
    -- let grouped = map ((take 2) . fst) $ groupsortBy snd results
    let grouped = map (map (\([x,y,_],_) -> (x,y))) $ groupsortBy snd results
    -- print $ map (\x -> (Dots, x))
    -- mapM_ print $ sortBy (compare `on` fst) (grouped!!0)
    plotListsStyle [] $ map (\x -> (defaultStyle{ plotType=Dots }, x)) grouped

main = do
    -- testMedinSchaffer
    -- testContinuous
    -- testZeithamova
    testTVTask

