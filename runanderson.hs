{-# LANGUAGE DeriveDataTypeable, RecordWildCards #-}

-- import qualified Control.Monad.State as State
import Data.List (sortBy)
import Data.Ord
import Data.Function (on)
import Data.Maybe
import qualified Data.Map as Map
import Data.Vector (freeze, thaw, MVector, Vector, (!))
import qualified Data.Vector as V
import Control.Applicative
import Control.Arrow
import Control.Monad.ST
import Control.Monad.Random
import System.Console.CmdArgs
import Text.CSV

import Utils

import Types
import Rational
import Anderson
import Tasks

-- Returns the likelihood of "0" as the first argument, and the most likely cluster as the second.
runTest :: (ClusterPrior, [PDFFromSample]) -> Stims -> Partition -> Stims -> Vector (Double, Int)
runTest prior stimuli assignments = V.map getpred
  where
    getpred = takeBestCluster . takePrediction . rawpred
    takeBestCluster = second (argmax . init)
    takePrediction = first head
    rawpred = infer prior stimuli assignments

-- * Run simulation on the various tasks

-- |Predictions for the Anderson model are given explicitly in Anderson (1991),
-- so they make a good test that the model is working as expected in the discrete case.
runMedinSchaffer :: IO ()
runMedinSchaffer = do
    let (task, dists) = medinSchafferTask [1,1]
    let couplingParam = dirichletProcess 1.0
    (partition, guesses) <- evalRandIO $ andersonSample EncodeActual (couplingParam, dists) task
    print partition

runContinuous = do
    let mu1 = 1
        sigma1 = 1
        mu2 = -1
        sigma2 = 1
        n = 3
    putStrLn "Testing on:" 
    putStrLn $ "mu1=" ++ (show mu1) ++ " sigma1=" ++ (show sigma1) ++ " n=" ++ (show n)
    putStrLn $ "mu2=" ++ (show mu2) ++ " sigma2=" ++ (show sigma2) ++ " n=" ++ (show n)
    (task, distpriors) <- evalRandIO $ twodtask [(mu1, sigma1), (mu2, sigma2)] n
    
    let prior  = (dirichletProcess 1.0, distpriors)
    (guesses, partition) <- evalRandIO $ andersonSample EncodeActual prior task
    V.forM_ (V.zip task (V.map (fromMaybe (-1)) partition)) print

runZeithamova = do
    (task, distpriors) <- evalRandIO $ zeithamovaMaddox (1, 1, Nothing, 0) 100
    
    let prior  = (dirichletProcess 1.0, distpriors)
    (guesses, partition) <- evalRandIO $ andersonSample EncodeActual prior task
    V.forM_ (V.zip task (V.map (fromMaybe (-1)) partition)) print

-- |These are the simulations for McDonnell et al. (2013)
runTVTask :: ModelArgs -> IO ()
runTVTask params@TVTask{..} = do
    let maxlab = 16
        (nlab, nounlab) = if nlabarg < 0 then (maxlab, True) else (nlabarg, False)
    
    -- Set up priors
    let filterfun = if nounlab then V.filter (isJust . V.last) else id
    (task, distpriors) <- evalRandIO $ first filterfun <$> mcdonnellTaskOrdered order (a0, lambda0, Just sigma0, bias) (28*4) nlab
    
    let prior  = (dirichletProcess alphaparam, distpriors)
    
    -- Now run the model
    (guesses, partition) <- evalRandIO $ andersonSample encoding prior task
    
    -- Dump all those mabies
    let demabify = V.map (fromMaybe (-9))
        demabified = V.map demabify task
        results = map (\(x,Just y) -> V.toList $ V.snoc x (fromIntegral y)) $ (filter (isJust . snd)) $ zip (V.toList demabified) (V.toList partition)
        resultstrings = map (map show) results
    
    -- Print it out
    putStrLn $ printCSV $ map ("STIM":) resultstrings
    
    putStrLn $ printCSV $ map (("CLUST":) . (map show)) $ summarizeClusters distpriors task partition
    
    -- If there were no labels, we have to assume they are associated with the largest clusters.
    let initialinferences = runTest prior task partition gridTest
        inferpartition = map snd $ V.toList initialinferences
        partitions = Map.toList $ countUnique inferpartition
        largestClusts = map fst $ take 2 $ sortBy (compare `on` (Utils.Down . snd)) partitions
        addlabel (i, clust) taskvec = inplacewrite index newstim taskvec
          where
            index = fromJust $ V.elemIndex (Just clust) partition 
            newstim = inplacewrite 2 (Just i) (taskvec!index)
        newtask = foldr addlabel task (zip [0..] largestClusts)
        modtask = if (nlab==0) then newtask else task
    
    -- Get inference at each point in a grid
    let testInferences = runTest prior modtask partition gridTest
        teststims = V.map (V.toList . (V.map (fromMaybe (-1))) . (V.take 2)) gridTest
        ret = V.zipWith (\stim (lab, clust) -> "INFER" : (map show (stim ++ [lab])) ++ [show clust]) teststims testInferences
    putStrLn $ printCSV $ V.toList ret

runVandistTask :: ModelArgs -> IO ()
runVandistTask params@Vandist{..} = do
    let priorparams = (a0, lambda0, if sigma0==0 then Nothing else Just sigma0, bias)
    (task, distpriors) <- evalRandIO $ vandistTask priorparams numtrials proplab
    
    let prior  = (dirichletProcess alphaparam, distpriors)
    
    -- Now run the model
    (guesses, partition) <- evalRandIO $ andersonSample encoding prior task
    
    -- Dump all those mabies
    let demabify = V.map (fromMaybe (-9))
        demabified = V.map demabify task
        withguesses = V.zipWith V.snoc demabified guesses
        results = map (\(x,Just y) -> V.toList $ V.snoc x (fromIntegral y)) $ (filter (isJust . snd)) $ zip (V.toList withguesses) (V.toList partition)
        resultstrings = map (map show) results
    
    -- Print it out
    putStrLn $ printCSV $ map ("STIM":) resultstrings
    
    putStrLn $ printCSV $ map (("CLUST":) . (map show)) $ summarizeClusters distpriors task partition

data ModelArgs = TVTask {
                   alphaparam :: Double,
                   order :: SortOrder,
                   encoding :: Encoding,
                   sigma0 :: Double,
                   a0 :: Double,
                   lambda0 :: Double,
                   bias :: Double,
                   nlabarg :: Int
                 } | Vandist {
                   alphaparam :: Double,
                   encoding :: Encoding,
                   sigma0 :: Double,
                   a0 :: Double,
                   lambda0 :: Double,
                   bias :: Double,
                   numtrials :: Int,
                   proplab :: Double
                 }
                 deriving (Data, Show, Typeable)

tvMode = TVTask {
           alphaparam = 1         &= help "Dirichlet parameter (alpha)",
           order = enum [Interspersed, 
                         LabeledFirst, 
                         LabeledLast],
           encoding = enum [EncodeActual, 
                            EncodeGuess, 
                            EncodeGuessSoft],
           sigma0 = 0             &= help "Prior over sd",  -- Could use .35 as well
           a0 = 1                 &= help "Strength of prior over mean",
           lambda0 = 1            &= help "Strength of prior over variance",
           bias = 0               &= help "Bias. Valence determines direction",
           nlabarg = 16           &= help "# of labeled items, if negative uses 16-all-labeled"
           } 
           &= help "Run the McDonnell (2013) TV task"

vandistMode = Vandist {
           alphaparam = 1         &= help "Dirichlet parameter (alpha)",
           encoding = enum [EncodeActual, 
                            EncodeGuess, 
                            EncodeGuessSoft],
           sigma0 = 0             &= help "Prior over sd",
           a0 = 1                 &= help "Strength of prior over mean",
           lambda0 = 1            &= help "Strength of prior over variance",
           bias = 0               &= help "Bias. Valence determines direction",
           numtrials = 800        &= help "Number of trials.",
           proplab = 0.5          &= help "Proporation of items labeled"
           } 
           &= help "Run with Vandist (2009) task"

modelArgs = cmdArgsMode 
           $ modes [tvMode, vandistMode]
           &= program "./runanderson" 
           &= summary "Anderson Rational Model" 
           -- &= details ["Runs a simulation of the Anderson Rational Model in the given task."]

main :: IO ()
main = do
    opts <- cmdArgsRun modelArgs
    case opts of TVTask{..} -> runTVTask opts
                 otherwise ->  runVandistTask opts
    -- runMedinSchaffer
    -- runContinuous
    -- runZeithamova
    -- runTVTask opts
    -- runVandistTask opts

