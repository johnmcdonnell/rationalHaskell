
-- import qualified Control.Monad.State as State
import Data.Function (on)
import Data.List (maximumBy)
import Data.Maybe
import System.Random
import Data.Random

import qualified Data.IORef as IORef
import Control.Monad (forM, forM_)
import Control.Applicative

import JohnFuns (debug)
import Stats
import Random
import Rational


-- Anderson sampling
sampleNext :: (Double, [PDFFromSample]) -> [Stim] -> Partition -> Stim -> Int
sampleNext (couplingParam, distributions) stimuli assignments newstim = assignment
  where 
    assignment = fst $ maximumBy (compare `on` snd) $ zip [0..] posterior
    posterior = clusterPosterior clusterPrior distributions stimuli assignments newstim
    clusterPrior = dirichletProcess couplingParam

andersonSample :: (Double, [PDFFromSample])  -- ^ (Coupling param, posterior estimators for each dimension)
               -> [Stim]                     -- ^ Task stimuli
               ->  IO Partition
andersonSample prior stimuli = do
    let n = length stimuli
    assignmentStore <- IORef.newIORef (replicate n Nothing)
    forM_ (zip [0..] stimuli) (\(i, newstim) -> do
        assignments <- IORef.readIORef assignmentStore
        -- print assignments
        let chosenclust = sampleNext prior stimuli assignments newstim
        IORef.modifyIORef assignmentStore (replaceAtIndex i (Just chosenclust))
        )
    final <- IORef.readIORef assignmentStore
    -- print "FINAL:"
    -- print final
    return final

medinSchaffer = do
    let medinSchafferTask = [[1,1,1,1,1], [1,0,1,0,1], [1,0,1,1,0], [0,0,0,0,0], [0,1,0,1,1], [0,1,0,0,0]]
    let multinom_prior = multinomialPosterior [1,1]
    let andersondists = replicate 5 multinom_prior
    partition <- andersonSample (1.0, andersondists) medinSchafferTask
    print partition

onedtask :: [(Double, Double)] -> Int -> IO [[Double]]
onedtask params n = do
    samples <- forM params (\(mu, sigma) -> (take n) <$> runRandomIO (normals mu sigma))
    let stims = map (\x -> [x]) $ concat samples
    runRandomIO $ shuffleNR stims (n * (length params))

continuous = do
    let mu1 = 3
    let sigma1 = 1
    let mu2 = -15
    let sigma2 = 1
    let n = 10
    putStrLn "Testing on:" 
    putStrLn $ "mu1=" ++ (show mu1) ++ " sigma1=" ++ (show sigma1) ++ " n=" ++ (show n)
    putStrLn $ "mu2=" ++ (show mu2) ++ " sigma2=" ++ (show sigma2) ++ " n=" ++ (show n)
    task <- onedtask [(mu1, sigma1), (mu2, sigma2)] n
    
    let tprior = (0, 1, 1, 1)
    let prior  = (0.5, [tPosterior tprior])
    partition <- andersonSample prior task
    forM_ (zip task (catMaybes partition)) print

main = do
    -- medinSchaffer
    continuous

