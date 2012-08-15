
-- import qualified Control.Monad.State as State
import Data.Function (on)
import Data.List (maximumBy)
import Data.Maybe

import qualified Data.IORef as IORef
import Control.Monad (forM_)

import JohnFuns (debug)
import Stats
import Rational


-- * Hardcoded presets
clusterPrior = dirichletProcess 1.0

stimuli = map (\x -> [x]) [0.2,0.3,0.8] :: [Stim]
andersonstims = [[1,1,1,1,1], [1,0,1,0,1], [1,0,1,1,0], [0,0,0,0,0], [0,1,0,1,1], [0,1,0,0,0]]
multinom_prior = multinomialDistribution [1,1]
andersondists = replicate 5 multinom_prior


-- Anderson sampling
sampleNext :: [PDFFromSample] -> [Stim] -> Partition -> Stim -> Int
sampleNext distributions stimuli assignments newstim = assignment
  where 
    assignment = fst $ maximumBy (compare `on` snd) $ zip [0..] posterior
    posterior = clusterPosterior clusterPrior distributions stimuli assignments newstim

andersonSample :: [PDFFromSample] -> [Stim] ->  IO Partition
andersonSample distributions stimuli = do
    assignmentStore <- IORef.newIORef [Nothing]
    forM_ (zip [0..] stimuli) (\(i, newstim) -> do
        assignments <- IORef.readIORef assignmentStore
        print assignments
        let chosenclust = sampleNext distributions stimuli assignments newstim
        IORef.modifyIORef assignmentStore (replaceAtIndex i (Just chosenclust))
        )
    final <- IORef.readIORef assignmentStore
    print "FINAL:"
    print final
    return final


main = do
    partition <- andersonSample andersondists andersonstims
    print partition

