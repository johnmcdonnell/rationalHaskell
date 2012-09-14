
module Anderson (andersonSample) where

import Data.Function
import Data.List
import Data.Maybe
import Control.Monad
import qualified Data.IORef as IORef

import Stats
import Rational

-- Anderson sampling
sampleNext :: (ClusterPrior, [PDFFromSample]) -> [Stim] -> Partition -> Stim -> Int
sampleNext (clusterPrior, distributions) stimuli assignments newstim = assignment
  where 
    assignment = fst $ maximumBy (compare `on` snd) $ zip [0..] posterior
    posterior = clusterPosterior clusterPrior distributions stimuli assignments newstim

andersonSample :: (ClusterPrior, [PDFFromSample])  -- ^ (Coupling param, estimators for each dimension)
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


