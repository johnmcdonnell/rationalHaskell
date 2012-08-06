
import Control.Monad.State (get)
import Data.Function (on)
import Data.List (maximumBy)
import Data.Maybe

import Stats
import Rational


-- * Hardcoded presets
dirichlet_alpha = 0.5

stimuli = map (\x -> [x]) [0.2,0.3,0.8] :: [Stim]
andersonstims = [[1,1,1,1,1], [1,0,1,0,1], [1,0,1,1,0], [0,0,0,0,0], [0,1,0,1,1], [0,1,0,0,0]]
multinom_prior = multinomialDistribution 1 [1,1]
andersondists = replicate 5 multinom_prior


-- Anderson sampling
sampleNext :: [PDFFromSample] -> [Stim] -> Partition -> Stim -> Int
sampleNext distributions stimuli assignments newstim = assignment
  where 
    assignment = fst $ maximumBy (compare `on` snd) $ zip [0..] posterior
    posterior = clusterPosterior distributions stimuli assignments newstim

post = clusterPosterior andersondists andersonstims [Just 0] [1,0,1,0,1]
selection = sampleNext andersondists andersonstims [Just 0] [1,0,1,0,1]

