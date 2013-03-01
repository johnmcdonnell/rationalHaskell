
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
import System.Environment (getArgs)
import Text.CSV

import Utils

import Stats
import Rational
import Anderson
import Tasks

-- Returns the likelihood of "0" as the first argument, and the most likely cluster as the second.
runTest :: (ClusterPrior, [PDFFromSample]) -> Stims -> Partition -> Stims -> Vector (Double, Int)
runTest prior stimuli assignments = V.map getpred
  where
    getpred = (second (argmax)) . (first head) . (infer prior stimuli assignments)

runTVTask :: [String] -> IO ()
runTVTask args = do
    let cparam = if length args > 0 then read (args!!0) else 1
        maxlab = 16
        (nlab, nounlab) = if length args > 1 then (\x -> if x<0 then (maxlab, True) else (x, False) ) $ read (args!!1) 
                                             else (maxlab, False)
        orderarg = if length args > 2 then (args!!2) else "interspersed"
        encodearg = if length args > 3 then (args!!3) else "actual"
    
        order = case orderarg of "interspersed" -> Interspersed
                                 "labfirst"     -> LabeledFirst
                                 "lablast"      -> LabeledLast
                                 otherwise      -> error $ "Inappropriate order: " ++ orderarg ++ "; order should be one of interspersed, labfirst, lablast."
        encoding = case encodearg of "guess"    -> EncodeGuess
                                     "softguess" -> EncodeGuessSoft
                                     otherwise  -> EncodeActual
    
    -- Set up priors
    let filterfun = if nounlab then V.filter (isJust . V.last) else id
    (task, distpriors) <- evalRandIO $ first filterfun <$> mcdonnellTaskOrdered order (1, 1, (-1)) (28*4) nlab
    -- print $ sortBy (compare `on` fst) $ map (\((Just bimod):(Just unimod):xs) -> (bimod, unimod)) task
    
    let prior  = (dirichletProcess cparam, distpriors)
    
    -- Now run the model
    (partition, guesses) <- evalRandIO $ andersonSample encoding prior task
    
    -- Dump all those mabies
    let demabify = V.map (fromMaybe (-9))
        demabified = V.map demabify task
        results = map (\(x,Just y) -> V.toList $ V.snoc x (fromIntegral y)) $ (filter (isJust . snd)) $ zip (V.toList demabified) (V.toList partition)
        resultstrings = map (map show) results
    
    -- Print it out
    putStrLn $ printCSV $ map ("STIM":) resultstrings
    
    putStrLn $ printCSV $ map (("CLUST":) . (map show)) $ summarizeClusters distpriors task partition
    
    -- If there were no labels, we have to assume they are associated with the largest clusters.
    let initialinferences =  runTest prior task partition gridTest
        inferpartition = map snd $ V.toList initialinferences
        partitions = Map.assocs $ countUnique inferpartition
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

main = do
    args <- getArgs
    -- runMedinSchaffer
    -- runContinuous
    -- runZeithamova
    runTVTask args
    -- runVandistTask args
