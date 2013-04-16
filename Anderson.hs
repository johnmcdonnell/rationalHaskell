{-# LANGUAGE DeriveDataTypeable #-}

module Anderson (andersonSample,
                ) where

import Data.Function
import Data.List
import Data.Maybe
import Data.Vector (freeze, thaw, MVector, Vector, (!))
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM
import Control.Monad
import Control.Applicative
import Control.Monad.Random
import Control.Monad.Reader
import Control.Monad.State
import Control.Monad.ST
import System.Random
import Statistics.Sample
import System.Console.CmdArgs

import Utils
import Types
import Rational

-- Anderson sampling
sampleNext :: (ClusterPrior, [PDFFromSample]) -> Stims -> Partition -> Stim -> Int
sampleNext (clusterPrior, distributions) stimuli assignments newstim = assignment
  where 
    assignment = (V.maxIndex . V.fromList) posterior
    posterior = clusterPosterior clusterPrior distributions stimuli assignments newstim

-- * Types and functions for changing the encoding, not currently implemented (always encodeActual)
encodeActual :: Stim -> [Double] -> Rand StdGen Stim
encodeActual newstim guess = return $ newstim

encodeGuess :: Stim -> [Double] -> Rand StdGen  Stim
encodeGuess newstim [x] = return $ V.snoc (V.init newstim) (if x>0.5 then Just 1 else (if x<0.5 then Just 0 else Nothing))
encodeGuess newstim _   = return newstim

-- TODO: Make sure the mapping is correct! Inference prob of .8 should mean .8 chance of being 1!
encodeGuessSoft :: Stim -> [Double] -> Rand StdGen  Stim
encodeGuessSoft newstim []        = return newstim
encodeGuessSoft newstim inference = encode <$> (getRandom :: Rand StdGen Double)
  where
    encode = (V.snoc (V.init newstim) . makechoice)
    makechoice needle = if prob > needle then Just 1 else Just 0
    prob = case inference of [x] -> x
                             otherwise -> 0.5


-- | Accept a particular stimulus and its index, update model and return a guess of what its label will be
andersonIterate :: Encoding -> (Int, Stim) -> Model Double
andersonIterate encoding  (i, newstim) = do
    (prior, stims) <- ask
    assignments <- get
    -- The ccode is here for this but it's not implemented!
    let encodefun = case encoding of EncodeGuess -> encodeGuess
                                     EncodeGuessSoft -> encodeGuessSoft 
                                     EncodeActual -> encodeActual
    
    let incomplete = V.any isNothing newstim
    
    -- Assume the last is the label
    let guessstim = V.snoc (V.init newstim) Nothing
        guess = fst $ infer prior stims assignments guessstim
    
    let chosenclust = sampleNext prior stims assignments newstim
    put $ V.modify (\vec -> VM.unsafeWrite vec i (Just chosenclust)) assignments
    return $ head guess

-- Not yet implemented
summarizeCluster :: Stims -> Partition -> Int -> String
summarizeCluster allstims part i = "hi"
  where
    stims = V.map (allstims!) indices
    indices = V.elemIndices (Just i) part

summarizeClusters :: Stims -> Partition -> String
summarizeClusters stims partition = concat $ intersperse "\n" $ map summfun [0..nclusts-1]
  where
    summfun = summarizeCluster stims partition
    nclusts = length $ (nub . catMaybes) $ partList
    partList = V.toList partition

andersonSample ::   Encoding                  -- ^ How to encode items (Use guess or not)
                 -> (ClusterPrior, [PDFFromSample])  -- ^ (Coupling param, estimators for each dimension)
                 -> Stims                     -- ^ Task stimuli
                 -> Rand StdGen (V.Vector Double, Partition) -- ^ Partition and list of guesses
andersonSample encoding prior stimuli = runStateT (runReaderT (runM model) (prior, stimuli)) assignmentStore
  where
    model = liftM V.fromList $ forM enumStimuli (andersonIterate encoding)
    enumStimuli = zip [0..] (V.toList stimuli)
    assignmentStore = V.replicate n Nothing
    n = V.length stimuli

