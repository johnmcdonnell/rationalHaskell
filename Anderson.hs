
module Anderson (andersonSample,
                Encoding (EncodeActual, EncodeGuess)) where

import Data.Function
import Data.List
import Data.Maybe
import Control.Monad
import qualified Data.IORef as IORef

import Stats
import Rational

import Control.Applicative
import Control.Monad.ST
import Data.Vector (freeze, thaw, MVector, Vector, (!))
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM

import JohnFuns

-- Anderson sampling
sampleNext :: (ClusterPrior, [PDFFromSample]) -> Stims -> Partition -> Stim -> Int
sampleNext (clusterPrior, distributions) stimuli assignments newstim = assignment
  where 
    assignment = (V.maxIndex . V.fromList) posterior
    posterior = clusterPosterior clusterPrior distributions stimuli assignments newstim

-- andersonSample :: (ClusterPrior, [PDFFromSample])  -- ^ (Coupling param, estimators for each dimension)
--                -> Stims                     -- ^ Task stimuli
--                -> Partition
-- andersonSample prior stimuli = runST $ do
--     let n = V.length stimuli
--     assignmentStore <- VM.replicate n Nothing
--     V.forM_ (V.indexed stimuli) (\(i, newstim) -> do
--         assignments <- V.freeze assignmentStore
--         let chosenclust = sampleNext prior stimuli assignments newstim
--         VM.write assignmentStore i (Just chosenclust)
--         )
--     V.freeze assignmentStore

encodeActual :: Stim -> [Double] -> Stim
encodeActual newstim guess = newstim

encodeGuess :: Stim -> [Double] -> Stim
encodeGuess newstim guess = case guess of [x] -> V.snoc (V.init newstim) (if x>0.5 then Just 1 else (if x<0.5 then Just 0 else Nothing))
                                          otherwise -> newstim
data Encoding = EncodeActual | EncodeGuess

andersonIterate :: (ClusterPrior, [PDFFromSample]) -> Encoding -> (Partition, Stims) -> (Int, Stim) -> (Partition, Stims)
andersonIterate prior encoding (assignments, stims) (i, newstim) = (retAssign, retStims)
  where
    retStims = V.snoc stims encodestim -- TODO would be better in place.
    retAssign = V.modify (\vec -> VM.unsafeWrite vec i (Just chosenclust)) assignments
    chosenclust = sampleNext prior stims assignments encodestim
    encodestim = (encodefun newstim guess)
    guess = infer (fst prior) (snd prior) stims assignments newstim
    encodefun = case encoding of EncodeGuess -> encodeGuess
                                 EncodeActual -> encodeActual

andersonSample ::    Encoding                  -- ^ How to encode items (Use guess or not)
                 -> (ClusterPrior, [PDFFromSample])  -- ^ (Coupling param, estimators for each dimension)
                 -> Stims                     -- ^ Task stimuli
                 -> Partition
andersonSample encoding prior stimuli = finalAssign
  where
    (finalAssign, finalStims) = V.foldl' (andersonIterate prior encoding) (assignmentStore, stimStore) (V.indexed stimuli)
    assignmentStore = V.replicate n Nothing
    stimStore = V.empty
    n = V.length stimuli

