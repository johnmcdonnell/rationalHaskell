
module Random (R,
               runRandom,
               runRandomIO,
               rand,
               normals,
               shuffleR,
               shuffleNR
               ) where

import Control.Monad
import Control.Monad.State (State, evalState, get, put)
import System.Random (Random, StdGen, getStdGen, mkStdGen, random)
import System.Random.Shuffle
import Control.Applicative ((<$>))

-- Thanks to: http://stackoverflow.com/questions/2110535/sampling-sequences-of-random-numbers-in-haskell/2125329#2125329
-- for a start on this.

type R a = State StdGen a

runRandom :: R a -> Int -> a
runRandom action seed = evalState action $ mkStdGen seed

runRandomIO :: R a -> IO a
runRandomIO action = (return . (evalState action)) =<< getStdGen

rand :: (Random a) => R a
rand = do
  gen <- get
  let (r, gen') = random gen
  put gen'
  return r

randPair :: (Random a) => R (a, a)
randPair = do
  x <- rand
  y <- rand
  return (x,y)

-- Random int between lower and (upper-1)
randInt :: Int -> Int ->  R Int
randInt lower upper = do
  let range = upper - lower
  x <- rand :: R Double
  return $ floor $ (x * (fromIntegral range)) + (fromIntegral lower)

boxMuller :: Double -> Double -> (Double, Double) -> Double
boxMuller mu sigma (r1,r2) =  mu + sigma * sqrt (-2 * log r1) * cos (2 * pi * r2)

normals :: Double -> Double -> R [Double]
normals mu sd = mapM (\_ -> boxMuller mu sd <$> randPair) $ repeat ()

shuffleNR :: [a] -> Int -> R [a]
shuffleNR xs n =
    (shuffle xs) <$> (sequence $ map (\i -> randInt 0 (n-i) ) [1..(n-1)])

shuffleR :: [a] -> R [a]
shuffleR xs = shuffleNR xs n
  where
    n = length xs



