{-# LANGUAGE NoMonomorphismRestriction #-}

    module Main (
  main
  ) where

import qualified Data.Vector.Unboxed as V
import Data.Random.Source.PureMT
import Data.Random
import Control.Monad.State


nSamples :: Int
nSamples = 100

beta :: Double
beta = 2.0

testXs :: Int -> V.Vector Double
testXs m =
  V.fromList $
  evalState (replicateM m (sample StdUniform))
  (pureMT 2)

testEpsilons :: Int -> V.Vector Double
testEpsilons m =
  V.fromList $
  evalState (replicateM m (sample StdNormal))
  (pureMT 2)

testYs = V.zipWith (\x e -> beta * x + e) (testXs nSamples) (testEpsilons nSamples)

main = putStrLn $ show testYs