{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

{-# LANGUAGE TypeFamilies                  #-}

module StochVolMain (
    bigM
  , runMC
  , barDiag
  , displayHeader
  , hist
  , main
  ) where

import Control.Applicative

import qualified Data.Vector.Unboxed as U

import Data.Histogram ( asList )
import Data.Histogram.Fill
import Data.Histogram.Generic ( Histogram )

import qualified Control.Foldl as F

import Data.Colour
import Diagrams.Prelude hiding ( sample, render )

import Diagrams.Backend.Cairo.CmdLine
import Diagrams.Backend.CmdLine

import Graphics.Rendering.Chart.Backend.Diagrams
import Graphics.Rendering.Chart

import Data.Default.Class

import Text.Printf

import System.IO.Unsafe

import StochVol


denv :: DEnv
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 500 500

numBins :: Int
numBins = 40

stats :: (F.Foldable f, Fractional a) =>
         f a -> (a, a, a)
stats v = F.fold stats v
  where
    stats = f <$> (F.premap (\x -> x * x) F.sum) <*> F.sum <*> F.genericLength
    f x2Sum xSum n = (var, mean, n)
      where
        mean = xSum / n
        mean2 = x2Sum / n
        var = n * (mean2 - mean * mean) / (n - 1)

hb :: F.Foldable f =>
      f Double -> HBuilder Double (Histogram U.Vector BinD Double)
hb xs = forceDouble -<< mkSimple (binD lower numBins upper)
  where
    (varX, xBar, _) = stats xs
    lower = xBar - 2.0 * sqrt varX
    upper = xBar + 2.0 * sqrt varX

hist :: F.Foldable f =>
        f Double -> Histogram U.Vector BinD Double
hist xs = fillBuilder (hb xs) xs

displayHeader :: FilePath -> Diagram B R2 -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

barChart :: String ->
            [(Double, Double)] ->
            Graphics.Rendering.Chart.Renderable ()
barChart title bvs = toRenderable layout
  where
    layout =
      layout_title .~ title
      $ layout_x_axis . laxis_generate .~ autoIndexAxis (map (printf "%4.3f" . fst) bvs)

      $ layout_y_axis . laxis_title .~ "Frequency"
      $ layout_plots .~ (map plotBars $ [bars1])
      $ def

    bars1 =
      plot_bars_titles .~ [title]
      $ plot_bars_values .~ addIndexes (map return $ map snd bvs)
      $ plot_bars_style .~ BarsClustered
      $ plot_bars_item_styles .~ [(solidFillStyle (blue `withOpacity` 0.25), Nothing)]
      $ def

barDiag :: String ->
           [(Double, Double)] ->
           Diagram B R2
barDiag title bvs = fst $ runBackend denv (render (barChart title bvs) (500, 500))

main :: IO ()
main = do
  let result = runMC
      mus = map (fst. fst) result
      phis = map (snd . fst) result
      taus = map snd result
  putStrLn $ show $ (/(fromIntegral bigM)) $ sum mus
  putStrLn $ show $ (/(fromIntegral bigM)) $ sum phis
  putStrLn $ show $ (/(fromIntegral bigM)) $ sum taus
  displayHeader "mus.png"
    (barDiag "Mu" (zip (map fst $ asList (hist mus)) (map snd $ asList (hist mus))))
  displayHeader "phis.png"
    (barDiag "Phi" (zip (map fst $ asList (hist phis)) (map snd $ asList (hist phis))))
  displayHeader "taus.png"
    (barDiag "Tau" (zip (map fst $ asList (hist taus)) (map snd $ asList (hist taus))))
