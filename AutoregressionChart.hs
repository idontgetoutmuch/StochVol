{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

module AutoregressionChart (
    diag
  ) where

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Diagrams
import Diagrams.Backend.Cairo.CmdLine
import Diagrams.Prelude hiding ( render, Renderable )
import Data.Default.Class

import System.IO.Unsafe


denv :: DEnv
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 600 500

diag :: String ->
        [(Double, Double)] ->
        [(Double, Double)] ->
        Diagram Cairo R2
diag t xs ys =
  fst $ runBackend denv (render (chart t xs ys) (600, 500))


chart :: String ->
         [(Double, Double)] ->
         [(Double, Double)] ->
         Renderable ()
chart t obs preds = toRenderable layout
  where

    actuals = plot_lines_values .~ [obs]
              $ plot_lines_style  . line_color .~ opaque blue
              $ plot_lines_title .~ "Actual"
              $ plot_lines_style  . line_width .~ 1.0
              $ def

    predicts = plot_lines_values .~ [preds]
              $ plot_lines_style  . line_color .~ opaque red
              $ plot_lines_title .~ "Prediction"
              $ plot_lines_style  . line_width .~ 1.0
              $ def

    layout = layout_title .~ t
           $ layout_plots .~ [toPlot actuals, toPlot predicts]
           $ layout_y_axis . laxis_title .~ "Cubic Metres per Second"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "Days since 19-Feb-13"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def
