{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

{-# LANGUAGE ViewPatterns                 #-}

module Main (
    main
    ) where

import Diagrams.Prelude
import Diagrams.Backend.CmdLine
import Diagrams.Backend.Cairo.CmdLine

import qualified Data.Vector as V
import qualified Data.ByteString.Lazy as BL
import Data.Csv
import Data.Time
import System.Locale
import Data.Char
import qualified Data.ByteString as B
import Control.Monad

import Autoregression
import AutoregressionChart

displayHeader :: FilePath -> Diagram B R2 -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

data RiverObs = RiverObs
               {  riverObsDate :: UTCTime
               ,  riverObsFlow :: LaxDouble
               ,  riverObsTemp :: Maybe LaxDouble
               , _empty1       :: Maybe String
               , _empty2       :: Maybe String
               , _empty3       :: Maybe String
               , _empty4       :: Maybe String
               }
  deriving Show

newtype LaxDouble = LaxDouble { laxDouble :: Double }
  deriving Show

instance FromField LaxDouble where
  parseField = fmap LaxDouble . parseField . addLeading

    where

      addLeading :: B.ByteString -> B.ByteString
      addLeading bytes =
            case B.uncons bytes of
              Just (c -> '.', _)    -> B.cons (o '0') bytes
              Just (c -> '-', rest) -> B.cons (o '-') (addLeading rest)
              _ -> bytes

      c = chr . fromIntegral
      o = fromIntegral . ord

instance FromField UTCTime where
  parseField s = do
    f <- parseField s
    case parseTime defaultTimeLocale "%d-%b-%y" f of
      Nothing -> fail "Unable to parse UTC time"
      Just g  -> return g

instance FromRecord RiverObs where
  parseRecord v
         | V.length v == 7
         = RiverObs <$>
           v .!  0 <*>
           v .!  1 <*>
           v .!  2 <*>
           v .!  3 <*>
           v .!  4 <*>
           v .!  5 <*>
           v .!  6
         | otherwise = mzero

main :: IO ()
main = do
  riverObsCsv <- BL.readFile "EnvAgencyFlowData.csv"
  let records :: Either String (V.Vector RiverObs)
      records = decode HasHeader riverObsCsv
  case records of
    Left err -> putStrLn err
    Right rivObs-> do
      let flows = V.map (log . laxDouble . riverObsFlow) rivObs
          preds = predictions flows
      displayHeader "diagrams/AutoregressionVary1.png"
                    (diag "Predicted Flow at Kingston Bridge (Varying Parameters)"
                           (zip [0..] (V.toList $ V.map exp $ V.take 800 flows))
                           (zip [0..] (V.toList $ V.map exp $ V.take 800 preds)))
      putStrLn "hello"
