import StochVol

main :: IO ()
main = do
  let result = runMC
      mus = map (fst. fst) result
      phis = map (snd . fst) result
      taus = map snd result
  putStrLn $ show $ (/(fromIntegral bigM)) $ sum mus
  putStrLn $ show $ (/(fromIntegral bigM)) $ sum phis
  putStrLn $ show $ (/(fromIntegral bigM)) $ sum taus

