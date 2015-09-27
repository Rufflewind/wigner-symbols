module Main (main) where
import Data.Foldable (for_)
import Criterion.Main
import WignerSymbols
import WignerSymbols.Internal

main :: IO ()
main = do
  print (count <$> [5, 10, 15, 20, 25])
  defaultMain
    [ bench "cg.5" (whnfIO (doIt 5))
    , bench "cg.10" (whnfIO (doIt 10))
    , bench "cg.15" (whnfIO (doIt 15))
    , bench "cg.20" (whnfIO (doIt 20))
    , bench "cg.25" (whnfIO (doIt 25))
    ]

doIt :: Int -> IO ()
doIt tjMax = for_ (get3tjms tjMax) $ \ tjms ->
  (clebschGordan tjms :: Double) `seq` return ()

count :: Int -> Int
count tjMax = length (get3tjms tjMax)

-- tjmax = 15: 0.21 sec -> 3.0μs
-- tjmax = 20: 0.91 sec -> 3.5μs
-- tjmax = 25: 2.92 sec -> 4.0μs
-- tjmax = 30: 7.63 sec -> 4.4μs
-- tjmax = 40: 37.5 sec -> 5.4μs
