module Main (main) where
import Criterion.Main

main :: IO ()
main =
  defaultMain
  [ bench "example" $ whnf (+ 1) 1
  ]

module Main (main) where
import Control.Monad (when)
import Data.Foldable (for_)
import Data.Functor
import Data.List (intercalate)
import Data.Monoid ((<>))
import Data.Ratio (numerator, denominator)
import System.Exit (exitFailure)
import System.IO (IOMode(WriteMode), hPutStrLn, stderr, withFile)
import Prelude
import Crypto.Hash (Digest, MD5, hashlazy)
import qualified Data.ByteString.Lazy as ByteStringL
import WignerSymbols

hexMD5 :: ByteStringL.ByteString -> String
hexMD5 s = show (hashlazy s :: Digest MD5)

main :: IO ()
main = do
  withFile filename WriteMode $ \ h ->
    for_ [0 .. maxTj] $ \ tj1 ->
    for_ [0 .. maxTj] $ \ tj2 ->
    for_ (filter (<= maxTj) (getTj12s tj1 tj2)) $ \ tj12 ->
    for_ (getTms tj1) $ \ tm1 ->
    for_ (getTms tj2) $ \ tm2 ->
      let tm12 = tm1 + tm2
          SignedSqrtRatio r = clebschGordanSq tj1 tm1 tj2 tm2 tj12 tm12
      in when (abs tm12 <= tj12) $
         hPutStrLn h $
         intercalate "\t" (show <$> [tj1, tm1, tj2, tm2, tj12, tm12]) <>
         "\t" <> show (numerator r) <> "/" <> show (denominator r)
  newHash <- hexMD5 <$> ByteStringL.readFile filename
  case lookup maxTj knownHashes of
    Nothing -> do
      hPutStrLn stderr "***[error] no previous hash available"
      hPutStrLn stderr ("current hash is: " <> newHash)
      exitFailure
    Just oldHash
      | oldHash == newHash -> putStrLn "[ok] hash matched"
      | otherwise -> do
        hPutStrLn stderr "***[error] hash does not match!"
        hPutStrLn stderr ("previous hash is: " <> oldHash)
        hPutStrLn stderr ("current hash is:  " <> newHash)
        exitFailure

  where

    maxTj :: Int
    maxTj = 25

    filename = "dist/cg-tj" <> show maxTj <> ".txt"

    knownHashes =
      [ (15, "9192023f26dae0eebcce11afa7372eb6")
      , (20, "75ef56391b61e1bb2336e36ac7834216")
      , (25, "5901128892a264b73b5479b70b331fd0")
      , (30, "75ef56391b61e1bb2336e36ac7834216")
      , (40, "2f9b936ea977249c1fea8a22d190a4cf")
      ]

-- tjmax = 15: 0.21 sec -> 3.0μs
-- tjmax = 20: 0.91 sec -> 3.5μs
-- tjmax = 25: 2.92 sec -> 4.0μs
-- tjmax = 30: 7.63 sec -> 4.4μs
-- tjmax = 40: 37.5 sec -> 5.4μs
