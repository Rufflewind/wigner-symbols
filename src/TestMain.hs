{-# LANGUAGE CPP #-}
module Main (main) where
#if !MIN_VERSION_base(4, 8, 0)
import Data.Functor ((<$>))
#endif
import Data.Foldable (for_)
import Data.List (intercalate)
import Data.Monoid ((<>))
import System.Exit (exitFailure)
import System.IO (IOMode(WriteMode), hPutStrLn, stderr, withFile)
import Crypto.Hash (Digest, MD5, hashlazy)
import qualified Data.ByteString.Lazy as ByteStringL
import WignerSymbols
import WignerSymbols.Internal

hexMD5 :: ByteStringL.ByteString -> String
hexMD5 s = show (hashlazy s :: Digest MD5)

main :: IO ()
main = do
  withFile filename WriteMode $ \ h ->
    for_ (get3tjms tjMax) $ \ (tj1, tm1, tj2, tm2, tj3, tm3) ->
      let r = clebschGordanSq (tj1, tm1, tj2, tm2, tj3, -tm3) in
      hPutStrLn h $
        intercalate "\t" (show <$> [tj1, tm1, tj2, tm2, tj3, -tm3]) <>
        "\t" <> show (ssr_signum r * ssr_numerator r) <>
        "/" <> show (ssr_denominator r)
  newHash <- hexMD5 <$> ByteStringL.readFile filename
  case lookup tjMax knownHashes of
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

    tjMax :: Int
    tjMax = 25

    filename = "dist/cg-tj" <> show tjMax <> ".txt"

    knownHashes =
      [ (15, "9192023f26dae0eebcce11afa7372eb6")
      , (20, "75ef56391b61e1bb2336e36ac7834216")
      , (25, "5901128892a264b73b5479b70b331fd0")
      , (30, "75ef56391b61e1bb2336e36ac7834216")
      , (40, "2f9b936ea977249c1fea8a22d190a4cf")
      ]
