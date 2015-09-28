{-# LANGUAGE CPP #-}
module Main (main) where
#if !MIN_VERSION_base(4, 8, 0)
import Data.Functor ((<$>))
#endif
import Data.Foldable (for_)
import Data.List (intercalate)
import Data.Monoid ((<>))
import Data.Ratio (denominator, numerator)
import System.Exit (exitFailure)
import System.IO (IOMode(WriteMode), hFlush, hPutStrLn,
                  stderr, stdout, withFile)
import Crypto.Hash (Digest, MD5, hashlazy)
import qualified Data.ByteString.Lazy as ByteStringL
import WignerSymbols
import WignerSymbols.Internal

hexMD5 :: ByteStringL.ByteString -> String
hexMD5 s = show (hashlazy s :: Digest MD5)

checkResults :: [(Int, String)]
             -> Int
             -> String
             -> (Int -> (String -> IO ()) -> IO ())
             -> IO ()
checkResults knownHashes tjMax name compute = do
  withFile filename WriteMode $ \ h ->
    compute tjMax (hPutStrLn h)
  newHash <- hexMD5 <$> ByteStringL.readFile filename
  case lookup tjMax knownHashes of
    Nothing -> do
      hPutStrLn stderr (errorPrefix <> "no previous hash available")
      hPutStrLn stderr ("current hash is: " <> newHash)
      hFlush stderr
      exitFailure
    Just oldHash
      | oldHash == newHash -> do
          putStrLn (okPrefix <> name <> ": hash matched (" <> oldHash <> ")")
          hFlush stdout
      | otherwise -> do
          hPutStrLn stderr (errorPrefix <> name <> ": hash does not match!")
          hPutStrLn stderr ("previous hash is: " <> oldHash)
          hPutStrLn stderr ("current hash is:  " <> newHash)
          hFlush stderr
          exitFailure
  where filename = "dist/" <> name <> "-tj" <> show tjMax <> ".txt"
        errorPrefix = "[\ESC[31;1merror\ESC[0m] "
        okPrefix = "[\ESC[32mok\ESC[0m] "

main :: IO ()
main = do

  checkResults knownHashes_cg 15 "cg" $ \ tjMax write ->
    for_ (get3tjms tjMax) $ \ (tj1, tm1, tj2, tm2, tj3, tm3) ->
      let r = clebschGordanSq (tj1, tm1, tj2, tm2, tj3, -tm3) in
      write $
        intercalate "\t" (show <$> [tj1, tm1, tj2, tm2, tj3, -tm3]) <>
        "\t" <> show (ssr_signum r * ssr_numerator r) <>
        "/" <> show (ssr_denominator r)

  checkResults knownHashes_cg 15 "w3j" $ \ tjMax write ->
    for_ (get3tjms tjMax) $ \ (tj1, tm1, tj2, tm2, tj3, tm3) ->
      let SignedSqrtRatio wr = wigner3jSq (tj1, tm1, tj2, tm2, tj3, tm3)
          r = wr
            * (fromIntegral tj3 + 1)
            / (-1) ^^ ((tj1 - tj2 - tm3) `div` 2) in
      write $
        intercalate "\t" (show <$> [tj1, tm1, tj2, tm2, tj3, -tm3]) <>
        "\t" <> show (numerator r) <>
        "/" <> show (denominator r)

knownHashes_cg :: [(Int, String)]
knownHashes_cg =
  [ (15, "9192023f26dae0eebcce11afa7372eb6")
  , (20, "75ef56391b61e1bb2336e36ac7834216")
  , (25, "5901128892a264b73b5479b70b331fd0")
  , (30, "75ef56391b61e1bb2336e36ac7834216")
  , (40, "2f9b936ea977249c1fea8a22d190a4cf")
  ]
