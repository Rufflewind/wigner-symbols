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

main :: IO ()
main = do

  checkResults knownHashes_cg 25 "cg" $ \ tjMax write ->
    for_ (get3tjms tjMax) $ \ (tj1, tm1, tj2, tm2, tj3, tm3) ->
      let r = clebschGordanSq (tj1, tm1, tj2, tm2, tj3, -tm3) in
      write . intercalate "\t" $
        (show <$> [tj1, tm1, tj2, tm2, tj3, -tm3]) <>
        [show (ssr_signum r * ssr_numerator r) <>
         "/" <> show (ssr_denominator r)]

  checkResults knownHashes_cg 5 "cgw" $ \ tjMax write ->
    for_ (get3tjms tjMax) $ \ (tj1, tm1, tj2, tm2, tj3, tm3) ->
      let SignedSqrtRatio wr = wigner3jSq (tj1, tm1, tj2, tm2, tj3, tm3)
          r = wr
            * (fromIntegral tj3 + 1)
            / (-1) ^^ ((tj1 - tj2 - tm3) `div` 2) in
      write . intercalate "\t" $
        (show <$> [tj1, tm1, tj2, tm2, tj3, -tm3]) <>
        [show (numerator r) <> "/" <> show (denominator r)]

  checkResults knownHashes_w6j 15 "w6j" $ \ tjMax write ->
    for_ (get6tjs tjMax) $ \ tjs ->
      let r = wigner6jSq tjs in
      write . intercalate "\t" $
        (show <$> tuple6ToList tjs) <>
        [show (ssr_signum r * ssr_numerator r) <>
         "/" <> show (ssr_denominator r)]

  checkResults knownHashes_w9j 7 "w9j" $ \ tjMax write ->
    for_ (get9tjs tjMax) $ \ tjs ->
      let r = wigner9jSq tjs in
      write . intercalate "\t" $
        (show <$> tuple9ToList tjs) <>
        [show (ssr_signum r * ssr_numerator r) <>
         "/" <> show (ssr_denominator r)]

knownHashes_cg :: [(Int, String)]
knownHashes_cg =
  [ (5,  "e74c501299b456a6cb29e4f5714e9061") --     681
  , (10, "b6d0770101f4ebdaa9a55d94f07b001f") --   11487
  , (15, "9192023f26dae0eebcce11afa7372eb6") --   69272
  , (20, "75ef56391b61e1bb2336e36ac7834216") --  259523
  , (25, "5901128892a264b73b5479b70b331fd0") --  737113
  , (30, "75ef56391b61e1bb2336e36ac7834216") -- 1747984
  , (40, "2f9b936ea977249c1fea8a22d190a4cf") -- 6931995
  ]

knownHashes_w6j :: [(Int, String)]
knownHashes_w6j =
  [ (5,  "26c24e568fc96f1732ebb3130a46f22a") --    1479
  , (10, "f892f4b466e0558179ca870941d0a456") --   42393
  , (15, "f50b0163194cef1699727b7064760ec0") --  363196
  , (20, "e1b5dad0f1469cc54b6139533f982815") -- 1766270
  , (25, "f326bf6e12a94120d2f46582e95e92f8") -- 6171698
  ]

knownHashes_w9j :: [(Int, String)]
knownHashes_w9j =
  [ (3,  "4005ef20e2ed8c789917dce99d027bc4") --    1616
  , (4,  "92cfc13320e7fd6a34b3970ebef58e06") --    9060
  , (5,  "d596fa3960aafae148754b6f3274507d") --   38031
  , (7,  "7b338708ef3aa4ba0a4f5bd8c8b4e6aa") --  401899
  , (10, "479c0a020eaceff5539e2dda2200c1ab") -- 5898846
  ]

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
      hPutStrLn stderr (errorPrefix <> "no known hash available")
      hPutStrLn stderr ("actual: " <> newHash)
      hFlush stderr
      exitFailure
    Just oldHash
      | oldHash == newHash -> do
          putStrLn (okPrefix <> name <> ": hash matched (" <> oldHash <> ")")
          hFlush stdout
      | otherwise -> do
          hPutStrLn stderr (errorPrefix <> name <> ": hash does not match!")
          hPutStrLn stderr ("expected: " <> oldHash)
          hPutStrLn stderr ("actual:   " <> newHash)
          hFlush stderr
          exitFailure
  where filename = "dist/" <> name <> "-tj" <> show tjMax <> ".txt"
        errorPrefix = "[\ESC[31;1merror\ESC[0m] "
        okPrefix = "[\ESC[32mok\ESC[0m] "

hexMD5 :: ByteStringL.ByteString -> String
hexMD5 s = show (hashlazy s :: Digest MD5)
