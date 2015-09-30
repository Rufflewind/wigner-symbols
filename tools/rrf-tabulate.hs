-- This code tabulates Wigner symbols using the Root Rational Fraction (RRF)
-- program by Anthony Stone and Charles Wood.  It is used here for
-- verification purposes.
--
-- http://www-stone.ch.cam.ac.uk/documentation/rrf
-- http://www-stone.ch.cam.ac.uk/pub/rrf-4.0.tgz
--
module Main (main) where
import Data.Char (isDigit)
import Data.Monoid ((<>))
import Data.List (intercalate)
import Data.Ratio ((%), denominator, numerator)
import Data.Foldable (for_)
import Foreign
import Foreign.C
import Prelude hiding (pi)
import System.IO (IOMode(WriteMode), hPutStrLn, withFile)
import Text.ParserCombinators.ReadP (ReadP)
import qualified Text.ParserCombinators.ReadP as P
import WignerSymbols.Internal

infixl 3 <++
(<++) :: ReadP a -> ReadP a -> ReadP a
(<++) = (P.<++)

foreign import ccall "__rrf_module_MOD_init_rrf" f_init_rrf
  :: Ptr CInt -> Ptr CInt -> IO ()

foreign import ccall "__wigner_MOD_ninej" f_ninej
  :: Ptr CInt -> Ptr CInt -> Ptr CInt
  -> Ptr CInt -> Ptr CInt -> Ptr CInt
  -> Ptr CInt -> Ptr CInt -> Ptr CInt
  -> Ptr CInt -> Ptr (Ptr ()) -> IO ()

-- sum, string, p, 100
foreign import ccall "__rrf_module_MOD_char4i" f_char4i
  :: Ptr (Ptr ()) -> Ptr CChar -> Ptr CInt -> Ptr CInt -> IO ()

main :: IO ()
main = do
  initRRFLib numOfPrimes
  rrf <- rrf_new maxStrLen

  tabulate 7 "w9j" $ \ tjMax write ->
    for_ (get9tjs tjMax) $ \ tjs -> do
      rrf_ninej rrf tjs
      r <- rrf_read rrf
      write . intercalate "\t" $
        (show <$> tuple9ToList tjs) <>
        [show (numerator r) <> "/" <> show (denominator r)]

  where numOfPrimes = 100
        maxStrLen   = 1024

tabulate :: Int
         -> String
         -> (Int -> (String -> IO ()) -> IO ())
         -> IO ()
tabulate tjMax name compute = do
  withFile filename WriteMode $ \ h ->
    compute tjMax (hPutStrLn h)
  where filename = "dist/rrf" <> name <> "-tj" <> show tjMax <> ".txt"

runReadP :: ReadP a -> String -> Maybe a
runReadP p s =
  case P.readP_to_S p s of
    [(x, _)] -> Just x
    _        -> Nothing

p_rrf :: ReadP Rational
p_rrf = combine <$> prefactor <*> radical
  where
    combine c r = signum c * r * c ^ (2 :: Int)
    prefactor = (*) <$> p_sign <*> (p_optionalParens p_frac <++ pure 1)
    times     = p_opt(P.char '*')
    radical   = times *> P.string "sqrt" *> p_parens p_frac <++ pure 1

p_opt :: ReadP a -> ReadP (Maybe a)
p_opt p = Just <$> p <++ pure Nothing

p_parens :: ReadP a -> ReadP a
p_parens p = P.char '(' *> p <* P.char ')'

p_optionalParens :: ReadP a -> ReadP a
p_optionalParens p = p_parens p <++ p

p_sign :: Num a => ReadP a
p_sign = (P.char '-' *> pure (-1) <++ pure 1)

p_frac :: ReadP Rational
p_frac = (%) <$> p_int <*> (P.char '/' *> p_int <++ pure 1)

p_int :: ReadP Integer
p_int = read <$> P.munch1 isDigit

withInt :: Int -> (Ptr CInt -> IO b) -> IO b
withInt = with . fromIntegral

-- | Be sure to specify enough primes or the library will crash your program.
initRRFLib :: Int -- ^ Number of prime numbers to precalculate.
           -> IO ()
initRRFLib primes =
  withInt 1 $ \ pquiet ->
  withInt primes $ \ pprimes ->
    f_init_rrf pquiet pprimes

data RRF = RRF !(Ptr (Ptr ())) !CStringLen

-- | We can't free 'RRF' because the Fortran program doesn't provide us with a
--   destructor! :(  Instead, try to reuse it for as long as you can to avoid
--   causing a massive memory leak.
rrf_new :: Int -- ^ maximum length of the string buffer used to marshal the RRF
        -> IO RRF
rrf_new maxStrLen =
  RRF <$> new nullPtr
      <*> newCStringLen (take maxStrLen (repeat '\NUL'))

rrf_ninej :: RRF -> (Int, Int, Int, Int, Int, Int, Int, Int, Int) -> IO ()
rrf_ninej (RRF pk _) (a, b, c, d, e, f, g, h, i) =
  withInt a $ \ pa ->
  withInt b $ \ pb ->
  withInt c $ \ pc ->
  withInt d $ \ pd ->
  withInt e $ \ pe ->
  withInt f $ \ pf ->
  withInt g $ \ pg ->
  withInt h $ \ ph ->
  withInt i $ \ pi ->
  withInt 2 $ \ px ->
    f_ninej pa pb pc pd pe pf pg ph pi px pk

rrf_show :: RRF -> IO String
rrf_show (RRF pk pStrLen@(pStr, maxLen)) =
  -- we don't know if f_char4i will overwrite maxLen-th element
  -- so just play it safe here
  withInt 1 $ \ pm1 ->
  withInt (maxLen - 2) $ \ pm2 -> do
    f_char4i pk pStr pm1 pm2
    trimEnd <$> peekCStringLen pStrLen
  where trimEnd = reverse . dropWhile (`elem` " \NUL") . reverse

rrf_read :: RRF -> IO Rational
rrf_read rrf = do
  s <- rrf_show rrf
  case runReadP (p_rrf <* P.eof) s of
    Nothing -> error ("rrf_get: cannot parse: " <> s)
    Just r  -> pure r
