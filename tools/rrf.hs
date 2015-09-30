-- This code tabulates Wigner symbols using the Root Rational Fraction (RRF)
-- program by Anthony Stone and Charles Wood.  It is used here for
-- verification purposes.
--
-- http://www-stone.ch.cam.ac.uk/documentation/rrf
-- http://www-stone.ch.cam.ac.uk/pub/rrf-4.0.tgz
module Main (main) where
import Data.Char (isDigit)
import Data.Ratio ((%), denominator, numerator)
import Data.Foldable (for_)
import Foreign
import Foreign.C
import Prelude hiding (pi)
import Text.ParserCombinators.ReadP (ReadP)
import qualified Text.ParserCombinators.ReadP as P
import WignerSymbols.Internal

infixl 3 <++
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
  let numOfPrimes = 100
      maxStrLen   = 1024
  initRRFLibrary numOfPrimes
  rrf <- rrf_new maxStrLen
  for_ (get9tjs 3) $ \ tjs -> do
    s <- rrf_ninej rrf tjs
    let Just r = runReadP (p_rrf <* P.eof) s
    putStrLn (show (numerator r) ++ "/" ++ show (denominator r))

fromJust (Just x) = x

runReadP :: ReadP a -> String -> Maybe a
runReadP p s =
  case P.readP_to_S p s of
    [(x, _)] -> Just x
    _        -> Nothing

p_rrf = combine <$> prefactor <*> radical
  where
    combine c r = signum c * r * c ^ (2 :: Int)
    prefactor = (*) <$> p_sign <*> p_optionalParens p_frac
    radical   = P.string "*sqrt" *> p_parens p_frac <++ pure 1

p_parens p = P.char '(' *> p <* P.char ')'

p_optionalParens p = p_parens p <++ p

p_sign :: Num a => ReadP a
p_sign = (P.char '-' *> pure (-1) <++ pure 1)

p_frac = (%) <$> p_int <*> (P.char '/' *> p_int <++ pure 1)

p_int :: ReadP Integer
p_int = read <$> P.munch1 isDigit

withInt :: Int -> (Ptr CInt -> IO b) -> IO b
withInt = with . fromIntegral

data RRF =
  RRF
  { rrf_ptr  :: !Ptr (Ptr ())
  , rrf_buf :: !CStringLen
  }

-- | We can't free 'RRF' because the Fortran program doesn't provide us with a
--   destructor! :(  Instead, try to reuse it for as long as you can to avoid
--   causing a massive memory leak.
rrf_new :: Int -- ^ maximum length of the string buffer used to marshal the RRF
        -> IO RRF
rrf_new maxStrLen =
  RRF <$> new nullPtr
      <*> newCStringLen (take maxStrLen (repeat 0))
      <*> newCStringLen (take maxStrLen (repeat 0))

-- | Be sure to specify enough primes or the library will crash your program.
initRRFLibrary :: Int -- ^ Number of prime numbers to precalculate.
               -> IO ()
initRRFLibrary primes =
  withInt 1 $ \ pquiet ->
  withInt primes $ \ pprimes ->
    f_init_rrf pquiet pprimes

rrf_ninej :: RRF -> (Int, Int, Int, Int, Int, Int, Int, Int, Int) -> IO ()
rrf_ninej rrf@(RRF pk _) (a, b, c, d, e, f, g, h, i) =
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
    char4i rrf

char4i :: RRF -> IO String
char4i (RRF pk pStrLen@(pStr, maxLen)) =
  -- we don't know if f_char4i will overwrite maxLen-th element
  -- so just play it safe here
  withInt (maxLen - 2) $ \ pm2 ->
  withInt 1 $ \ pm1 -> do
    f_char4i pk pStr pm1 pm2
    trimEnd <$> peekCStringLen pStrLen
  where trimEnd = reverse . dropWhile (`elem` " \NUL") . reverse
