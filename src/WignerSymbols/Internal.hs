-- | Stability : unstable
--
-- /Functions in this module are subject to change without notice./
--
{-# LANGUAGE BangPatterns, CPP #-}
module WignerSymbols.Internal where
#if !MIN_VERSION_base(4, 8, 0)
import Control.Applicative (pure)
#endif
import Control.Monad (guard)
import Data.List (sort)
import Data.Foldable (foldl')
import Data.Ratio ((%), numerator, denominator)

------------------------------------------------------------------------------

-- | Represents a mathematical expression of the form:
--
-- @
-- s √(n / d)
-- @
--
-- where
--
-- * @s@ is a sign (@+@, @-@, or @0@),
-- * @n@ is a nonnegative numerator, and
-- * @d@ is a positive denominator.
--
newtype SignedSqrtRational
  = SignedSqrtRational Rational
  deriving (Eq, Ord)

-- | / /
instance Read SignedSqrtRational where
  readsPrec p =
    readParen (p >= 11) $ \ s1 -> do
      ("ssr_new", s2) <- lex s1
      (x, s3) <- readsPrec 11 s2
      pure (ssr_new x, s3)

-- | / /
instance Show SignedSqrtRational where
  showsPrec p x =
    showParen (p >= 11) $
      showString "ssr_new " .
      showsPrec 11 (ssr_split x)

-- | Construct a 'SignedSqrtRational' equal to @c √r@.
{-# INLINE ssr_new #-}
ssr_new :: (Integer, Rational)          -- ^ @(c, r)@ / /
        -> SignedSqrtRational
ssr_new (c, r) = SignedSqrtRational (signum c % 1 * r * (c ^ (2 :: Int) % 1))

-- | Deconstruct a 'SignedSqrtRational'.
{-# INLINE ssr_split #-}
ssr_split :: SignedSqrtRational -> (Integer, Rational)
ssr_split (SignedSqrtRational x) = (signum (numerator x), abs x)

-- | Extract the sign of a 'SignedSqrtRational'.
{-# INLINE ssr_signum #-}
ssr_signum :: SignedSqrtRational -> Integer
ssr_signum (SignedSqrtRational r) = signum (numerator r)

-- | Extract the numerator of a 'SignedSqrtRational'.
{-# INLINE ssr_numerator #-}
ssr_numerator :: SignedSqrtRational -> Integer
ssr_numerator (SignedSqrtRational r) = abs (numerator r)

-- | Extract the denominator of a 'SignedSqrtRational'.
{-# INLINE ssr_denominator #-}
ssr_denominator :: SignedSqrtRational -> Integer
ssr_denominator (SignedSqrtRational r) = denominator r

-- | Approximate a 'SignedSqrtRational' as a floating-point number.
{-# INLINE ssr_approx #-}
ssr_approx :: Floating b => SignedSqrtRational -> b
ssr_approx x =
  case ssr_split x of
    (s, r) -> fromInteger s * sqrt (realToFrac r)

------------------------------------------------------------------------------

-- | Calculate a Clebsch-Gordan coefficient:
--
-- @
-- ⟨j1 j2 m1 m2|j1 j2 j12 m12⟩
-- @
{-# INLINABLE clebschGordan #-}
clebschGordan :: (Int, Int, Int, Int, Int, Int)
              -- ^ @(tj1, tm1, tj2, tm2, tj12, tm12)@ / /
              -> Double
clebschGordan = ssr_approx . clebschGordanSq

-- | Similar to 'clebschGordan' but exact.
{-# INLINABLE clebschGordanSq #-}
clebschGordanSq :: (Int, Int, Int, Int, Int, Int)
                -- ^ @(tj1, tm1, tj2, tm2, tj12, tm12)@ / /
                -> SignedSqrtRational
clebschGordanSq (tj1, tm1, tj2, tm2, tj12, tm12) =
  SignedSqrtRational (z * fromIntegral (tj12 + 1))
  where SignedSqrtRational z = wigner3jSqRawC (tj1, tm1, tj2, tm2, tj12, -tm12)

-- | Used only as a reference, it implements the formula from Wikipedia,
--   which comes from equation (2.41) on page 172 of
--   /Quantum Mechanics: Foundations and Applications/ (1993)
--   by A. Bohm and M. Loewe (ISBN 0-387-95330-2).
--
-- (Note: 'clebschGordan' is /not/ implemented using this function.)
{-# INLINABLE clebschGordanSqSlow #-}
clebschGordanSqSlow :: (Int, Int, Int, Int, Int, Int) -> SignedSqrtRational
clebschGordanSqSlow (tj1, tm1, tj2, tm2, tj12, tm12)
  | selectionRuleSatisfied = ssr_new (numerator (signum c), r * c ^ (2 :: Int))
  | otherwise              = ssr_new (0, 0)
  where

    selectionRuleSatisfied =
      tm1 + tm2 == tm12 &&
      abs tm1 <= tj1 &&
      abs tm2 <= tj2 &&
      abs tm12 <= tj12 &&
      (tj1 + tm1) `rem` 2 == 0 &&
      (tj2 + tm2) `rem` 2 == 0 &&
      triangleCondition (tj1, tj2, tj12)

    tkmin = -minimum [0, tj12 - tj2 + tm1, tj12 - tj1 - tm2]
    tkmax = minimum [tj1 + tj2 - tj12, tj1 - tm1, tj2 + tm2]

    facHalf n = factorial (n `quot` 2)

    c = sum [ toInteger (minusOnePow (tk `quot` 2))
              % ( facHalf tk
                * facHalf (tj1 + tj2 - tj12 - tk)
                * facHalf (tj1 - tm1 - tk)
                * facHalf (tj2 + tm2 - tk)
                * facHalf (tj12 - tj2 + tm1 + tk)
                * facHalf (tj12 - tj1 - tm2 + tk) )
            | tk <- [tkmin, tkmin + 2 .. tkmax] ]

    r = ( toInteger (tj12 + 1)
        * facHalf(tj12 + tj1 - tj2)
        * facHalf(tj12 - tj1 + tj2)
        * facHalf(tj1 + tj2 - tj12)
        * facHalf(tj12 + tm12)
        * facHalf(tj12 - tm12)
        * facHalf(tj1 - tm1)
        * facHalf(tj1 + tm1)
        * facHalf(tj2 - tm2)
        * facHalf(tj2 + tm2)
        ) % facHalf(tj1 + tj2 + tj12 + 2)

-- | Calculate a Wigner 3-j symbol:
--
-- @
-- ⎛j1 j2 j3⎞
-- ⎝m1 m2 m3⎠
-- @
{-# INLINABLE wigner3j #-}
wigner3j :: (Int, Int, Int, Int, Int, Int)
         -- ^ @(tj1, tm1, tj2, tm2, tj3, tm3)@ / /
         -> Double
wigner3j = ssr_approx . wigner3jSq

-- | Similar to 'wigner3j' but exact.
{-# INLINABLE wigner3jSq #-}
wigner3jSq :: (Int, Int, Int, Int, Int, Int)
           -- ^ @(tj1, tm1, tj2, tm2, tj3, tm3)@ / /
           -> SignedSqrtRational
wigner3jSq (tj1, tm1, tj2, tm2, tj3, tm3) = SignedSqrtRational (s * z)
  where s = fromIntegral (minusOnePow ((tj1 - tj2 - tm3) `quot` 2))
        SignedSqrtRational z = wigner3jSqRawC (tj1, tm1, tj2, tm2, tj3, tm3)

-- | Calculate the Wigner 3-j symbol times @(−1) ^ (j1 − j2 − m3)@.
{-# INLINABLE wigner3jSqRawC #-}
wigner3jSqRawC :: (Int, Int, Int, Int, Int, Int)
               -- ^ @(tj1, tm1, tj2, tm2, tj3, tm3)@ / /
               -> SignedSqrtRational
wigner3jSqRawC tjms@(tj1, tm1, tj2, tm2, tj3, tm3) =
  if tm1 + tm2 + tm3 == 0 &&
     abs tm1 <= tj1 &&
     abs tm2 <= tj2 &&
     abs tm3 <= tj3 &&
     jmr1 == 0 &&
     jmr2 == 0 &&
     triangleCondition (tj1, tj2, tj3)
  then wigner3jSqRaw (jm1, jm2) tjms
  else ssr_new (0, 0)
  where
    (!jm1, !jmr1) = (tj1 + tm1) `quotRem` 2
    (!jm2, !jmr2) = (tj2 + tm2) `quotRem` 2

-- | Calculate the Wigner 3-j symbol times @(−1) ^ (j1 − j2 − m3)@.
--   The selection rules are not checked.
{-# INLINE wigner3jSqRaw #-}
wigner3jSqRaw :: (Int, Int)
              -- ^ @(j1 + m1, j2 + m2)@ / /
              -> (Int, Int, Int, Int, Int, Int)
              -- ^ @(tj1, tm1, tj2, tm2, tj3, tm3)@ / /
              -> SignedSqrtRational
wigner3jSqRaw (jm1, jm2) (tj1, tm1, tj2, tm2, tj3, tm3) = ssr_new (z2, z1)
  where

    !z1 = (binomial tj1 jjj1 * binomial tj2 jjj2 * binomial tj3 jjj3)
        % (binomial tj1 jm1 * binomial tj2 jm2 * binomial tj3 jm3)
        * triangularFactorRaw jjj (jjj1, jjj2, jjj3)

    !z2 = if kmin > kmax
          then 0
          else let !c0 = toInteger (minusOnePow kmin)
                       * binomial jjj2 kmin
                       * binomial jjj1 (jsm1 - kmin)
                       * binomial jjj3 (jm2 - kmin)
               in fst (foldl' f (c0, c0) [succ kmin .. kmax])

    f (s, c) k = (s', -c')
      where
        !c' = c
            * toInteger (jjj2 - k + 1) `quot` toInteger k
            * toInteger (jsm1 - k + 1) `quot` toInteger (jjj1 - (jsm1 - k))
            * toInteger (jm2  - k + 1) `quot` toInteger (jjj3 - (jm2  - k))
        !s' = s - c'

    !kmin = maximum [0, tj1 - tj3 + tm2, tj2 - tj3 - tm1] `quot` 2
    !kmax = minimum [jjj2, jsm1, jm2]

    !jjj1 = (tj1 - tj2 + tj3) `quot` 2
    !jjj2 = (tj2 - tj3 + tj1) `quot` 2
    !jjj3 = (tj3 - tj1 + tj2) `quot` 2
    !jjj  = (tj1 + tj2 + tj3) `quot` 2 + 1

    !jsm1 = (tj1 - tm1) `quot` 2
    !jm3  = (tj3 + tm3) `quot` 2

-- | Calculate a Wigner 6-j symbol:
--
-- @
-- ⎧j11 j12 j13⎫
-- ⎩j21 j22 j23⎭
-- @
{-# INLINABLE wigner6j #-}
wigner6j :: (Int, Int, Int, Int, Int, Int)
         -- ^ @(tj11, tj12, tj13, tj21, tj22, tj23)@ / /
         -> Double
wigner6j = ssr_approx . wigner6jSq

-- | Similar to 'wigner6j' but exact.
{-# INLINABLE wigner6jSq #-}
wigner6jSq :: (Int, Int, Int, Int, Int, Int)
           -- ^ @(tj11, tj12, tj13, tj21, tj22, tj23)@ / /
           -> SignedSqrtRational
wigner6jSq tjs@(tja, tjb, tjc, tjd, tje, tjf) =
  if triangleCondition (tja, tjb, tjc) &&
     triangleCondition (tja, tje, tjf) &&
     triangleCondition (tjd, tjb, tjf) &&
     triangleCondition (tjd, tje, tjc)
  then wigner6jSqRaw tjs
  else ssr_new (0, 0)

-- | Calculate the Wigner 6-j symbol.  The selection rules are not checked.
{-# INLINE wigner6jSqRaw #-}
wigner6jSqRaw :: (Int, Int, Int, Int, Int, Int)
              -- ^ @(tj11, tj12, tj13, tj21, tj22, tj23)@ / /
              -> SignedSqrtRational
wigner6jSqRaw (tja, tjb, tjc, tjd, tje, tjf) = ssr_new (z2, z1)
  where

    !z1 = triangularFactor (tja, tje, tjf)
        * triangularFactor (tjd, tjb, tjf)
        * triangularFactor (tjd, tje, tjc)
        / triangularFactor (tja, tjb, tjc)

    !z2 = tetrahedralSum (tja, tje, tjf, tjd, tjb, tjc)

-- | Calculate a Wigner 9-j symbol:
--
-- @
-- ⎧j11 j12 j13⎫
-- ⎨j21 j22 j23⎬
-- ⎩j31 j32 j33⎭
-- @
{-# INLINABLE wigner9j #-}
wigner9j :: (Int, Int, Int, Int, Int, Int, Int, Int, Int)
         -- ^ @(tj11, tj12, tj13, tj21, tj22, tj23, tj31, tj32, tj33)@ / /
         -> Double
wigner9j = ssr_approx . wigner9jSq

-- | Similar to 'wigner9j' but exact.
{-# INLINABLE wigner9jSq #-}
wigner9jSq :: (Int, Int, Int, Int, Int, Int, Int, Int, Int)
           -- ^ @(tj11, tj12, tj13, tj21, tj22, tj23, tj31, tj32, tj33)@ / /
           -> SignedSqrtRational
wigner9jSq tjs@(tja, tjb, tjc, tjd, tje, tjf, tjg, tjh, tji) =
  if triangleCondition (tja, tjb, tjc) &&
     triangleCondition (tjd, tje, tjf) &&
     triangleCondition (tjg, tjh, tji) &&
     triangleCondition (tja, tjd, tjg) &&
     triangleCondition (tjb, tje, tjh) &&
     triangleCondition (tjc, tjf, tji)
  then wigner9jSqRaw tjs
  else ssr_new (0, 0)

-- | Calculate the Wigner 9-j symbol.  The selection rules are not checked.
{-# INLINE wigner9jSqRaw #-}
wigner9jSqRaw :: (Int, Int, Int, Int, Int, Int, Int, Int, Int)
              -- ^ @(tj11, tj12, tj13, tj21, tj22, tj23, tj31, tj32, tj33)@ / /
              -> SignedSqrtRational
wigner9jSqRaw (tja, tjb, tjc, tjd, tje, tjf, tjg, tjh, tji) = ssr_new (z2, z1)
  where

    !z1 =
      triangularFactor (tja, tjb, tjc) *
      triangularFactor (tjd, tje, tjf) *
      triangularFactor (tjg, tjh, tji) *
      triangularFactor (tja, tjd, tjg) *
      triangularFactor (tjb, tje, tjh) *
      triangularFactor (tjc, tjf, tji)

    !z2 =
      sum [ toInteger (minusOnePow tk * (tk + 1))
          * tetrahedralSum (tja, tjb, tjc, tjf, tji, tk)
          * tetrahedralSum (tjf, tjd, tje, tjh, tjb, tk)
          * tetrahedralSum (tjh, tji, tjg, tja, tjd, tk)
          | tk <- [tkmin, tkmin + 2 .. tkmax] ]

    !tkmin =
      maximum
      [ abs (tjh - tjd)
      , abs (tjb - tjf)
      , abs (tja - tji) ]

    !tkmax =
      minimum
      [ tjh + tjd
      , tjb + tjf
      , tja + tji ]

------------------------------------------------------------------------------

-- | Calculate the factorial @n!@.
{-# INLINABLE factorial #-}
factorial :: Int                        -- ^ @n@ / /
          -> Integer
factorial = go 1
  where go r n
          | n <= 1    = r
          | otherwise =
              let !r' = r * toInteger n
                  !n' = n - 1
              in go r' n'

-- | Calculate the falling factorial, i.e. the product of the integers @[n, k)@.
{-# INLINABLE fallingFactorial #-}
fallingFactorial :: Int                 -- ^ @n@ / /
                 -> Int                 -- ^ @k@ / /
                 -> Integer
fallingFactorial = go 1 1
  where go r i n k
          | i > k     = r
          | otherwise =
              let !r' = r * toInteger n
                  !i' = i + 1
                  !n' = n - 1
              in go r' i' n' k

-- | Calculate the binomial coefficient @C(n, k)@.
{-# INLINABLE binomial #-}
binomial :: Int                         -- ^ @n@ / /
         -> Int                         -- ^ @k@ / /
         -> Integer
binomial = go 1 1
  where go r i n k
          | i > k     = r
          | otherwise =
              let !r' = r * toInteger n `quot` toInteger i
                  !i' = i + 1
                  !n' = n - 1
              in go r' i' n' k

-- | Calculate @(−1) ^ n@.
{-# INLINABLE minusOnePow #-}
minusOnePow :: Int                      -- ^ @n@ / /
            -> Int
minusOnePow n = 1 - n `mod` 2 * 2

-- | Check @|j1 − j2| ≤ j3 ≤ j1 + j2@ and @j1 + j2 + j3 ∈ ℤ@.
{-# INLINABLE triangleCondition #-}
triangleCondition :: (Int, Int, Int)    -- ^ @(tj1, tj2, tj3)@.
                  -> Bool
triangleCondition (a, b, c) = d >= 0 && d `rem` 2 == 0 && c - abs (a - b) >= 0
  where !d = a + b - c

-- | Calculate the triangular factor:
--
-- @
-- Δ(j1, j2, j3) = (−j1 + j2 _ j3)! (j1 − j2 + j3)! (j1 + j2 − j3)!
--               / (j1 + j2 + j3 + 1)!
-- @
--
{-# INLINABLE triangularFactor #-}
triangularFactor :: (Int, Int, Int)     -- ^ @(tj1, tj2, tj3)@ / /
                 -> Rational
triangularFactor (tja, tjb, tjc) = triangularFactorRaw jjj (jjja, jjjb, jjjc)
  where !jjja = (tjc - tja + tjb) `quot` 2
        !jjjb = (tja - tjb + tjc) `quot` 2
        !jjjc = (tjb - tjc + tja) `quot` 2
        !jjj  = (tja + tjb + tjc) `quot` 2 + 1

-- | Calculate @ja! jb! jc! / jd!@.
{-# INLINABLE triangularFactorRaw #-}
triangularFactorRaw :: Int              -- ^ @jd@ / /
                    -> (Int, Int, Int)  -- ^ @(ja, jb, jc)@ / /
                    -> Rational
triangularFactorRaw jjj (jjja, jjjb, jjjc) =
  factorial jjju * factorial jjjv % fallingFactorial jjj (jjj - jjjw)
  where [!jjju, !jjjv, !jjjw] = sort [jjja, jjjb, jjjc]

-- | Calculate the symbol in the paper by L. Wei that is enclosed in square
--   brackets:
--
-- @
-- ⎡j11 j12 j13⎤
-- ⎣j21 j22 j23⎦
-- @
--
-- This is essentially a Wigner 6-j symbol without the triangular factors,
-- although the ordering of the arguments is a bit funky here.
--
{-# INLINABLE tetrahedralSum #-}
tetrahedralSum :: (Int, Int, Int, Int, Int, Int)
               -- ^ @(tj11, tj12, tj13, tj21, tj22, tj23)@ / /
               -> Integer
tetrahedralSum (tja, tje, tjf, tjd, tjb, tjc) =
  sum [ toInteger (minusOnePow k)
      * binomial (k + 1) (k - jabc)
      * binomial jjja (k - jaef)
      * binomial jjjb (k - jdbf)
      * binomial jjjc (k - jdec)
      | k <- [kmin .. kmax] ]
  where

    !jjja = (tjc - tja + tjb) `quot` 2
    !jjjb = (tja - tjb + tjc) `quot` 2
    !jjjc = (tjb - tjc + tja) `quot` 2

    !jabc = (tja + tjb + tjc) `quot` 2
    !jaef = (tja + tje + tjf) `quot` 2
    !jdbf = (tjd + tjb + tjf) `quot` 2
    !jdec = (tjd + tje + tjc) `quot` 2

    !kmin = maximum [jabc, jdec, jdbf, jaef]

    !kmax = minimum [ tja + tjd + tjb + tje
                    , tjb + tje + tjc + tjf
                    , tja + tjd + tjc + tjf ] `quot` 2

------------------------------------------------------------------------------

-- | Get all angular momenta that satisfy the triangle condition with the
--   given pair of angular momenta, up to a maximum of @jmax@.
{-# INLINABLE getTriangularTjs #-}
getTriangularTjs :: Int                 -- ^ @tjmax@ / /
                 -> (Int, Int)          -- ^ @(tj1, tj2)@ / /
                 -> [Int]
getTriangularTjs tjMax (tja, tjb) = [tjmin, tjmin + 2 .. tjmax]
  where tjmin = abs (tja - tjb)
        tjmax = min tjMax (tja + tjb)

-- | Get all angular momenta that satisfy the triangle condition with each of
--   the given two pairs of angular momenta, up to a maximum of @jmax@.
{-# INLINABLE getBitriangularTjs #-}
getBitriangularTjs :: Int   -- ^ @tjmax@ / /
                   -> ((Int, Int), (Int, Int))
                   -- ^ @((tj11, tj12), (tj21, tj22))@ / /
                   -> [Int]
getBitriangularTjs tjMax ((tja, tjb), (tjc, tjd)) = [tjmin, tjmin + 2 .. tjmax]
  where tjmin = max (abs (tja - tjb)) (abs (tjc - tjd))
        tjmax = minimum [tjMax, tja + tjb, tjc + tjd]

-- | Get all projection quantum numbers that lie within the multiplet of @j@.
{-# INLINABLE getTms #-}
getTms :: Int                           -- ^ @tj@ / /
       -> [Int]
getTms tj = [-tj, -tj + 2 .. tj]

-- | Get all possible arguments of the Wigner 3-j symbol that satisfy the
--   selection rules up to a maximum of @jmax@.
{-# INLINABLE get3tjms #-}
get3tjms :: Int                         -- ^ @tjmax@ / /
         -> [(Int, Int, Int, Int, Int, Int)]
get3tjms tjMax = do
  tj1 <- [0 .. tjMax]
  tj2 <- [0 .. tjMax]
  tj3 <- getTriangularTjs tjMax (tj1, tj2)
  tm1 <- getTms tj1
  tm2 <- getTms tj2
  let tm3 = -(tm1 + tm2)
  guard (abs tm3 <= tj3)
  pure (tj1, tm1, tj2, tm2, tj3, tm3)

-- | Get all possible arguments of the Wigner 6-j symbol that satisfy the
--   selection rules up to a maximum of @jmax@.
{-# INLINABLE get6tjs #-}
get6tjs :: Int -> [(Int, Int, Int, Int, Int, Int)]
get6tjs tjMax = do
  tja <- [0 .. tjMax]
  tjb <- [0 .. tjMax]
  tjc <- getTriangularTjs tjMax (tja, tjb)
  tjd <- [0 .. tjMax]
  tje <- getTriangularTjs tjMax (tjd, tjc)
  tjf <- getBitriangularTjs tjMax ((tja, tje), (tjd, tjb))
  pure (tja, tjb, tjc, tjd, tje, tjf)

-- | Get all possible arguments of the Wigner 9-j symbol that satisfy the
--   selection rules up to a maximum of @jmax@.
{-# INLINABLE get9tjs #-}
get9tjs :: Int -> [(Int, Int, Int, Int, Int, Int, Int, Int, Int)]
get9tjs tjMax = do
  tja <- [0 .. tjMax]
  tjb <- [0 .. tjMax]
  tjc <- getTriangularTjs tjMax (tja, tjb)
  tjd <- [0 .. tjMax]
  tje <- [0 .. tjMax]
  tjf <- getTriangularTjs tjMax (tjd, tje)
  tjg <- getTriangularTjs tjMax (tja, tjd)
  tjh <- getTriangularTjs tjMax (tjb, tje)
  tji <- getBitriangularTjs tjMax ((tjg, tjh), (tjc, tjf))
  pure (tja, tjb, tjc, tjd, tje, tjf, tjg, tjh, tji)

-- | Convert a 6-tuple into a list.
{-# INLINABLE tuple6ToList #-}
tuple6ToList :: (a, a, a, a, a, a) -> [a]
tuple6ToList (a, b, c, d, e, f) = [a, b, c, d, e, f]

-- | Convert a 9-tuple into a list.
{-# INLINABLE tuple9ToList #-}
tuple9ToList :: (a, a, a, a, a, a, a, a, a) -> [a]
tuple9ToList (a, b, c, d, e, f, g, h, i) = [a, b, c, d, e, f, g, h, i]
