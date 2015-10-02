{-# LANGUAGE BangPatterns, CPP #-}
module WignerSymbols.Internal where
#if !MIN_VERSION_base(4, 8, 0)
import Control.Applicative (pure)
#endif
import Control.Monad (guard)
import Data.List (sort)
import Data.Foldable (foldl')
import Data.Ratio (Ratio, (%), numerator, denominator)

-- | Represents a mathematical expression of the form:
--
-- @
--     s × √(n / d)
-- @
--
-- where
--
-- * @s@ is a sign (@+@, @-@, or @0@),
-- * @n@ is a nonnegative numerator, and
-- * @d@ is a positive denominator.

newtype SignedSqrtRatio a = SignedSqrtRatio (Ratio a)
                          deriving (Eq, Read, Show)

type SignedSqrtRational = SignedSqrtRatio Integer

{-# INLINABLE ssr_split #-}
ssr_split :: (Integral a, Num b) => SignedSqrtRatio a -> (b, Ratio a)
ssr_split (SignedSqrtRatio x) =
  (fromIntegral (signum (numerator x)), abs x)

{-# INLINABLE ssr_approx #-}
ssr_approx :: (Integral a, Floating b) => SignedSqrtRatio a -> b
ssr_approx x = case ssr_split x of
  (s, r) -> s * sqrt (realToFrac r)

{-# INLINABLE ssr_numerator #-}
ssr_numerator :: Integral a => SignedSqrtRatio a -> a
ssr_numerator (SignedSqrtRatio r) = abs (numerator r)

{-# INLINABLE ssr_denominator #-}
ssr_denominator :: Integral a => SignedSqrtRatio a -> a
ssr_denominator (SignedSqrtRatio r) = denominator r

{-# INLINABLE ssr_signum #-}
ssr_signum :: Integral a => SignedSqrtRatio a -> a
ssr_signum (SignedSqrtRatio r) = fromIntegral (signum (numerator r))

{-# INLINABLE factorial #-}
factorial :: Int -> Integer
factorial = go 1
  -- difference between toInteger first or toInteger per iteration seems negligible
  where go r n
          | n <= 1    = r
          | otherwise =
              let !r' = r * toInteger n
                  !n' = n - 1
              in go r' n'

{-# INLINABLE factorial_simple #-}
-- Difference between _simple and the other one seems negligible
factorial_simple :: Int -> Integer
factorial_simple n = product [1 .. toInteger n]

{-# INLINABLE binomial #-}
binomial :: (Fractional a, Integral b) => a -> b -> a
binomial = binomial_descending (\ x y -> x / fromIntegral y)

-- binomial_descending is about twice as slow as the binomialI
-- \ n k -> numerator (binomial_descending (%) n k)
{-# INLINABLE binomial_descending #-}
binomial_descending :: (Num a, Integral b, Num c) => (a -> b -> c) -> a -> b -> c
binomial_descending divide n k
  | k > 0     = divide n k * binomial_descending divide (n - 1) (k - 1)
  | otherwise = 1

{-# INLINABLE binomialI #-}
binomialI :: Int -> Int -> Integer
binomialI = go 1 1
  where go r i n k
          | i > k     = r
          | otherwise =
              let !r' = r * toInteger n `quot` toInteger i
                  !i' = i + 1
                  !n' = n - 1
              in go r' i' n' k

{-# INLINABLE fallingFactorialI #-}
fallingFactorialI :: Int -> Int -> Integer
fallingFactorialI = go 1 1
  where go r i n k
          | i > k     = r
          | otherwise =
              let !r' = r * toInteger n
                  !i' = i + 1
                  !n' = n - 1
              in go r' i' n' k

{-# INLINABLE minusOnePow #-}
minusOnePow :: Int -> Int
minusOnePow n = 1 - n `mod` 2 * 2

-- note: -fllvm makes this faster (~10%)

{-# INLINABLE triangleCondition #-}
-- Int -> Int -> Int -> Bool
triangleCondition :: (Num a, Ord a) => a -> a -> a -> Bool
triangleCondition a b c = abs (a - b) <= c && c <= a + b

{-# INLINABLE triangleConditionI #-}
--triangleConditionI :: Int -> Int -> Int -> Bool
triangleConditionI :: Integral a => a -> a -> a -> Bool
triangleConditionI a b c =
  c - abs (a - b) >= 0 &&
  a + b - c >= 0 &&
  (a + b - c) `rem` 2 == 0

-- | Calculate a Clebsch-Gordan coefficient.
{-# INLINABLE clebschGordan #-}
clebschGordan :: (Int, Int, Int, Int, Int, Int)
              -- ^ @(tj1, tm1, tj2, tm2, tj12, tm12)@.
              -> Double
clebschGordan = ssr_approx . clebschGordanSq

{-# INLINABLE clebschGordanSq #-}
clebschGordanSq :: (Int, Int, Int, Int, Int, Int)
                -- ^ @(tj1, tm1, tj2, tm2, tj12, tm12)@.
                -> SignedSqrtRational
clebschGordanSq (tj1, tm1, tj2, tm2, tj12, tm12) =
  SignedSqrtRatio (z * fromIntegral (tj12 + 1))
  where SignedSqrtRatio z = wigner3jSqRaw (tj1, tm1, tj2, tm2, tj12, -tm12)

-- | Calculate a Wigner 3-j symbol.
{-# INLINABLE wigner3j #-}
wigner3j :: (Int, Int, Int, Int, Int, Int)
         -- ^ @(tj1, tm1, tj2, tm2, tj3, tm3)@.
         -> Double
wigner3j = ssr_approx . wigner3jSq

{-# INLINABLE wigner3jSq #-}
wigner3jSq :: (Int, Int, Int, Int, Int, Int)
           -- ^ @(tj1, tm1, tj2, tm2, tj3, tm3)@.
           -> SignedSqrtRational
wigner3jSq (tj1, tm1, tj2, tm2, tj3, tm3) = SignedSqrtRatio (s * z)
  where s = fromIntegral (minusOnePow ((tj1 - tj2 - tm3) `quot` 2))
        SignedSqrtRatio z = wigner3jSqRaw (tj1, tm1, tj2, tm2, tj3, tm3)

-- | This uses the formula described in
-- Wei1999 doi:10.1016/S0010-4655(99)00232-5
-- http://meghnad.iucaa.ernet.in/~tarun/pprnt/compute/ClebADKL.pdf
{-# INLINABLE wigner3jSqRaw #-}
wigner3jSqRaw :: (Int, Int, Int, Int, Int, Int) -> SignedSqrtRational
wigner3jSqRaw (tj1, tm1, tj2, tm2, tj3, tm3)
  | satisfiesSelectionRule = SignedSqrtRatio z
  | otherwise              = SignedSqrtRatio 0
  where

    satisfiesSelectionRule =
      tm1 + tm2 + tm3 == 0 &&
      abs tm1 <= tj1 &&
      abs tm2 <= tj2 &&
      abs tm3 <= tj3 &&
      jmr1 == 0 &&
      jmr2 == 0 &&
      jmr3 == 0 &&
      triangleCondition tj1 tj2 tj3

    z :: Rational
    z = fromIntegral (signum z3) * z1 * z2 * fromIntegral (z3 ^ (2 :: Int))

    -- merging z1 and z2 makes it slower
    z1 :: Rational
    z1 = factorial jjja * factorial jjjb %
         -- using fallingFactorialI gives ~10% improvement for large j
         fallingFactorialI (succ jjj) (succ jjj - jjjc)
      where [jjja, jjjb, jjjc] = sort [jjj1, jjj2, jjj3]

    z2 :: Rational
    z2 =
      -- splitting z2 will make it slower
      (binomialI tj1 jjj1 * binomialI tj2 jjj2 * binomialI tj3 jjj3) %
      (binomialI tj1 jm1 * binomialI tj2 jm2 * binomialI tj3 jm3)

    z3 :: Integer
    z3 | kmin > kmax = 0
       | otherwise   =
           let !c0 = toInteger (minusOnePow kmin)
                   * binomialI jjj2 kmin
                   * binomialI jjj1 (jsm1 - kmin)
                   * binomialI jjj3 (jm2 - kmin)
           in fst (foldl' f (c0, c0) [succ kmin .. kmax])

    f (s, c) k = (s', -c')
      where
        !c' = c
            * toInteger (jjj2 - k + 1) `quot` toInteger k
            * toInteger (jsm1 - k + 1) `quot` toInteger (jjj1 - (jsm1 - k))
            * toInteger (jm2 - k + 1) `quot` toInteger (jjj3 - (jm2 - k))
        !s' = s - c'

    kmin = maximum [0, tj1 - tj3 + tm2, tj2 - tj3 - tm1] `quot` 2
    kmax = minimum [jjj2, jsm1, jm2]

    jjj1 = (tj1 - tj2 + tj3) `quot` 2
    jjj2 = (tj2 - tj3 + tj1) `quot` 2
    jjj3 = (tj3 - tj1 + tj2) `quot` 2
    jjj  = (tj1 + tj2 + tj3) `quot` 2

    (jm1, jmr1) = (tj1 + tm1) `quotRem` 2
    (jm2, jmr2) = (tj2 + tm2) `quotRem` 2
    (jm3, jmr3) = (tj3 + tm3) `quotRem` 2

    jsm1 = (tj1 - tm1) `quot` 2

{-# INLINABLE clebschGordanSqSlow #-}
clebschGordanSqSlow :: (Int, Int, Int, Int, Int, Int) -> SignedSqrtRational
clebschGordanSqSlow (tj1, tm1, tj2, tm2, tj12, tm12)
  | conservationLaw &&
    triangleCondition tj1 tj2 tj12 &&
    jParity = SignedSqrtRatio (sign * surd)
  | otherwise = SignedSqrtRatio 0
  where

    conservationLaw = tm1 + tm2 == tm12

    jParity =
      (tj1 + tj2 + tj12) `rem` 2 == 0 &&
      (tj1 + tm1) `rem` 2 == 0 &&
      (tj2 + tm2) `rem` 2 == 0 &&
      (tj12 + tm12) `rem` 2 == 0

    tkmin = -minimum [0, tj12 - tj2 + tm1, tj12 - tj1 - tm2]
    tkmax = minimum [tj1 + tj2 - tj12, tj1 - tm1, tj2 + tm2]

    facHalf n = factorial (n `quot` 2)

    -- r = sum [minusOnePow (tk `quot` 2) %
    --          (facHalf(tk) *
    --           facHalf(tj1 + tj2 - tj12 - tk) *
    --           facHalf(tj1 - tm1 - tk) *
    --           facHalf(tj2 + tm2 - tk) *
    --           facHalf(tj12 - tj2 + tm1 + tk) *
    --           facHalf(tj12 - tj1 - tm2 + tk)
    --          ) | tk <- [tkmin, tkmin + 2 .. tkmax]]

    -- here we calculate the summation; we accumulate the factorials to avoid
    -- repeating the same calculations; this gives a ~10% speed increase

    r | tkmin > tkmax = 0
      | otherwise     = s
      where (s, _, _, _, _, _, _, _) =
              foldl' adder (ic, ic, ic1, ic2, ic3, ic4, ic5, ic6)
                     [tkmin + 2 .. tkmax]

    ic = toInteger (minusOnePow (tkmin `quot` 2)) %
         product (facHalf <$> [ic1, ic2, ic3, ic4, ic5, ic6])

    ic1 = tkmin
    ic2 = tj1 + tj2 - tj12 - tkmin
    ic3 = tj1 - tm1 - tkmin
    ic4 = tj2 + tm2 - tkmin
    ic5 = tj12 - tj2 + tm1 + tkmin
    ic6 = tj12 - tj1 - tm2 + tkmin

    adder (s, c, c1, c2, c3, c4, c5, c6) _ =
      ( s + c'
      , c'
      , c1 + 2
      , c2 - 2
      , c3 - 2
      , c4 - 2
      , c5 + 2
      , c6 + 2
      )
      where c' = c * negate (toInteger (c2 * c3 * c4) %
                             toInteger ((c1 + 2) * (c5 + 2) * (c6 + 2)))

    sign = fromIntegral (numerator (signum r))
    surd | r == 0    = 0
         | otherwise = q

    q = (
            (
                toInteger (tj12 + 1) *
                facHalf(tj12 + tj1 - tj2) *
                facHalf(tj12 - tj1 + tj2) *
                facHalf(tj1 + tj2 - tj12) *
                facHalf(tj12 + tm12) *
                facHalf(tj12 - tm12) *
                facHalf(tj1 - tm1) *
                facHalf(tj1 + tm1) *
                facHalf(tj2 - tm2) *
                facHalf(tj2 + tm2)
            ) %
            facHalf(tj1 + tj2 + tj12 + 2)
        ) * r ^ (2 :: Int)

-- | Calculate a Wigner 6-j symbol.
{-# INLINABLE wigner6j #-}
wigner6j :: (Int, Int, Int, Int, Int, Int)
         -- ^ @(tj11, tj12, tj13, tj21, tj22, tj23)@.
         -> Double
wigner6j = ssr_approx . wigner6jSq

{-# INLINABLE wigner6jSq #-}
wigner6jSq :: (Int, Int, Int, Int, Int, Int) -> SignedSqrtRational
wigner6jSq (tja, tjb, tjc, tjd, tje, tjf)
  | satisfiesSelectionRule = SignedSqrtRatio z
  | otherwise              = SignedSqrtRatio 0
  where

    satisfiesSelectionRule =
      triangleConditionI tja tjb tjc &&
      triangleConditionI tja tje tjf &&
      triangleConditionI tjd tjb tjf &&
      triangleConditionI tjd tje tjc

    z = fromInteger (signum z2) * z1 * fromInteger (z2 ^ (2 :: Int))

    z1 =
      triangleFactor (tja, tje, tjf) *
      triangleFactor (tjd, tjb, tjf) *
      triangleFactor (tjd, tje, tjc) /
      triangleFactor (tja, tjb, tjc)

    z2 = tetrahedralSum (tja, tje, tjf, tjd, tjb, tjc)

-- | Calculate a Wigner 9-j symbol.
{-# INLINABLE wigner9j #-}
wigner9j :: (Int, Int, Int, Int, Int, Int, Int, Int, Int)
         -- ^ @(tj11, tj12, tj13, tj21, tj22, tj23, tj31, tj32, tj33)@.
         -> Double
wigner9j = ssr_approx . wigner9jSq

{-# INLINABLE wigner9jSq #-}
wigner9jSq :: (Int, Int, Int, Int, Int, Int, Int, Int, Int)
           -> SignedSqrtRational
wigner9jSq tjs@(tja, tjb, tjc, tjd, tje, tjf, tjg, tjh, tji) =
  SignedSqrtRatio $
  if triangleConditionI tja tjb tjc &&
     triangleConditionI tjd tje tjf &&
     triangleConditionI tjg tjh tji &&
     triangleConditionI tja tjd tjg &&
     triangleConditionI tjb tje tjh &&
     triangleConditionI tjc tjf tji
  then wigner9jSqRaw tjs
  else 0

{-# INLINABLE wigner9jSqRaw #-}
wigner9jSqRaw :: (Int, Int, Int, Int, Int, Int, Int, Int, Int) -> Rational
wigner9jSqRaw (tja, tjb, tjc, tjd, tje, tjf, tjg, tjh, tji) = z
  where

    !z = fromInteger (signum z2) * z1 * fromInteger (z2 ^ (2 :: Int))

    !z1 =
      triangleFactor (tja, tjb, tjc) *
      triangleFactor (tjd, tje, tjf) *
      triangleFactor (tjg, tjh, tji) *
      triangleFactor (tja, tjd, tjg) *
      triangleFactor (tjb, tje, tjh) *
      triangleFactor (tjc, tjf, tji)

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

{-# INLINABLE triangleFactor #-}
triangleFactor :: (Int, Int, Int) -> Rational
triangleFactor (tja, tjb, tjc) = triangleFactorRaw jjj (jjja, jjjb, jjjc)
  where !jjja = (tjc - tja + tjb) `quot` 2
        !jjjb = (tja - tjb + tjc) `quot` 2
        !jjjc = (tjb - tjc + tja) `quot` 2
        !jjj  = (tja + tjb + tjc) `quot` 2 + 1

{-# INLINABLE triangleFactorRaw #-}
triangleFactorRaw :: Int -> (Int, Int, Int) -> Rational
triangleFactorRaw jjj (jjja, jjjb, jjjc) =
  factorial jjju * factorial jjjv % fallingFactorialI jjj (jjj - jjjw)
  where [!jjju, !jjjv, !jjjw] = sort [jjja, jjjb, jjjc]

{-# INLINABLE getTriangularTjs #-}
getTriangularTjs :: Int -> (Int, Int) -> [Int]
getTriangularTjs tjMax (tja, tjb) = [tjmin, tjmin + 2 .. tjmax]
  where tjmin = abs (tja - tjb)
        tjmax = min tjMax (tja + tjb)

{-# INLINABLE tetrahedralSum #-}
tetrahedralSum :: (Int, Int, Int, Int, Int, Int) -> Integer
tetrahedralSum (tja, tje, tjf, tjd, tjb, tjc) =
  sum [ toInteger (minusOnePow k)
      * binomialI (k + 1) (k - (tja + tjb + tjc) `quot` 2)
      * binomialI ((tjc - tja + tjb) `quot` 2) (k - (tja + tje + tjf) `quot` 2)
      * binomialI ((tja - tjb + tjc) `quot` 2) (k - (tjd + tjb + tjf) `quot` 2)
      * binomialI ((tjb - tjc + tja) `quot` 2) (k - (tjd + tje + tjc) `quot` 2)
      | k <- [kmin .. kmax] ]

  where

    kmin =
      maximum
      [ tja + tjb + tjc
      , tjd + tje + tjc
      , tjd + tjb + tjf
      , tja + tje + tjf ] `quot` 2

    kmax =
      minimum
      [ tja + tjd + tjb + tje
      , tjb + tje + tjc + tjf
      , tja + tjd + tjc + tjf ] `quot` 2

{-# INLINABLE getBitriangularTjs #-}
getBitriangularTjs :: Int -> ((Int, Int), (Int, Int)) -> [Int]
getBitriangularTjs tjMax ((tja, tjb), (tjc, tjd)) = [tjmin, tjmin + 2 .. tjmax]
  where tjmin = max (abs (tja - tjb)) (abs (tjc - tjd))
        tjmax = minimum [tjMax, tja + tjb, tjc + tjd]

{-# INLINABLE getTms #-}
getTms :: Int -> [Int]
getTms tj = [-tj, -tj + 2 .. tj]

{-# INLINABLE get3tjms #-}
get3tjms :: Int -> [(Int, Int, Int, Int, Int, Int)]
get3tjms tjMax = do
  tj1 <- [0 .. tjMax]
  tj2 <- [0 .. tjMax]
  tj3 <- getTriangularTjs tjMax (tj1, tj2)
  tm1 <- getTms tj1
  tm2 <- getTms tj2
  let tm3 = -(tm1 + tm2)
  guard (abs tm3 <= tj3)
  pure (tj1, tm1, tj2, tm2, tj3, tm3)

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

{-# INLINABLE tuple6ToList #-}
tuple6ToList :: (a, a, a, a, a, a) -> [a]
tuple6ToList (a, b, c, d, e, f) = [a, b, c, d, e, f]

{-# INLINABLE tuple9ToList #-}
tuple9ToList :: (a, a, a, a, a, a, a, a, a) -> [a]
tuple9ToList (a, b, c, d, e, f, g, h, i) = [a, b, c, d, e, f, g, h, i]
