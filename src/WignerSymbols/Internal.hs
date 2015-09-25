-- TODO: Remove NMR later
{-# LANGUAGE NoMonomorphismRestriction #-}
module WignerSymbols.Internal where
import Data.Foldable (foldl')
import Data.Ratio (Ratio, (%), numerator)

newtype SignedSqrtRatio a = SignedSqrtRatio (Ratio a)
                          deriving (Eq, Read, Show)

type SignedSqrtRational = SignedSqrtRatio Integer

splitSignedSqrtRatio :: (Integral a, Num b) => SignedSqrtRatio a -> (b, Ratio a)
splitSignedSqrtRatio (SignedSqrtRatio x) =
  (fromIntegral (signum (numerator x)), abs x)

approxSignedSqrtRatio :: (Integral a, Floating b) => SignedSqrtRatio a -> b
approxSignedSqrtRatio x = case splitSignedSqrtRatio x of
  (s, r) -> s * sqrt (realToFrac r)

factorial :: Int -> Integer
factorial = go 1 . toInteger
  -- difference between toInteger first or toInteger per iteration seems negligible
  where go r n
          | n <= 1    = r
          | otherwise = go (r * n) (pred n)

-- Difference between _simple and the other one seems negligible
factorial_simple :: Int -> Integer
factorial_simple n = product [1 .. toInteger n]

binomial :: (Fractional a, Integral b) => a -> b -> a
binomial = binomial_descending (\ x y -> x / fromIntegral y)

binomialI :: (Integral a, Integral b, Num c) => a -> b -> c
binomialI n k = fromInteger (binomial_ascendingI (fromIntegral n) (fromIntegral k))
  -- don't use the one below as it's ~twice as slow:
  -- \ n k -> numerator (binomial_descending (%) n k)

binomial_descending :: (Num a, Integral b, Num c) => (a -> b -> c) -> a -> b -> c
binomial_descending divide n k
  | k > 0     = divide n k * binomial_descending divide (n - 1) (k - 1)
  | otherwise = 1

binomial_ascending :: Integral a => a -> a -> a
binomial_ascending = go 1 1
  where go r i n k
          | i > k     = r
          | otherwise = go (r * n `quot` i) (succ i) (pred n) k

binomial_ascendingI :: Int -> Int -> Integer
binomial_ascendingI = go 1 1
  where go r i n k
          | i > k     = r
          | otherwise = go (r * toInteger n `quot` toInteger i) (succ i) (pred n) k

minusOnePow :: Integral a => a -> a
minusOnePow n = 1 - n `mod` 2 * 2

-- note: -fllvm makes this faster (~10%)

triangleCondition :: (Num a, Ord a) => a -> a -> a -> Bool
triangleCondition a b c = abs (a - b) <= c && c <= a + b

-- | Compute the sign and the square of a Clebsch-Gordan coefficient.
--   Input arguments are /twice/ the usual quantum numbers.
clebschGordanSq :: Int -> Int -> Int -> Int -> Int -> Int -> SignedSqrtRational
clebschGordanSq = clebschGordanSq_b

clebschGordanSq_b :: Int -> Int -> Int -> Int -> Int -> Int -> SignedSqrtRational
clebschGordanSq_b tj1 tm1 tj2 tm2 tj12 tm12 =
  SignedSqrtRatio (z * fromIntegral (tj12 + 1))
  where SignedSqrtRatio z = wigner3jSqRaw tj1 tm1 tj2 tm2 tj12 (-tm12)

wigner3jSq :: Int -> Int -> Int -> Int -> Int -> Int -> SignedSqrtRational
wigner3jSq tj1 tm1 tj2 tm2 tj3 tm3 = SignedSqrtRatio (s * z)
  where s = fromIntegral (minusOnePow ((tj1 - tj2 - tm3) `quot` 2))
        SignedSqrtRatio z = wigner3jSqRaw tj1 tm1 tj2 tm2 tj3 tm3

{-# SPECIALIZE wigner3j :: Int -> Int -> Int -> Int -> Int -> Int -> Double #-}
{-# SPECIALIZE wigner3j :: Int -> Int -> Int -> Int -> Int -> Int -> Float #-}
wigner3j :: Floating a => Int -> Int -> Int -> Int -> Int -> Int -> a
wigner3j tj1 tm1 tj2 tm2 tj3 tm3 =
  approxSignedSqrtRatio (wigner3jSq tj1 tm1 tj2 tm2 tj3 tm3)

wigner3jSqRaw :: Int -> Int -> Int -> Int -> Int -> Int -> SignedSqrtRational
wigner3jSqRaw tj1 tm1 tj2 tm2 tj3 tm3
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
    z1 = factorial jjj1 * factorial jjj2 * factorial jjj3 % factorial (succ jjj)

    z2 :: Rational
    z2 =
      -- splitting z2 will make it slower
      (binomialI tj1 jjj1 * binomialI tj2 jjj2 * binomialI tj3 jjj3) %
      (binomialI tj1 jm1 * binomialI tj2 jm2 * binomialI tj3 jm3)

    z3 :: Integer
    z3 | kmin > kmax = 0
       | otherwise   = r
      where
        (r, _) = foldl' f (c0, c0) [succ kmin .. kmax]
        c0 = toInteger (minusOnePow kmin) *
             binomialI jjj2 kmin *
             binomialI jjj1 (jsm1 - kmin) *
             binomialI jjj3 (jm2 - kmin)
        f (s, c) k =
          let c' = c
                 * toInteger (jjj2 - k + 1) `quot` toInteger k
                 * toInteger (jsm1 - k + 1) `quot` toInteger (jjj1 - (jsm1 - k))
                 * toInteger (jm2 - k + 1) `quot` toInteger (jjj3 - (jm2 - k))
          in (s - c', -c')

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

clebschGordanSq_a :: Int -> Int -> Int -> Int -> Int -> Int -> SignedSqrtRational
clebschGordanSq_a tj1 tm1 tj2 tm2 tj12 tm12
  | conservationLaw && triangleCondition tj1 tj2 tj12 && jParity = SignedSqrtRatio (sign * surd)
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

-- | Compute the sign and the square of a Clebsch-Gordan coefficient.
--   Input arguments are //twice// the usual quantum numbers.
{-# SPECIALIZE clebschGordan :: Int -> Int -> Int -> Int -> Int -> Int -> Double #-}
{-# SPECIALIZE clebschGordan :: Int -> Int -> Int -> Int -> Int -> Int -> Float #-}
clebschGordan :: Floating a => Int -> Int -> Int -> Int -> Int -> Int -> a
clebschGordan tj1 tm1 tj2 tm2 tj12 tm12 =
  approxSignedSqrtRatio (clebschGordanSq tj1 tm1 tj2 tm2 tj12 tm12)

getTj12s :: (Enum a, Num a) => a -> a -> [a]
getTj12s tj1 tj2 = [tj12min, tj12min + 2 .. tj12max]
  where tj12min = abs (tj1 - tj2)
        tj12max = tj1 + tj2

getTms :: (Enum a, Num a) => a -> [a]
getTms tj = [-tj, -tj + 2 .. tj]
