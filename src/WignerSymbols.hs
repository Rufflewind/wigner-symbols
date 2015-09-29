-- | Clebsch-Gordan coefficients and Wigner /n/-j symbols.
--
-- Note that all @j@ or @m@ arguments are represented via integers equal to
-- /twice/ their mathematical values.  To make this distinction clear, we name
-- these variables @tj@ or @tm@.
--
module WignerSymbols
  (
    -- * 'SignedSqrtRatio'
    SignedSqrtRatio
  , SignedSqrtRational
  , ssr_signum
  , ssr_numerator
  , ssr_denominator
  , ssr_approx
  , ssr_split

    -- * Coupling/uncoupling coefficients
  , clebschGordan
  , clebschGordanSq
  , wigner3j
  , wigner3jSq

    -- * Recoupling coefficients
  , wigner6j
  , wigner6jSq
  ) where
import WignerSymbols.Internal
