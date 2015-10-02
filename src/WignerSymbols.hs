-- | Stability : stable
--
-- Clebsch-Gordan coefficients and Wigner /n/-j symbols.
--
-- Note that all @j@ or @m@ arguments are represented via integers equal to
-- /twice/ their mathematical values.  To make this distinction clear, we
-- label these variables @tj@ or @tm@.
--
-- The current implementation uses the exact formulas described by
-- <https://doi.org/10.1016/S0010-4655(99)00232-5 L. Wei (1999)>
-- (<http://meghnad.iucaa.ernet.in/~tarun/pprnt/compute/ClebADKL.pdf PDF>).
--
module WignerSymbols
  (

    -- * 'SignedSqrtRational'
    SignedSqrtRational
  , ssr_new
  , ssr_split
  , ssr_signum
  , ssr_numerator
  , ssr_denominator
  , ssr_approx

    -- * Coupling/uncoupling coefficients
  , clebschGordan
  , clebschGordanSq
  , wigner3j
  , wigner3jSq

    -- * Recoupling coefficients
  , wigner6j
  , wigner6jSq
  , wigner9j
  , wigner9jSq

  ) where
import WignerSymbols.Internal
