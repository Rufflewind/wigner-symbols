{-# LANGUAGE CPP #-}
module Common
  ( Monoid
  , (<>)
  , (<$>)
  , pure
  ) where

#if !MIN_VERSION_base(4, 8, 0)
import Control.Applicative (pure)
import Data.Functor ((<$>))
#endif

#if MIN_VERSION_base(4, 8, 0)
import Data.Monoid ((<>))
#elif MIN_VERSION_base(4, 5, 0)
import Data.Monoid (Monoid, (<>))
#else
import Data.Monoid (Monoid, mappend)
#endif

#if !MIN_VERSION_base(4, 5, 0)
(<>) :: Monoid m => m -> m -> m
(<>) = mappend
#endif
