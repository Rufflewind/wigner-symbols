{-# LANGUAGE BangPatterns, CPP, ScopedTypeVariables #-}
module Main (main) where
#if !MIN_VERSION_base(4, 8, 0)
import Control.Applicative (pure)
import Data.Functor ((<$>))
#endif
import Control.Monad.Primitive (PrimMonad, PrimState)
import Control.Monad.ST (runST)
import Data.Foldable (for_)
import Data.Monoid ((<>))
import Data.Vector (Vector)
import Data.Vector.Generic (Mutable)
import Data.Vector.Unboxed (Unbox)
import qualified Data.Vector.Generic as Vector
import qualified Data.Vector.Generic.Mutable as MVector
import qualified Data.Vector.Unboxed as Vector_Unboxed
import System.Random (randomRIO)
import Criterion.Main
import WignerSymbols
import WignerSymbols.Internal

type Vector_Unboxed a = Vector_Unboxed.Vector a

mvector_traverse :: (PrimMonad m, MVector.MVector v a) =>
                    Int -> (a -> m a) -> v (PrimState m) a -> m ()
mvector_traverse i f v = MVector.read v i >>= f >>= MVector.write v i

------------------------------------------------------------------------------

main :: IO ()
main =
  defaultMain
  [ envTable 25 getTable3tjms $ \ t ->
    bgroup ""
    [ bgroup "cg" (benchPart 5 t flippedClebschGordan)
    , bgroup "w3j" (benchPart 5 t wigner3j)
    ]
  , envTable 20 getTable6tjs $ \ t ->
    bgroup "w6j" (benchPart 5 t wigner6jSq)
  , bench "cg.tj=8" (whnf clebschGordanSq (8, 0, 8, 0, 8, 0))
  , bench "cg.tj=100" (whnf clebschGordanSq (100, 0, 100, 0, 100, 0))
  , bench "w6j.tj=8" (whnf wigner6jSq (10, 10, 10, 10, 10, 10))
  ]

------------------------------------------------------------------------------

type Partition a = (Vector_Unboxed Int, Vector_Unboxed a)

categorize :: forall a . Unbox a => Int -> [a] -> (a -> Int) -> Partition a
categorize numCategories xs classify = (offsets, table)
  where

    defaultCap = 25

    tableS :: Vector (Vector_Unboxed a)
    tableS = runST $ do
      mtableG <- MVector.replicateM numCategories (gmv_new defaultCap)
      for_ xs $ \ x ->
        mvector_traverse (classify x) (`gmv_append` x) mtableG
      tableG <- Vector.freeze mtableG
      Vector.mapM gmv_unsafeFreeze tableG

    table :: Vector_Unboxed a
    table  = Vector.concat (Vector.toList tableS)

    offsetsB :: Vector Int
    offsetsB = Vector.postscanl' (+) 0 (Vector.length <$> tableS)

    offsets :: Vector_Unboxed Int
    offsets  = Vector.convert offsetsB

partRange :: Partition a -> (Int, Int) -> (Int, Int)
partRange (offsets, _) (partMin, partMax) = (lowerIndex, upperIndex)
  where
    lowerIndex | partMin == 0 = 0
               | otherwise    = offsets Vector.! (partMin - 1)
    upperIndex = offsets Vector.! partMax - 1

randomPartElem :: Unbox a => Partition a -> (Int, Int) -> IO a
randomPartElem (_, indices) range =
  (indices Vector.!) <$> randomRIO range

benchPart :: Unbox a => Int -> (Int, Partition a) -> (a -> b) -> [Benchmark]
benchPart tjStep (tableTjMax, inputTable) f =
  [ let tjMin = tjMax - tjStep
        range = partRange inputTable (tjMin, tjMax) in
    bench ("tj=" <> show tjMin <> ".." <> show tjMax)
          (whnfIO (f <$> randomPartElem inputTable range))
  | tjMax <- [tjStep, tjStep * 2 .. tableTjMax] ]

------------------------------------------------------------------------------

type Growable v =
  ( Int -- actual length
  , v )

gmv_new :: (PrimMonad m, MVector.MVector v a) =>
           Int -> m (Growable (v (PrimState m) a))
gmv_new cap = do
  v <- MVector.new cap
  pure (0, v)

gmv_append :: (PrimMonad m, MVector.MVector v a) =>
              Growable (v (PrimState m) a)
           -> a
           -> m (Growable (v (PrimState m) a))
gmv_append (len, v) x = do
  v' <-
    let !cap = MVector.length v in
    if len < cap
      then pure v
      else MVector.grow v (cap * 2)
  MVector.write v' len x
  pure (len + 1, v')

gmv_unsafeFreeze :: (PrimMonad m, Vector.Vector v a) =>
                    Growable (Mutable v (PrimState m) a) -> m (v a)
gmv_unsafeFreeze (len, v) = Vector.unsafeFreeze (MVector.slice 0 len v)

------------------------------------------------------------------------------

type SixTk = (Int, Int, Int, Int, Int, Int)

flippedClebschGordan :: SixTk -> Double
flippedClebschGordan (tj1, tm1, tj2, tm2, tj3, tm3) =
  clebschGordan (tj1, tm1, tj2, tm2, tj3, -tm3)

getTable3tjms :: Int -> Partition SixTk
getTable3tjms tableTjMax =
  categorize (tableTjMax + 1) (get3tjms tableTjMax) $
  \ (tj1, _, tj2, _, tj3, _) ->
    maximum [tj1, tj2, tj3]

getTable6tjs :: Int -> Partition SixTk
getTable6tjs tableTjMax =
  categorize (tableTjMax + 1) (get6tjs tableTjMax) $
  \ (tj1, tj2, tj3, tj4, tj5, tj6) ->
    maximum [tj1, tj2, tj3, tj4, tj5, tj6]

envTable :: Int
         -> (Int -> Partition SixTk)
         -> ((Int, Partition SixTk) -> Benchmark)
         -> Benchmark
envTable tableTjMax getTable benchmark =
  env (pure (getTable tableTjMax)) $
  \ e -> benchmark (tableTjMax, e)
