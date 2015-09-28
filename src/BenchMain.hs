{-# LANGUAGE BangPatterns, CPP, ScopedTypeVariables #-}
module Main (main) where
#if !MIN_VERSION_base(4, 8, 0)
import Control.Applicative (pure)
import Data.Functor ((<$>))
#endif
import Control.Monad.ST -- (runST)
import Data.Foldable (for_)
import Data.Monoid ((<>))
import Data.Vector (Vector)
import Data.Vector.Unboxed (Unbox)
import qualified Data.Vector.Generic as Vector
import qualified Data.Vector.Generic.Mutable as MVector
import qualified Data.Vector.Unboxed as Vector_Unboxed
import System.Random (randomRIO)
import Criterion.Main
import WignerSymbols
import WignerSymbols.Internal
type Vector_Unboxed a = Vector_Unboxed.Vector a

type ThreeTjm = (Int, Int, Int, Int, Int, Int)

type Partition a = (Vector_Unboxed Int, Vector_Unboxed a)

main :: IO ()
main =
  defaultMain
    [ env (pure (getTable3tjms tableTjMax)) $ \ ~table3tjms ->
      bgroup ""
      [ bgroup "cg" (bench3tjms table3tjms flippedClebschGordan)
      , bgroup "3j" (bench3tjms table3tjms wigner3j)
      ]
    ]
  where

    tableTjMax = 25

    tjStep = 5

    bench3tjms table3tjms f =
      [ let tjMin = tjMax - tjStep
            iRange = indexRange3tjm table3tjms tjMin tjMax in
        bench ("tj=" <> show tjMin <> ".." <> show tjMax)
              (whnfIO (f <$> random3tjm table3tjms iRange))
      | tjMax <- [tjStep, tjStep * 2 .. tableTjMax] ]

flippedClebschGordan :: ThreeTjm -> Double
flippedClebschGordan (tj1, tm1, tj2, tm2, tj3, tm3) =
  clebschGordan (tj1, tm1, tj2, tm2, tj3, -tm3)

categorize :: forall a . Unbox a => Int -> [a] -> (a -> Int) -> Partition a
categorize numCategories xs classify = (offsets, table)
  where

    tableL :: Vector [a]
    tableL = runST $ do
      mtableL <- MVector.replicate numCategories []
      for_ xs $ \ x ->
        MVector.modify mtableL (x :) (classify x)
      Vector.freeze mtableL

    tableG :: Vector (Vector_Unboxed a)
    tableG = Vector.fromList <$> tableL

    table :: Vector_Unboxed a
    table  = Vector.concat (Vector.toList tableG)

    offsetsB :: Vector Int
    offsetsB = Vector.postscanl' (+) 0 (Vector.length <$> tableG)

    offsets :: Vector_Unboxed Int
    offsets  = Vector.convert offsetsB

indexRange3tjm :: Partition ThreeTjm -> Int -> Int -> (Int, Int)
indexRange3tjm (offsets, _) tjMin tjMax = (lowerIndex, upperIndex)
  where
    lowerIndex | tjMin == 0 = 0
               | otherwise  = offsets Vector.! (tjMin - 1)
    upperIndex = offsets Vector.! tjMax - 1

random3tjm :: Partition ThreeTjm -> (Int, Int) -> IO ThreeTjm
random3tjm (_, indices) indexRange =
  (indices Vector.!) <$> randomRIO indexRange

getTable3tjms :: Int -> Partition ThreeTjm
getTable3tjms tableTjMax =
  categorize (tableTjMax + 1) (get3tjms tableTjMax) $
  \ (tj1, _, tj2, _, tj3, _) -> maximum [tj1, tj2, tj3]
