`wigner-symbols`
================

This library calculates:

  - Clebsch-Gordan coefficients
  - Wigner 3-j symbols
  - Wigner 6-j symbols
  - Wigner 9-j symbols

These sets of numbers are commonly used in physics to couple, uncouple, and
recouple eigenstates of angular momentum and similar operators.
Mathematically, they describe the relationship between the bases of
irreducible representations of SU(2) or SO(3) and those of their tensor
products.

The library uses the Condon-Shortley phase convention as typical in physics.

Each function has a `…Sq` variant that returns the exact value as a
`SignedSqrtRational`, which represents a mathematical expression of the form:

    s × √(n / d)

where

  - `s` is a sign (`+`, `-`, or `0`),
  - `n` is a nonnegative numerator, and
  - `d` is a positive denominator.
