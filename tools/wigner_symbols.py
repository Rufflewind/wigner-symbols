#!/usr/bin/env python
#
# This is a Python implementation of Clebsch-Gordan coefficients.  We use it
# to verify the results and also to compare with Sympy as well as our Haskell
# implementation.

from fractions import Fraction

try:
    long
except NameError:
    _int = int
else:
    _int = long

def factorial(n):
    r = _int(1)
    for i in range(2, int(n + 1)):
        r *= i
    return r

def clebschgordansq(tj1, tm1, tj2, tm2, tj12, tm12):
    '''Compute the sign and the square of a Clebsch-Gordan coefficient.
    Each argument must be an integer or half-integer.'''
    j1 = tj1 / 2.
    m1 = tm1 / 2.
    j2 = tj2 / 2.
    m2 = tm2 / 2.
    j12 = tj12 / 2.
    m12 = tm12 / 2.
    if (m1 + m2 != m12 or
        j1 + j2 < j12 or
        abs(j1 - j2) > j12 or
        int(2 * (j1 + j2 + j12)) % 2):
        return 0, Fraction(0, 1)
    kmin = -int(min(
        0,
        j12 - j2 + m1,
        j12 - j1 - m2,
    ))
    kmax = int(min(
        j1 + j2 - j12,
        j1 - m1,
        j2 + m2,
    ))

    # # old way of doing summation:
    # r = sum(Fraction(
    #     (-1) ** k,
    #     (
    #         factorial(k) *
    #         factorial(j1 + j2 - j12 - k) *
    #         factorial(j1 - m1 - k) *
    #         factorial(j2 + m2 - k) *
    #         factorial(j12 - j2 + m1 + k) *
    #         factorial(j12 - j1 - m2 + k)
    #     )
    # ) for k in range(kmin, kmax + 1))

    # here we calculate the summation; we accumulate the factorials to avoid
    # repeating the same calculations; this gives a 5% speed increase
    if kmin > kmax:
        return 0, Fraction(0, 1)
    else:
        c1 = int(kmin)
        c2 = int(j1 + j2 - j12 - kmin)
        c3 = int(j1 - m1 - kmin)
        c4 = int(j2 + m2 - kmin)
        c5 = int(j12 - j2 + m1 + kmin)
        c6 = int(j12 - j1 - m2 + kmin)
        c = Fraction(
            (-1) ** kmin,
            (
                factorial(c1) *
                factorial(c2) *
                factorial(c3) *
                factorial(c4) *
                factorial(c5) *
                factorial(c6)
            )
        )
        r = c
        for k in range(kmin + 1, kmax + 1):
            c1 += 1
            c5 += 1
            c6 += 1
            c *= Fraction(-c2 * c3 * c4, c1 * c5 * c6)
            c2 -= 1
            c3 -= 1
            c4 -= 1
            r += c

    return 1 if r > 0 else -1 if r < 0 else 0, (
        Fraction(
            (
                int(2 * j12 + 1) *
                factorial(j12 + j1 - j2) *
                factorial(j12 - j1 + j2) *
                factorial(j1 + j2 - j12)
            ),
            factorial(j1 + j2 + j12 + 1)
        ) *
        factorial(j12 + m12) *
        factorial(j12 - m12) *
        factorial(j1 - m1) *
        factorial(j1 + m1) *
        factorial(j2 - m2) *
        factorial(j2 + m2) *
        r ** 2
    )

def clebschgordan(tj1, tm1, tj2, tm2, tj12, tm12):
    from math import sqrt
    s, r = clebschgordansq(tj1, tm1, tj2, tm2, tj12, tm12)
    return s * sqrt(r)

def md5sum(filename, blocksize=65536):
    import hashlib
    h = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            h.update(block)
    return h.hexdigest()

def mkdirs(path):
    import os
    try:
        os.makedirs(path)
    except OSError:
        pass

def show_signedsqrtfrac(r):
    s = str((1 if r >= 0 else -1) * r ** 2)
    return s if "/" in s else s + "/1"

def iter_12tjms(tjmax):
    for tj1 in range(tjmax + 1):
        for tj2 in range(tjmax + 1):
            for tj12 in range(abs(tj1 - tj2), tj1 + tj2 + 1, 2):
                if tj12 > tjmax:
                    continue
                for tm1 in range(-tj1, tj1 + 1, 2):
                    for tm2 in range(-tj2, tj2 + 1, 2):
                        tm12 = tm1 + tm2
                        if (tj1 + tj2 + tj12) % 2 or abs(tm12) > tj12:
                            continue
                        yield tj1, tm1, tj2, tm2, tj12, tm12

def triangular_range(tja, tjb, tjmax):
    return range(abs(tja - tjb), min(tja + tjb, tjmax) + 1, 2)

def bitriangular_range(tja, tjb, tjc, tjd, tjmax):
    # (tja, tjb) and (tjc, tjd) must be of the same parity;
    # this automatically holds if the two triangles share a common edge
    assert (tja + tjb) % 2 == (tjc + tjd) % 2
    return range(
        max(abs(tja - tjb), abs(tjc - tjd)),
        min(tja + tjb, tjc + tjd, tjmax) + 1,
        2
    )

def iter_6tjs(tjmax):
    return [(tja, tjb, tjc, tjd, tje, tjf)
            for tja in range(tjmax + 1)
            for tjb in range(tjmax + 1)
            for tjc in triangular_range(tja, tjb, tjmax)
            for tjd in range(tjmax + 1)
            for tje in triangular_range(tjd, tjc, tjmax)
            for tjf in bitriangular_range(tja, tje, tjd, tjb, tjmax)]

def iter_9tjs(tjmax):
    return [(tja, tjb, tjc, tjd, tje, tjf, tjg, tjh, tji)
            for tja in range(tjmax + 1)
            for tjb in range(tjmax + 1)
            for tjc in triangular_range(tja, tjb, tjmax)
            for tjd in range(tjmax + 1)
            for tje in range(tjmax + 1)
            for tjf in triangular_range(tjd, tje, tjmax)
            for tjg in triangular_range(tja, tjd, tjmax)
            for tjh in triangular_range(tjb, tje, tjmax)
            for tji in bitriangular_range(tjg, tjh, tjc, tjf, tjmax)]

def dump_output(name, tjmax):
    mkdirs("dist")
    filename = "dist/{0}-tj{1}.txt".format(name, tjmax)
    with open(filename, "w") as f:
        yield f.write
    print(md5sum(filename) + "  " + filename)

def dump_cg(tjmax, use_sympy=False):
    name = "cg"
    if use_sympy:
        import sympy.physics.wigner as spw
        name = "sym" + name
    for write in dump_output(name, tjmax):
        for tj1, tm1, tj2, tm2, tj12, tm12 in iter_12tjms(tjmax):
            if use_sympy:
                z = spw.clebsch_gordan(
                    Fraction(tj1, 2), Fraction(tj2, 2), Fraction(tj12, 2),
                    Fraction(tm1, 2), Fraction(tm2, 2), Fraction(tm12, 2),
                )
                write("\t".join(map(str, (
                    tj1, tm1, tj2, tm2, tj12, tm12, show_signedsqrtfrac(z),
                ))) + "\n")
            else:
                s, r = clebschgordansq(tj1, tm1, tj2, tm2, tj12, tm12)
                write("\t".join(map(str, (
                    tj1, tm1, tj2, tm2, tj12, tm12,
                    "{0}/{1}".format(s * r.numerator, r.denominator)
                ))) + "\n")

def dump_w6j(tjmax):
    import sympy.physics.wigner as spw
    for write in dump_output("symw6j", tjmax):
        for tjs in iter_6tjs(tjmax):
            z = spw.wigner_6j(*(Fraction(tj, 2) for tj in tjs))
            write("\t".join(map(str, tjs + (show_signedsqrtfrac(z),))) + "\n")

def dump_w9j(tjmax):
    # note: sympy's wigner_9j is buggy for half-integer arguments
    #       (as of 0.7.6.1)
    import sympy.physics.wigner as spw
    for write in dump_output("symw9j", tjmax):
        for tjs in iter_9tjs(tjmax):
            z = spw.wigner_9j(*(Fraction(tj, 2) for tj in tjs))
            write("\t".join(map(str, tjs + (show_signedsqrtfrac(z),))) + "\n")

def main():
    # at tjmax = 25,
    #
    # - Clesbch-Gordan coefficients take ~2min in our
    #   implementation, while sympy's will be ~10x slower
    # - sympy's wigner_6j is even slower (~8 hours)
    # - sympy's wigner_9j is extremely slow; tjmax = 10 takes 26 hours;
    #   don't think tjmax = 25 is possible in a reasonable timeframe
    #
    tjmax = 25
    use_sympy = False

    dump_cg(tjmax, use_sympy)
    dump_w6j(tjmax)
    #dump_w9j(tjmax)

if __name__ == "__main__":
    main()
