# predicates

Robust geometric predicates for Go: orientation and in-circle tests on
`float64` coordinates whose **sign is exact**, at close to plain
floating-point cost.

An independent Go implementation of the adaptive-precision method of
J. R. Shewchuk, *Adaptive Precision Floating-Point Arithmetic and Fast Robust
Geometric Predicates*, Discrete & Computational Geometry 18(3):305-363, 1997,
written from the paper and from his public-domain
[`predicates.c`](https://www.cs.cmu.edu/~quake/robust.html).

```go
import "github.com/go-geometry/predicates"

predicates.Orient2D(a, b, c)    // > 0 counterclockwise, < 0 clockwise, 0 collinear
predicates.InCircle(a, b, c, d) // > 0 if d is inside the circle through a, b, c (ccw)
predicates.Orient3D(a, b, c, d) // > 0 if d is below the plane of a, b, c (ccw seen from above)
predicates.InSphere(a, b, c, d, e) // > 0 if e is inside the sphere through a..d (Orient3D > 0)
```

| Function | |
|---|---|
| `Orient2D`, `InCircle`, `Orient3D`, `InSphere` | adaptive: a floating-point filter (stage A), then stages B, C, D only when the filter cannot certify the sign |
| `...Exact` | exact evaluation without the adaptive stages |
| `...Fast` | plain floating-point evaluation, **not robust**, for comparison |

- **Exact sign** for every finite input. Coordinates outside the range in
  which floating-point expansions cannot overflow or underflow (see the
  package documentation) are answered by an exact big-integer evaluation.
  NaN or ±Inf gives NaN.
- **Identical results on every architecture**: every product is rounded by an
  explicit `float64(...)` conversion, so Go never fuses it into an FMA.

## Validation

The sign of every predicate is compared with an exact `math/big` rational
evaluation on random, near-degenerate (rounded onto a line or circle, then
nudged by a few ulps), exactly degenerate (collinear lattices, Pythagorean
circles, repeated points under inexact translations), Shewchuk's tilted
grid, and huge / tiny / mixed-magnitude inputs: zero mismatches. Planted
breaks (stage A or stage B error bound removed) must produce wrong signs on
the same corpus, or the test fails. 100 % statement coverage.

Part of [go-geometry](https://go-geometry.github.io/); design:
[go-geometry design, section 4](https://go-geometry.github.io/docs/main/design/).

## Licence

BSD-3-Clause. `predicates.c`, which this package was written from, is in the
public domain.
