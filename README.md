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
| `OrientN`, `InSphereN` | exact orientation and in-sphere tests in any dimension 1 ≤ d ≤ 16: floating-point expansion determinants up to order 8, big-integer elimination above |
| `...Exact` | exact evaluation without the adaptive stages |
| `...Fast` | plain floating-point evaluation, **not robust**, for comparison |

- **Exact sign** for every finite input. Coordinates outside the range in
  which floating-point expansions cannot overflow or underflow (see the
  package documentation) are answered by an exact big-integer evaluation.
  NaN or ±Inf gives NaN.
- **Identical results on every architecture**: every product is rounded by an
  explicit `float64(...)` conversion, so Go never fuses it into an FMA.
- **Bounded cost** for untrusted input: `OrientN`/`InSphereN` accept
  d ≤ 16 (`MaxDim`); NaN/±Inf are refused with `ErrNotFinite`.

## Validation

- **Sign against exact arithmetic.** Every predicate's sign is compared with
  a `math/big` rational evaluation on random, near-degenerate (rounded onto a
  line, circle, plane or sphere, then nudged by a few ulps), exactly
  degenerate (collinear and coplanar lattices, Pythagorean circles and
  spheres, repeated points under inexact translations), Shewchuk's tilted
  grid, and huge / tiny / mixed-magnitude inputs: zero mismatches. The
  non-robust `Fast` evaluation is wrong on hundreds of these cases.
- **Planted breaks.** With the stage A or stage B error bound removed, the
  same corpus must produce wrong signs, or the test fails.
- **Fuzzing.** `FuzzOrient2D`, `FuzzInCircle`, `FuzzOrient3D`,
  `FuzzInSphere`, `FuzzNDim` check the sign against `math/big` on arbitrary
  bits (NaN, ±Inf, subnormals included); the corpus is committed under
  `testdata/fuzz`.
- **Same bits everywhere.** `TestBitsIdenticalAcrossArchitectures` hashes the
  bits of 5 600 results; CI checks the hash on amd64, arm64, riscv64,
  loong64, ppc64le, s390x, 386 and arm. Removing one `float64(...)` from
  `Orient2DFast` changes the hash on arm64, where Go then fuses an FMA.
- **Against `predicates.c`.** `bench/` (a separate module, cgo) compiles
  Shewchuk's code without FMA contraction and checks equal signs on every
  benchmark input.
- 100 % statement coverage, `-race`, govulncheck.

## Performance

```sh
go test -run '^$' -bench . -benchmem                  # fast vs adaptive vs exact
cd bench && go test -run TestPerf .                   # CI mode: allocations, table, public functions
cd bench && go test -run TestPerf -cpu 1 -perf.run .  # against predicates.c, load-gated
```

The comparison with `predicates.c` runs through
[refbench](https://github.com/go-scikit/refbench): interleaved rounds (C, Go,
then Go, C), one core, the 1-minute load recorded for every round and a refusal
to measure at 8 or more; method and status in [docs/perf.md](docs/perf.md), the
generated table in [docs/perf-table.md](docs/perf-table.md). Timing with the new
tool is deferred (the host's load was above the gate); the table below is from
the previous harness, which ran all the C rounds and then all the Go rounds.

Against `predicates.c` (clang `-O2`, no FMA contraction), same records,
direct calls on both sides, Apple M4 Max, Go 1.27.3, single core, median of
6 rounds. *Random*: points in [-1, 1), decided by the filter. *Degenerate*:
one point repeated under inexact translations, exact zero, run to stage D.

| | `predicates.c` | Go | Go allocs |
|---|---:|---:|---:|
| Orient2D random | 1.79 ns | 4.37 ns | 0 |
| InCircle random | 3.97 ns | 9.08 ns | 0 |
| Orient3D random | 3.57 ns | 11.7 ns | 0 |
| InSphere random | 10.5 ns | 17.2 ns | 0 |
| Orient2D degenerate | 1.29 ns | 4.19 ns | 0 |
| InCircle degenerate | 3.52 µs | 5.57 µs | 18 |
| Orient3D degenerate | 0.32 µs | 2.2 µs | 7 |
| InSphere degenerate | 5.5 µs | 18 µs | 35 |

⚠ The 1-minute load was 9.7–10.7 during this run, above the threshold the
fleet measures at: times are indicative (the degenerate rows vary by up to
±40 %), allocation counts are exact. **Go is slower than `predicates.c`
on every row**, an open defect,
[#4](https://github.com/go-geometry/predicates/issues/4): the fast path pays a
range check and by-value array arguments; stage D evaluates a full
expansion determinant instead of Shewchuk's incremental tail products.

Part of [go-geometry](https://go-geometry.github.io/); design:
[go-geometry design, section 4](https://go-geometry.github.io/docs/latest/design/);
documentation: [predicates](https://go-geometry.github.io/docs/latest/predicates/).

## Licence

BSD-3-Clause. `predicates.c`, which this package was written from, is in the
public domain.
