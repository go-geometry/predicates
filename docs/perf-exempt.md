# Public functions without a performance case

Read by `refbench` in the CI mode (`cd bench && go test -run TestPerf .`): every
public function of the module has a case in `bench/perf_test.go` or is named
here with its reason. The method and the table are in [perf.md](perf.md).

- `Orient2DFast`: the filter-only form; predicates.c's `orient2dfast` exists but `bench/loops.c` has no loop for it yet.
- `InCircleFast`: as `Orient2DFast`, `incirclefast`.
- `Orient3DFast`: as `Orient2DFast`, `orient3dfast`.
- `InSphereFast`: as `Orient2DFast`, `inspherefast`.
- `Orient2DExact`: the exact stage alone; predicates.c has `orient2dexact` but no loop calls it. The adaptive function's degenerate cases (`Orient2D/degenerate`) time the exact stage in context.
- `InCircleExact`: as `Orient2DExact`, `incircleexact`.
- `Orient3DExact`: as `Orient2DExact`, `orient3dexact`.
- `InSphereExact`: as `Orient2DExact`, `insphereexact`.
- `OrientN`: no reference; predicates.c has the 2-D and 3-D predicates only, and the exact n-D determinant is checked against math/big in `fuzz_test.go`.
- `InSphereN`: as `OrientN`.
