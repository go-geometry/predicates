// Package bench measures github.com/go-geometry/predicates against
// Shewchuk's public-domain predicates.c (copied verbatim beside this file,
// from https://www.cs.cmu.edu/afs/cs/project/quake/public/code/predicates.c),
// compiled without FMA contraction (pragmas in shewchuk.c), which the C
// code's error-free transformations forbid (the paper's section 5).
//
// It is a separate module so that cgo never enters the library's build.
//
//	cd bench && go test -run '^$' -bench . -count 6
package bench

/*
#cgo CFLAGS: -O2 -w
double loop_orient2d(double *pts, int n, long iters);
double loop_incircle(double *pts, int n, long iters);
double loop_orient3d(double *pts, int n, long iters);
double loop_insphere(double *pts, int n, long iters);
void exactinit(void);
*/
import "C"

import "unsafe"

func init() { C.exactinit() }

// Kind selects a C loop.
type Kind int

const (
	Orient2D Kind = iota
	InCircle
	Orient3D
	InSphere
)

// Loop runs the C predicate iters times over pts (records of 8 or 15
// float64) and returns the sum of the results, so nothing is optimised away.
func Loop(k Kind, pts []float64, records int, iters int) float64 {
	p := (*C.double)(unsafe.Pointer(&pts[0]))
	n, it := C.int(records), C.long(iters)
	switch k {
	case Orient2D:
		return float64(C.loop_orient2d(p, n, it))
	case InCircle:
		return float64(C.loop_incircle(p, n, it))
	case Orient3D:
		return float64(C.loop_orient3d(p, n, it))
	}
	return float64(C.loop_insphere(p, n, it))
}
