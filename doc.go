// Package predicates provides robust geometric predicates for float64
// coordinates: orientation and in-circle / in-sphere tests whose sign is
// exact, at close to plain floating-point cost.
//
// It is an independent Go implementation of the method of
//
//	J. R. Shewchuk, Adaptive Precision Floating-Point Arithmetic and Fast
//	Robust Geometric Predicates, Discrete & Computational Geometry
//	18(3):305-363, 1997 (CMU-CS-96-140R),
//
// written from that paper and from Shewchuk's public-domain predicates.c
// (https://www.cs.cmu.edu/~quake/robust.html): floating-point expansions,
// error-free sums and products, and adaptive evaluation in stages A (a
// floating-point filter with a proven error bound), B, C and D (exact).
//
// # Exactness
//
// Every function except the ...Fast ones returns a value whose sign is the
// sign of the exact determinant of its float64 inputs. Products are wrapped in
// explicit float64 conversions so that the Go compiler never fuses them into
// fused multiply-adds: results are identical on every architecture.
//
// # Range
//
// Floating-point expansions are exact only while nothing overflows or
// underflows. Each predicate therefore checks that every nonzero coordinate
// lies in a range that depends on the degree of its determinant:
//
//	Orient2D  [2^-400, 2^480]
//	Orient3D  [2^-250, 2^320]
//	InCircle  [2^-170, 2^240]
//	InSphere  [2^-130, 2^190]
//
// Inputs outside the range are still answered exactly, by a slower
// big-integer evaluation. NaN and ±Inf yield NaN.
package predicates
