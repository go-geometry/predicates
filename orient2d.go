package predicates

// Orient2D reports the orientation of the triangle a, b, c: a positive value
// if the points occur in counterclockwise order, a negative value if they
// occur in clockwise order, and zero if they are collinear. The sign is exact
// for all finite inputs; the magnitude approximates twice the signed area,
// det[[ax-cx, ay-cy], [bx-cx, by-cy]].
//
// The evaluation is adaptive (Shewchuk 1997, section 4.3): a floating-point
// filter (stage A) decides almost every call; stages B, C and D add precision
// only when the filter cannot certify the sign. A coordinate outside the
// range where the error analysis holds (see the package documentation) sends
// the call to an exact big-integer evaluation instead. NaN or ±Inf in any
// coordinate yields NaN.
func Orient2D(a, b, c [2]float64) float64 {
	if !inRange2(a, b, c, lo2, hi2) {
		return orient2dBig(a, b, c)
	}
	detleft := float64((a[0] - c[0]) * (b[1] - c[1]))
	detright := float64((a[1] - c[1]) * (b[0] - c[0]))
	det := detleft - detright
	var detsum float64
	switch {
	case detleft > 0:
		if detright <= 0 {
			return det
		}
		detsum = detleft + detright
	case detleft < 0:
		if detright >= 0 {
			return det
		}
		detsum = -detleft - detright
	default:
		return det
	}
	errbound := ccwErrBoundA * detsum
	if det >= errbound || -det >= errbound {
		return det
	}
	return orient2dAdapt(a, b, c, detsum)
}

// Orient2DFast is the plain floating-point evaluation of Orient2D's
// determinant. It is not robust: its sign may be wrong for nearly collinear
// points. It exists for comparison and benchmarks.
func Orient2DFast(a, b, c [2]float64) float64 {
	return float64((a[0]-c[0])*(b[1]-c[1])) - float64((a[1]-c[1])*(b[0]-c[0]))
}

// Orient2DExact evaluates Orient2D's determinant exactly, without the
// adaptive stages, and returns a value with its sign.
func Orient2DExact(a, b, c [2]float64) float64 {
	if !inRange2(a, b, c, lo2, hi2) {
		return orient2dBig(a, b, c)
	}
	return orient2dExact(a, b, c)
}

// orient2dAdapt runs stages B, C and D (paper, Figure 21).
func orient2dAdapt(a, b, c [2]float64, detsum float64) float64 {
	acx := a[0] - c[0]
	bcx := b[0] - c[0]
	acy := a[1] - c[1]
	bcy := b[1] - c[1]

	// Stage B: the exact determinant of the rounded differences.
	l1, l0 := twoProduct(acx, bcy)
	r1, r0 := twoProduct(acy, bcx)
	bexp := twoTwoDiff(l1, l0, r1, r0)
	det := estimate(bexp[:])
	errbound := ccwErrBoundB * detsum
	if det >= errbound || -det >= errbound {
		return det
	}

	acxt := twoDiffTail(a[0], c[0], acx)
	bcxt := twoDiffTail(b[0], c[0], bcx)
	acyt := twoDiffTail(a[1], c[1], acy)
	bcyt := twoDiffTail(b[1], c[1], bcy)
	if acxt == 0 && acyt == 0 && bcxt == 0 && bcyt == 0 {
		// The differences were exact: B is the determinant.
		return det
	}

	// Stage C: first-order correction by the tails.
	errbound = ccwErrBoundC*detsum + resultErrBound*abs(det)
	det += (float64(acx*bcyt) + float64(bcy*acxt)) - (float64(acy*bcxt) + float64(bcx*acyt))
	if det >= errbound || -det >= errbound {
		return det
	}

	// Stage D: B plus the three exact tail products.
	var c1, c2, d [16]float64
	u := cross2(acxt, bcy, acyt, bcx)
	s := sumExp(c1[:0], bexp[:], u[:])
	u = cross2(acx, bcyt, acy, bcxt)
	s = sumExp(c2[:0], s, u[:])
	u = cross2(acxt, bcyt, acyt, bcxt)
	s = sumExp(d[:0], s, u[:])
	return s[len(s)-1]
}

// orient2dExact is stage D computed directly from the exact differences.
func orient2dExact(a, b, c [2]float64) float64 {
	var back [64]float64
	ar := &arena{buf: back[:0]}
	acx := diffExp(ar, a[0], c[0])
	bcx := diffExp(ar, b[0], c[0])
	acy := diffExp(ar, a[1], c[1])
	bcy := diffExp(ar, b[1], c[1])
	det := addExp(ar, mulExp(ar, acx, bcy), mulExp(ar, negExp(ar, acy), bcx))
	return det[len(det)-1]
}

// orient2dBig is the out-of-range path.
func orient2dBig(a, b, c [2]float64) float64 {
	if !finite(a[0], a[1], b[0], b[1], c[0], c[1]) {
		return nan
	}
	s := scaleOf(a[0], a[1], b[0], b[1], c[0], c[1])
	m := [][]*bigInt{
		{bigDiff(a[0], c[0], s), bigDiff(a[1], c[1], s)},
		{bigDiff(b[0], c[0], s), bigDiff(b[1], c[1], s)},
	}
	return bigDet(m, 2*s)
}
