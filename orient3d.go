package predicates

import "math"

// Orient3D reports the orientation of the tetrahedron a, b, c, d: a positive
// value if d lies below the plane through a, b and c ("below" being the side
// from which a, b, c appear clockwise; equivalently a, b, c appear
// counterclockwise from above), a negative value if d lies above, zero if
// the four points are coplanar. The sign is exact for all finite inputs; the
// magnitude approximates six times the signed volume,
// det[[a-d], [b-d], [c-d]].
//
// Evaluation is adaptive, as for Orient2D. NaN or ±Inf yields NaN.
func Orient3D(a, b, c, d [3]float64) float64 {
	if !inRange4x3(a, b, c, d, lo3, hi3) {
		return orient3dBig(a, b, c, d)
	}
	adx, ady, adz := a[0]-d[0], a[1]-d[1], a[2]-d[2]
	bdx, bdy, bdz := b[0]-d[0], b[1]-d[1], b[2]-d[2]
	cdx, cdy, cdz := c[0]-d[0], c[1]-d[1], c[2]-d[2]

	bdxcdy := float64(bdx * cdy)
	cdxbdy := float64(cdx * bdy)
	cdxady := float64(cdx * ady)
	adxcdy := float64(adx * cdy)
	adxbdy := float64(adx * bdy)
	bdxady := float64(bdx * ady)

	det := float64(adz*(bdxcdy-cdxbdy)) + float64(bdz*(cdxady-adxcdy)) + float64(cdz*(adxbdy-bdxady))
	permanent := float64((abs(bdxcdy)+abs(cdxbdy))*abs(adz)) +
		float64((abs(cdxady)+abs(adxcdy))*abs(bdz)) +
		float64((abs(adxbdy)+abs(bdxady))*abs(cdz))
	errbound := float64(o3dErrBoundA * permanent)
	if det > errbound || -det > errbound {
		return det
	}
	return orient3dAdapt(a, b, c, d, permanent)
}

// Orient3DFast is the plain floating-point evaluation of Orient3D's
// determinant. Not robust.
func Orient3DFast(a, b, c, d [3]float64) float64 {
	adx, ady, adz := a[0]-d[0], a[1]-d[1], a[2]-d[2]
	bdx, bdy, bdz := b[0]-d[0], b[1]-d[1], b[2]-d[2]
	cdx, cdy, cdz := c[0]-d[0], c[1]-d[1], c[2]-d[2]
	return float64(adz*(float64(bdx*cdy)-float64(cdx*bdy))) +
		float64(bdz*(float64(cdx*ady)-float64(adx*cdy))) +
		float64(cdz*(float64(adx*bdy)-float64(bdx*ady)))
}

// Orient3DExact evaluates Orient3D's determinant exactly, without the
// adaptive stages.
func Orient3DExact(a, b, c, d [3]float64) float64 {
	if !inRange4x3(a, b, c, d, lo3, hi3) {
		return orient3dBig(a, b, c, d)
	}
	return orient3dExact(a, b, c, d)
}

func rangeAcc3(mx, mn uint64, p [3]float64) (uint64, uint64) {
	mx, mn = rangeAcc(mx, mn, p[0])
	mx, mn = rangeAcc(mx, mn, p[1])
	return rangeAcc(mx, mn, p[2])
}

func inRange4x3(a, b, c, d [3]float64, lo, hi float64) bool {
	mx, mn := rangeAcc3(0, math.MaxUint64, a)
	mx, mn = rangeAcc3(mx, mn, b)
	mx, mn = rangeAcc3(mx, mn, c)
	mx, mn = rangeAcc3(mx, mn, d)
	return rangeOK(mx, mn, lo, hi)
}

func inRange5x3(a, b, c, d, e [3]float64, lo, hi float64) bool {
	mx, mn := rangeAcc3(0, math.MaxUint64, a)
	mx, mn = rangeAcc3(mx, mn, b)
	mx, mn = rangeAcc3(mx, mn, c)
	mx, mn = rangeAcc3(mx, mn, d)
	mx, mn = rangeAcc3(mx, mn, e)
	return rangeOK(mx, mn, lo, hi)
}

// orient3dAdapt runs stages B, C and D.
func orient3dAdapt(a, b, c, d [3]float64, permanent float64) float64 {
	adx, ady, adz := a[0]-d[0], a[1]-d[1], a[2]-d[2]
	bdx, bdy, bdz := b[0]-d[0], b[1]-d[1], b[2]-d[2]
	cdx, cdy, cdz := c[0]-d[0], c[1]-d[1], c[2]-d[2]

	// Stage B: Σ z·(exact 2×2 minor), exactly.
	var t1, t2, t3 [8]float64
	var s1 [16]float64
	var fin [24]float64
	bc := cross2(bdx, cdy, cdx, bdy)
	ca := cross2(cdx, ady, adx, cdy)
	ab := cross2(adx, bdy, bdx, ady)
	e := sumExp(s1[:0], scaleExp(t1[:0], bc[:], adz), scaleExp(t2[:0], ca[:], bdz))
	e = sumExp(fin[:0], e, scaleExp(t3[:0], ab[:], cdz))
	det := estimate(e)
	errbound := o3dErrBoundB * permanent
	if det >= errbound || -det >= errbound {
		return det
	}

	adxt := twoDiffTail(a[0], d[0], adx)
	bdxt := twoDiffTail(b[0], d[0], bdx)
	cdxt := twoDiffTail(c[0], d[0], cdx)
	adyt := twoDiffTail(a[1], d[1], ady)
	bdyt := twoDiffTail(b[1], d[1], bdy)
	cdyt := twoDiffTail(c[1], d[1], cdy)
	adzt := twoDiffTail(a[2], d[2], adz)
	bdzt := twoDiffTail(b[2], d[2], bdz)
	cdzt := twoDiffTail(c[2], d[2], cdz)
	if adxt == 0 && bdxt == 0 && cdxt == 0 && adyt == 0 && bdyt == 0 && cdyt == 0 &&
		adzt == 0 && bdzt == 0 && cdzt == 0 {
		return det
	}

	// Stage C: first-order correction by the tails. One term per row:
	// z·(tail correction of the minor) + ztail·(minor).
	term := func(z, zt, qx, qy, qxt, qyt, rx, ry, rxt, ryt float64) float64 {
		minorT := (float64(qx*ryt) + float64(ry*qxt)) - (float64(qy*rxt) + float64(rx*qyt))
		minor := float64(qx*ry) - float64(qy*rx)
		return float64(z*minorT) + float64(zt*minor)
	}
	errbound = float64(o3dErrBoundC*permanent) + float64(resultErrBound*abs(det))
	det += term(adz, adzt, bdx, bdy, bdxt, bdyt, cdx, cdy, cdxt, cdyt) +
		term(bdz, bdzt, cdx, cdy, cdxt, cdyt, adx, ady, adxt, adyt) +
		term(cdz, cdzt, adx, ady, adxt, adyt, bdx, bdy, bdxt, bdyt)
	if det >= errbound || -det >= errbound {
		return det
	}

	// Stage D: the exact determinant.
	return orient3dExact(a, b, c, d)
}

// orient3dExact evaluates the determinant from the exact differences.
func orient3dExact(a, b, c, d [3]float64) float64 {
	ar := &arena{}
	row := func(p [3]float64) [][]float64 {
		return [][]float64{diffExp(ar, p[0], d[0]), diffExp(ar, p[1], d[1]), diffExp(ar, p[2], d[2])}
	}
	det := detExp(ar, [][][]float64{row(a), row(b), row(c)})
	return det[len(det)-1]
}

func orient3dBig(a, b, c, d [3]float64) float64 {
	if !finite(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], d[0], d[1], d[2]) {
		return nan
	}
	s := scaleOf(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], d[0], d[1], d[2])
	row := func(p [3]float64) []*bigInt {
		return []*bigInt{bigDiff(p[0], d[0], s), bigDiff(p[1], d[1], s), bigDiff(p[2], d[2], s)}
	}
	return bigDet([][]*bigInt{row(a), row(b), row(c)}, 3*s)
}
