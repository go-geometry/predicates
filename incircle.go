package predicates

// InCircle reports whether d lies inside the circle through a, b and c: a
// positive value if it does, a negative value if it lies outside, zero if the
// four points are cocircular. a, b, c must be in counterclockwise order
// (Orient2D(a, b, c) > 0), or the sign is reversed. The sign is exact for all
// finite inputs; the magnitude approximates the determinant
//
//	| ax-dx  ay-dy  (ax-dx)²+(ay-dy)² |
//	| bx-dx  by-dy  (bx-dx)²+(by-dy)² |
//	| cx-dx  cy-dy  (cx-dx)²+(cy-dy)² |
//
// Evaluation is adaptive, as for Orient2D. NaN or ±Inf yields NaN.
func InCircle(a, b, c, d [2]float64) float64 {
	if !inRange4x2(a, b, c, d, lo4, hi4) {
		return incircleBig(a, b, c, d)
	}
	adx, ady := a[0]-d[0], a[1]-d[1]
	bdx, bdy := b[0]-d[0], b[1]-d[1]
	cdx, cdy := c[0]-d[0], c[1]-d[1]

	bdxcdy := float64(bdx * cdy)
	cdxbdy := float64(cdx * bdy)
	alift := float64(adx*adx) + float64(ady*ady)
	cdxady := float64(cdx * ady)
	adxcdy := float64(adx * cdy)
	blift := float64(bdx*bdx) + float64(bdy*bdy)
	adxbdy := float64(adx * bdy)
	bdxady := float64(bdx * ady)
	clift := float64(cdx*cdx) + float64(cdy*cdy)

	det := float64(alift*(bdxcdy-cdxbdy)) + float64(blift*(cdxady-adxcdy)) + float64(clift*(adxbdy-bdxady))
	permanent := float64((abs(bdxcdy)+abs(cdxbdy))*alift) +
		float64((abs(cdxady)+abs(adxcdy))*blift) +
		float64((abs(adxbdy)+abs(bdxady))*clift)
	errbound := float64(iccErrBoundA * permanent)
	if det > errbound || -det > errbound {
		return det
	}
	return incircleAdapt(a, b, c, d, permanent)
}

// InCircleFast is the plain floating-point evaluation of InCircle's
// determinant. Not robust.
func InCircleFast(a, b, c, d [2]float64) float64 {
	adx, ady := a[0]-d[0], a[1]-d[1]
	bdx, bdy := b[0]-d[0], b[1]-d[1]
	cdx, cdy := c[0]-d[0], c[1]-d[1]
	alift := float64(adx*adx) + float64(ady*ady)
	blift := float64(bdx*bdx) + float64(bdy*bdy)
	clift := float64(cdx*cdx) + float64(cdy*cdy)
	return float64(alift*(float64(bdx*cdy)-float64(cdx*bdy))) +
		float64(blift*(float64(cdx*ady)-float64(adx*cdy))) +
		float64(clift*(float64(adx*bdy)-float64(bdx*ady)))
}

// InCircleExact evaluates InCircle's determinant exactly, without the
// adaptive stages.
func InCircleExact(a, b, c, d [2]float64) float64 {
	if !inRange4x2(a, b, c, d, lo4, hi4) {
		return incircleBig(a, b, c, d)
	}
	return incircleExact(a, b, c, d)
}

func inRange4x2(a, b, c, d [2]float64, lo, hi float64) bool {
	return inRange2(a, b, c, lo, hi) && inRange(d[0], lo, hi) && inRange(d[1], lo, hi)
}

// incircleAdapt runs stages B, C and D.
func incircleAdapt(a, b, c, d [2]float64, permanent float64) float64 {
	adx, ady := a[0]-d[0], a[1]-d[1]
	bdx, bdy := b[0]-d[0], b[1]-d[1]
	cdx, cdy := c[0]-d[0], c[1]-d[1]

	// Stage B: the exact determinant of the rounded differences.
	m := [][][]float64{
		{{adx}, {ady}, squareSumExp([]float64{adx}, []float64{ady})},
		{{bdx}, {bdy}, squareSumExp([]float64{bdx}, []float64{bdy})},
		{{cdx}, {cdy}, squareSumExp([]float64{cdx}, []float64{cdy})},
	}
	det := estimate(detExp(m))
	errbound := iccErrBoundB * permanent
	if det >= errbound || -det >= errbound {
		return det
	}

	adxt := twoDiffTail(a[0], d[0], adx)
	adyt := twoDiffTail(a[1], d[1], ady)
	bdxt := twoDiffTail(b[0], d[0], bdx)
	bdyt := twoDiffTail(b[1], d[1], bdy)
	cdxt := twoDiffTail(c[0], d[0], cdx)
	cdyt := twoDiffTail(c[1], d[1], cdy)
	if adxt == 0 && bdxt == 0 && cdxt == 0 && adyt == 0 && bdyt == 0 && cdyt == 0 {
		return det
	}

	// Stage C: first-order correction by the tails.
	errbound = float64(iccErrBoundC*permanent) + float64(resultErrBound*abs(det))
	// One term per row: lift·(tail correction of the 2×2 minor) +
	// 2·(tail correction of the lift)·(minor).
	term := func(px, py, pxt, pyt, qx, qy, qxt, qyt, rx, ry, rxt, ryt float64) float64 {
		lift := float64(px*px) + float64(py*py)
		minorT := (float64(qx*ryt) + float64(ry*qxt)) - (float64(qy*rxt) + float64(rx*qyt))
		liftT := float64(px*pxt) + float64(py*pyt)
		minor := float64(qx*ry) - float64(qy*rx)
		return float64(lift*minorT) + float64(float64(2*liftT)*minor)
	}
	det += term(adx, ady, adxt, adyt, bdx, bdy, bdxt, bdyt, cdx, cdy, cdxt, cdyt) +
		term(bdx, bdy, bdxt, bdyt, cdx, cdy, cdxt, cdyt, adx, ady, adxt, adyt) +
		term(cdx, cdy, cdxt, cdyt, adx, ady, adxt, adyt, bdx, bdy, bdxt, bdyt)
	if det >= errbound || -det >= errbound {
		return det
	}

	// Stage D: the exact determinant.
	return incircleExact(a, b, c, d)
}

// incircleExact evaluates the determinant from the exact differences.
func incircleExact(a, b, c, d [2]float64) float64 {
	row := func(p [2]float64) [][]float64 {
		x, y := diffExp(p[0], d[0]), diffExp(p[1], d[1])
		return [][]float64{x, y, squareSumExp(x, y)}
	}
	det := detExp([][][]float64{row(a), row(b), row(c)})
	return det[len(det)-1]
}

func incircleBig(a, b, c, d [2]float64) float64 {
	if !finite(a[0], a[1], b[0], b[1], c[0], c[1], d[0], d[1]) {
		return nan
	}
	row := func(p [2]float64) []*bigInt {
		x, y := bigDiff(p[0], d[0]), bigDiff(p[1], d[1])
		return []*bigInt{x, y, bigSquareSum(x, y)}
	}
	return bigDet([][]*bigInt{row(a), row(b), row(c)}, 4)
}
