package predicates

// InSphere reports whether e lies inside the sphere through a, b, c and d: a
// positive value if it does, a negative value if it lies outside, zero if the
// five points are cospherical. a, b, c, d must have positive orientation
// (Orient3D(a, b, c, d) > 0), or the sign is reversed. The sign is exact for
// all finite inputs; the magnitude approximates the determinant whose rows
// are [p-e, |p-e|²] for p = a, b, c, d.
//
// Evaluation is adaptive, as for Orient2D. NaN or ±Inf yields NaN.
func InSphere(a, b, c, d, e [3]float64) float64 {
	if !inRange4x3(a, b, c, d, lo5, hi5) || !inRange3(e, lo5, hi5) {
		return insphereBig(a, b, c, d, e)
	}
	aex, aey, aez := a[0]-e[0], a[1]-e[1], a[2]-e[2]
	bex, bey, bez := b[0]-e[0], b[1]-e[1], b[2]-e[2]
	cex, cey, cez := c[0]-e[0], c[1]-e[1], c[2]-e[2]
	dex, dey, dez := d[0]-e[0], d[1]-e[1], d[2]-e[2]

	aexbey, bexaey := float64(aex*bey), float64(bex*aey)
	bexcey, cexbey := float64(bex*cey), float64(cex*bey)
	cexdey, dexcey := float64(cex*dey), float64(dex*cey)
	dexaey, aexdey := float64(dex*aey), float64(aex*dey)
	aexcey, cexaey := float64(aex*cey), float64(cex*aey)
	bexdey, dexbey := float64(bex*dey), float64(dex*bey)
	ab := aexbey - bexaey
	bc := bexcey - cexbey
	cd := cexdey - dexcey
	da := dexaey - aexdey
	ac := aexcey - cexaey
	bd := bexdey - dexbey

	abc := float64(aez*bc) - float64(bez*ac) + float64(cez*ab)
	bcd := float64(bez*cd) - float64(cez*bd) + float64(dez*bc)
	cda := float64(cez*da) + float64(dez*ac) + float64(aez*cd)
	dab := float64(dez*ab) + float64(aez*bd) + float64(bez*da)

	alift := lift3(aex, aey, aez)
	blift := lift3(bex, bey, bez)
	clift := lift3(cex, cey, cez)
	dlift := lift3(dex, dey, dez)

	det := (float64(dlift*abc) - float64(clift*dab)) + (float64(blift*cda) - float64(alift*bcd))

	aezp, bezp, cezp, dezp := abs(aez), abs(bez), abs(cez), abs(dez)
	abp := abs(aexbey) + abs(bexaey)
	bcp := abs(bexcey) + abs(cexbey)
	cdp := abs(cexdey) + abs(dexcey)
	dap := abs(dexaey) + abs(aexdey)
	acp := abs(aexcey) + abs(cexaey)
	bdp := abs(bexdey) + abs(dexbey)
	permanent := float64((float64(cdp*bezp)+float64(bdp*cezp)+float64(bcp*dezp))*alift) +
		float64((float64(dap*cezp)+float64(acp*dezp)+float64(cdp*aezp))*blift) +
		float64((float64(abp*dezp)+float64(bdp*aezp)+float64(dap*bezp))*clift) +
		float64((float64(bcp*aezp)+float64(acp*bezp)+float64(abp*cezp))*dlift)
	errbound := float64(ispErrBoundA * permanent)
	if det > errbound || -det > errbound {
		return det
	}
	return insphereAdapt(a, b, c, d, e, permanent)
}

// lift3 is the rounded squared norm x²+y²+z².
func lift3(x, y, z float64) float64 {
	return float64(x*x) + float64(y*y) + float64(z*z)
}

// InSphereFast is the plain floating-point evaluation of InSphere's
// determinant. Not robust.
func InSphereFast(a, b, c, d, e [3]float64) float64 {
	aex, aey, aez := a[0]-e[0], a[1]-e[1], a[2]-e[2]
	bex, bey, bez := b[0]-e[0], b[1]-e[1], b[2]-e[2]
	cex, cey, cez := c[0]-e[0], c[1]-e[1], c[2]-e[2]
	dex, dey, dez := d[0]-e[0], d[1]-e[1], d[2]-e[2]
	ab := float64(aex*bey) - float64(bex*aey)
	bc := float64(bex*cey) - float64(cex*bey)
	cd := float64(cex*dey) - float64(dex*cey)
	da := float64(dex*aey) - float64(aex*dey)
	ac := float64(aex*cey) - float64(cex*aey)
	bd := float64(bex*dey) - float64(dex*bey)
	abc := float64(aez*bc) - float64(bez*ac) + float64(cez*ab)
	bcd := float64(bez*cd) - float64(cez*bd) + float64(dez*bc)
	cda := float64(cez*da) + float64(dez*ac) + float64(aez*cd)
	dab := float64(dez*ab) + float64(aez*bd) + float64(bez*da)
	return (float64(lift3(dex, dey, dez)*abc) - float64(lift3(cex, cey, cez)*dab)) +
		(float64(lift3(bex, bey, bez)*cda) - float64(lift3(aex, aey, aez)*bcd))
}

// InSphereExact evaluates InSphere's determinant exactly, without the
// adaptive stages.
func InSphereExact(a, b, c, d, e [3]float64) float64 {
	if !inRange4x3(a, b, c, d, lo5, hi5) || !inRange3(e, lo5, hi5) {
		return insphereBig(a, b, c, d, e)
	}
	return insphereExact(a, b, c, d, e)
}

// insphereAdapt runs stages B, C and D.
func insphereAdapt(a, b, c, d, e [3]float64, permanent float64) float64 {
	aex, aey, aez := a[0]-e[0], a[1]-e[1], a[2]-e[2]
	bex, bey, bez := b[0]-e[0], b[1]-e[1], b[2]-e[2]
	cex, cey, cez := c[0]-e[0], c[1]-e[1], c[2]-e[2]
	dex, dey, dez := d[0]-e[0], d[1]-e[1], d[2]-e[2]

	// Stage B: the exact determinant of the rounded differences.
	row := func(x, y, z float64) [][]float64 {
		return [][]float64{{x}, {y}, {z}, squareSumExp([]float64{x}, []float64{y}, []float64{z})}
	}
	det := estimate(detExp([][][]float64{row(aex, aey, aez), row(bex, bey, bez), row(cex, cey, cez), row(dex, dey, dez)}))
	errbound := ispErrBoundB * permanent
	if det >= errbound || -det >= errbound {
		return det
	}

	aext, aeyt, aezt := twoDiffTail(a[0], e[0], aex), twoDiffTail(a[1], e[1], aey), twoDiffTail(a[2], e[2], aez)
	bext, beyt, bezt := twoDiffTail(b[0], e[0], bex), twoDiffTail(b[1], e[1], bey), twoDiffTail(b[2], e[2], bez)
	cext, ceyt, cezt := twoDiffTail(c[0], e[0], cex), twoDiffTail(c[1], e[1], cey), twoDiffTail(c[2], e[2], cez)
	dext, deyt, dezt := twoDiffTail(d[0], e[0], dex), twoDiffTail(d[1], e[1], dey), twoDiffTail(d[2], e[2], dez)
	if aext == 0 && aeyt == 0 && aezt == 0 && bext == 0 && beyt == 0 && bezt == 0 &&
		cext == 0 && ceyt == 0 && cezt == 0 && dext == 0 && deyt == 0 && dezt == 0 {
		return det
	}

	// Stage C: first-order correction by the tails. top(q, r) is the leading
	// component of the exact 2×2 minor qx·ry - rx·qy; eps(q, r) its tail
	// correction.
	top := func(qx, qy, rx, ry float64) float64 { return cross2(qx, ry, rx, qy)[3] }
	eps := func(qx, qy, qxt, qyt, rx, ry, rxt, ryt float64) float64 {
		return (float64(qx*ryt) + float64(ry*qxt)) - (float64(qy*rxt) + float64(rx*qyt))
	}
	ab3, bc3, cd3 := top(aex, aey, bex, bey), top(bex, bey, cex, cey), top(cex, cey, dex, dey)
	da3, ac3, bd3 := top(dex, dey, aex, aey), top(aex, aey, cex, cey), top(bex, bey, dex, dey)
	abeps := eps(aex, aey, aext, aeyt, bex, bey, bext, beyt)
	bceps := eps(bex, bey, bext, beyt, cex, cey, cext, ceyt)
	cdeps := eps(cex, cey, cext, ceyt, dex, dey, dext, deyt)
	daeps := eps(dex, dey, dext, deyt, aex, aey, aext, aeyt)
	aceps := eps(aex, aey, aext, aeyt, cex, cey, cext, ceyt)
	bdeps := eps(bex, bey, bext, beyt, dex, dey, dext, deyt)
	// Sums of three products, each product rounded.
	s3 := func(p, q, r, s, t, u float64) float64 { return float64(p*q) + float64(r*s) + float64(t*u) }
	dot := func(x, xt, y, yt, z, zt float64) float64 { return s3(x, xt, y, yt, z, zt) }

	errbound = float64(ispErrBoundC*permanent) + float64(resultErrBound*abs(det))
	first := (float64(lift3(bex, bey, bez)*(s3(cez, daeps, dez, aceps, aez, cdeps)+s3(cezt, da3, dezt, ac3, aezt, cd3))) +
		float64(lift3(dex, dey, dez)*((float64(aez*bceps)-float64(bez*aceps)+float64(cez*abeps))+
			(float64(aezt*bc3)-float64(bezt*ac3)+float64(cezt*ab3))))) -
		(float64(lift3(aex, aey, aez)*((float64(bez*cdeps)-float64(cez*bdeps)+float64(dez*bceps))+
			(float64(bezt*cd3)-float64(cezt*bd3)+float64(dezt*bc3)))) +
			float64(lift3(cex, cey, cez)*(s3(dez, abeps, aez, bdeps, bez, daeps)+s3(dezt, ab3, aezt, bd3, bezt, da3))))
	second := (float64(dot(bex, bext, bey, beyt, bez, bezt)*s3(cez, da3, dez, ac3, aez, cd3)) +
		float64(dot(dex, dext, dey, deyt, dez, dezt)*(float64(aez*bc3)-float64(bez*ac3)+float64(cez*ab3)))) -
		(float64(dot(aex, aext, aey, aeyt, aez, aezt)*(float64(bez*cd3)-float64(cez*bd3)+float64(dez*bc3))) +
			float64(dot(cex, cext, cey, ceyt, cez, cezt)*s3(dez, ab3, aez, bd3, bez, da3)))
	det += first + float64(2*second)
	if det >= errbound || -det >= errbound {
		return det
	}

	// Stage D: the exact determinant.
	return insphereExact(a, b, c, d, e)
}

// insphereExact evaluates the determinant from the exact differences.
func insphereExact(a, b, c, d, e [3]float64) float64 {
	row := func(p [3]float64) [][]float64 {
		x, y, z := diffExp(p[0], e[0]), diffExp(p[1], e[1]), diffExp(p[2], e[2])
		return [][]float64{x, y, z, squareSumExp(x, y, z)}
	}
	det := detExp([][][]float64{row(a), row(b), row(c), row(d)})
	return det[len(det)-1]
}

func insphereBig(a, b, c, d, e [3]float64) float64 {
	if !finite(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], d[0], d[1], d[2], e[0], e[1], e[2]) {
		return nan
	}
	row := func(p [3]float64) []*bigInt {
		x, y, z := bigDiff(p[0], e[0]), bigDiff(p[1], e[1]), bigDiff(p[2], e[2])
		return []*bigInt{x, y, z, bigSquareSum(x, y, z)}
	}
	return bigDet([][]*bigInt{row(a), row(b), row(c), row(d)}, 5)
}
