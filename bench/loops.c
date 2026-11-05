// Timing loops over Shewchuk's predicates.c, so that one cgo call covers many
// predicate calls and the cgo transition does not enter the measurement.

double orient2d(double *pa, double *pb, double *pc);
double incircle(double *pa, double *pb, double *pc, double *pd);
double orient3d(double *pa, double *pb, double *pc, double *pd);
double insphere(double *pa, double *pb, double *pc, double *pd, double *pe);
void exactinit(void);

// pts holds n records of 4 points × 2 coordinates (2-D) or 5 × 3 (3-D).
double loop_orient2d(double *pts, int n, long iters) {
	double s = 0;
	for (long i = 0; i < iters; i++) {
		double *p = pts + 8 * (i % n);
		s += orient2d(p, p + 2, p + 4);
	}
	return s;
}

double loop_incircle(double *pts, int n, long iters) {
	double s = 0;
	for (long i = 0; i < iters; i++) {
		double *p = pts + 8 * (i % n);
		s += incircle(p, p + 2, p + 4, p + 6);
	}
	return s;
}

double loop_orient3d(double *pts, int n, long iters) {
	double s = 0;
	for (long i = 0; i < iters; i++) {
		double *p = pts + 15 * (i % n);
		s += orient3d(p, p + 3, p + 6, p + 9);
	}
	return s;
}

double loop_insphere(double *pts, int n, long iters) {
	double s = 0;
	for (long i = 0; i < iters; i++) {
		double *p = pts + 15 * (i % n);
		s += insphere(p, p + 3, p + 6, p + 9, p + 12);
	}
	return s;
}
