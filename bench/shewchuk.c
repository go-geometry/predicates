// predicates.c (public domain, J. R. Shewchuk), included verbatim from
// predicates_c.inc, compiled without contraction into fused multiply-adds:
// cgo refuses -ffp-contract=off in #cgo CFLAGS, so the pragmas say it here.
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif
#include "predicates_c.inc"
