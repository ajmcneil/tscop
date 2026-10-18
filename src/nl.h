// Nonlinear (v2p) IVT pair copulas for the d-vine kernel: the v2p transform, the pair log density
// c_V(T1(a), T2(b)) and the conditional distribution functions as graded Gauss-Legendre integrals of that
// density over the u scale (24 nodes, grading 3). Mirrors the R reference R/dvine-nl.R; see
// docs/sdvinecop-nonlinear-vt-spec.md. Templates over the scalar type (double or a dual number).
#pragma once
#include "pair.h"
namespace {

const int NL_NODES = 24;
const double NL_GRADE = 3.0;

// Gauss-Legendre nodes and weights on (0, 1), and the graded maps of the two panels:
//   lower panel [0, m]:  x = m yq,  dx = m wq   (nodes clustered at 0)
//   upper panel [d, x]:  x = d + len zq, dx = len vq   (nodes clustered at x)
struct NlRule {
  double yq[NL_NODES], wq[NL_NODES], zq[NL_NODES], vq[NL_NODES];
  NlRule() {
    const int n = NL_NODES;
    double xs[NL_NODES], ws[NL_NODES];
    for (int i = 0; i < (n + 1) / 2; i++) {
      double z = std::cos(M_PI * (i + 0.75) / (n + 0.5)), pp = 0.0;
      for (int it = 0; it < 100; it++) {
        double p1 = 1.0, p2 = 0.0;
        for (int j = 1; j <= n; j++) { const double p3 = p2; p2 = p1; p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j; }
        pp = n * (z * p1 - p2) / (z * z - 1.0);
        const double z1 = z;
        z = z1 - p1 / pp;
        if (std::fabs(z - z1) < 1e-15) break;
      }
      xs[i] = -z; xs[n - 1 - i] = z;
      ws[i] = ws[n - 1 - i] = 2.0 / ((1.0 - z * z) * pp * pp);
    }
    for (int k = 0; k < n; k++) {
      const double y = 0.5 * (xs[k] + 1.0), w = 0.5 * ws[k];
      yq[k] = std::pow(y, NL_GRADE);
      wq[k] = NL_GRADE * std::pow(y, NL_GRADE - 1.0) * w;
      zq[k] = 1.0 - std::pow(1.0 - y, NL_GRADE);
      vq[k] = NL_GRADE * std::pow(1.0 - y, NL_GRADE - 1.0) * w;
    }
  }
};
inline const NlRule &nl_rule() { static const NlRule r; return r; }

// no derivative direction active (always true for doubles)
inline bool dzero(double) { return true; }
template <int P> inline bool dzero(const Dual<P> &x) { for (int k = 0; k < P; k++) if (x.d[k] != 0.0) return false; return true; }

// v2p transform with fulcrum d and shape k, clamped to the open unit interval (zero derivative where clamped)
template <class T> inline T v2p_T(const T &x, const T &d, const T &k) {
  using std::exp; using std::log;
  T t;
  if (val(x) <= val(d)) t = 1.0 - x - (1.0 - d) * exp(k * log(x / d));
  else t = x - d * exp(log((1.0 - x) / (1.0 - d)) / k);
  return clampin(t);
}

// log density of the base copula (the pair copula `bfam` with rotation `brot` and parameters bp, bp2) at (a, b)
template <class T> inline T nl_base_ld(int bfam, int brot, const T &bp, const T &bp2, const T &a, const T &b) {
  T lc, h1, h2;
  pair_eval(bfam, brot, bp, bp2, a, b, false, false, lc, h1, h2);
  return lc;
}

// P(X <= x | Y = g) for the IVT pair with transform (dx, kx) of X and a base-scale value ty = T_y(g) of the
// conditioning variable; x_first is true when X is the earlier variable of the pair.
template <class T>
inline T nl_cdf(int bfam, int brot, const T &bp, const T &bp2, const T &dx, const T &kx, const T &x, const T &ty, bool x_first) {
  using std::exp;
  if (val(kx) == 1.0 && dzero(kx)) {
    // linear transform of the integrated variable (and no kappa derivative needed): closed form
    const T tx = vlin(x, dx);
    T lc, h1, h2;
    if (x_first) { pair_eval(bfam, brot, bp, bp2, tx, ty, true, true, lc, h1, h2); return clampu(dx - alin(x, dx) * h2); }
    pair_eval(bfam, brot, bp, bp2, ty, tx, true, true, lc, h1, h2);
    return clampu(dx - alin(x, dx) * h1);
  }
  const NlRule &r = nl_rule();
  const T m1 = val(x) <= val(dx) ? x : dx;   // lower panel [0, min(x, dx)]
  T sum(0.0);
  for (int k = 0; k < NL_NODES; k++) {
    const T s = v2p_T(m1 * r.yq[k], dx, kx);
    const T lc = x_first ? nl_base_ld(bfam, brot, bp, bp2, s, ty) : nl_base_ld(bfam, brot, bp, bp2, ty, s);
    sum = sum + exp(lc) * (m1 * r.wq[k]);
  }
  if (val(x) > val(dx)) {                    // upper panel [dx, x]
    const T len = x - dx;
    for (int k = 0; k < NL_NODES; k++) {
      const T s = v2p_T(dx + len * r.zq[k], dx, kx);
      const T lc = x_first ? nl_base_ld(bfam, brot, bp, bp2, s, ty) : nl_base_ld(bfam, brot, bp, bp2, ty, s);
      sum = sum + exp(lc) * (len * r.vq[k]);
    }
  }
  return clampu(sum);
}

}  // namespace
