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
template <class T> inline T v2p_T(const T &x, const T &d, const T &k, double eps = EPS_IN) {
  using std::exp; using std::log;
  T t;
  if (val(x) <= val(d)) t = 1.0 - x - (1.0 - d) * exp(k * log(x / d));
  else t = x - d * exp(log((1.0 - x) / (1.0 - d)) / k);
  const double v = val(t);
  return v < eps ? T(eps) : (v > 1.0 - eps ? T(1.0 - eps) : t);   // zero derivative where clamped
}

// log density of the base copula (the pair copula `bfam` with rotation `brot` and parameters bp, bp2) at (a, b)
template <class T> inline T nl_base_ld(int bfam, int brot, const T &bp, const T &bp2, const T &a, const T &b) {
  T lc, h1, h2;
  pair_eval(bfam, brot, bp, bp2, a, b, false, false, lc, h1, h2);
  return lc;
}

// Generic form: the base copula enters through ld(a, b), its log density at the base-scale values (a, b), and
// hb(a, b, h1, h2), its h-functions h1 = P(V2 <= b | V1 = a) and h2 = P(V1 <= a | V2 = b).
//
// P(X <= x | Y = g) for the IVT pair with transform (dx, kx) of X and a base-scale value ty = T_y(g) of the
// conditioning variable; x_first is true when X is the earlier variable of the pair.
template <class T, class LD, class HB>
inline T nl_cdf_g(const T &dx, const T &kx, const T &x, const T &ty, bool x_first, const LD &ld, const HB &hb,
                  double eps = EPS_IN) {
  using std::exp;
  if (val(kx) == 1.0 && dzero(kx)) {
    // linear transform of the integrated variable (and no kappa derivative needed): closed form
    const T tx = vlin(x, dx);
    T h1, h2;
    if (x_first) { hb(tx, ty, h1, h2); return clampu(dx - alin(x, dx) * h2); }
    hb(ty, tx, h1, h2);
    return clampu(dx - alin(x, dx) * h1);
  }
  const NlRule &r = nl_rule();
  const T m1 = val(x) <= val(dx) ? x : dx;   // lower panel [0, min(x, dx)]
  T sum(0.0);
  for (int k = 0; k < NL_NODES; k++) {
    const T s = v2p_T(m1 * r.yq[k], dx, kx, eps);
    const T lc = x_first ? ld(s, ty) : ld(ty, s);
    sum = sum + exp(lc) * (m1 * r.wq[k]);
  }
  if (val(x) > val(dx)) {                    // upper panel [dx, x]
    const T len = x - dx;
    for (int k = 0; k < NL_NODES; k++) {
      const T s = v2p_T(dx + len * r.zq[k], dx, kx, eps);
      const T lc = x_first ? ld(s, ty) : ld(ty, s);
      sum = sum + exp(lc) * (len * r.vq[k]);
    }
  }
  return clampu(sum);
}

// pair-copula bases (Clayton, Gumbel, Joe, Frank, BB1 with rotations)
template <class T>
inline T nl_cdf(int bfam, int brot, const T &bp, const T &bp2, const T &dx, const T &kx, const T &x, const T &ty, bool x_first) {
  return nl_cdf_g(dx, kx, x, ty, x_first,
                  [&](const T &a, const T &b) { return nl_base_ld(bfam, brot, bp, bp2, a, b); },
                  [&](const T &a, const T &b, T &h1, T &h2) { T lc; pair_eval(bfam, brot, bp, bp2, a, b, true, true, lc, h1, h2); });
}

// ---- absolute spherical t (ast) base: quantile table -----------------------------------------------------------
//
// The ast density needs the t quantile z(V) = qt_upper((1 - V) / 2, nu) of every base-scale value, i.e. of every
// quadrature node. For one nu the function is tabulated once on a uniform grid of r = -log(1 - V) as
// psi(r) = log(1 + z) (smooth, close to r / nu in the tail) with its exact r-derivative, and evaluated by cubic
// Hermite interpolation; for dual numbers the derivative with respect to nu is tabulated too.

inline double nl_pt_dk(double z, double k) {   // d pt_upper(z, k) / dk by a central difference
  const double h = 1e-4 * k;
  return (R::pt(z, k + h, 0, 0) - R::pt(z, k - h, 0, 0)) / (2.0 * h);
}

struct AstZTab {
  static const int N = 640;
  static constexpr double RMAX = 28.0;
  double h = RMAX / N, nu = 0.0;
  bool withNu = false;
  double psi[N + 1], dpsi[N + 1], psin[N + 1], dpsin[N + 1];
  void build(double nu_, bool withNu_) {
    nu = nu_; withNu = withNu_;
    for (int k = 0; k <= N; k++) {
      const double e = std::exp(-k * h);
      const double z = k == 0 ? 0.0 : R::qt(0.5 * e, nu, 0, 0);
      const double ft = R::dt(z, nu, 0);
      psi[k] = std::log1p(z);
      dpsi[k] = e / (2.0 * ft) / (1.0 + z);                     // d psi / dr
      if (withNu) psin[k] = (nl_pt_dk(z, nu) / ft) / (1.0 + z);   // d psi / d nu
    }
    if (withNu) {
      for (int k = 0; k <= N; k++) {
        const int a = k == 0 ? 0 : k - 1, b = k == N ? N : k + 1;
        dpsin[k] = (psin[b] - psin[a]) / ((b - a) * h);
      }
    }
  }
  // z, dz/dV and dz/dnu at the base-scale value v
  inline void eval(double v, double &z, double &dzdv, double &dzdnu) const {
    v = v < 1e-12 ? 1e-12 : (v > 1.0 - 1e-12 ? 1.0 - 1e-12 : v);
    const double r = -std::log1p(-v);
    const double sc = r / h;
    int k = static_cast<int>(sc);
    if (k >= N) k = N - 1;
    const double t = sc - k, t2 = t * t, t3 = t2 * t;
    const double h00 = 2 * t3 - 3 * t2 + 1, h10 = t3 - 2 * t2 + t, h01 = -2 * t3 + 3 * t2, h11 = t3 - t2;
    const double p = h00 * psi[k] + h10 * h * dpsi[k] + h01 * psi[k + 1] + h11 * h * dpsi[k + 1];
    const double dp = ((6 * t2 - 6 * t) * psi[k] + (3 * t2 - 4 * t + 1) * h * dpsi[k] + (-6 * t2 + 6 * t) * psi[k + 1] +
                       (3 * t2 - 2 * t) * h * dpsi[k + 1]) / h;
    z = std::expm1(p);
    dzdv = (1.0 + z) * dp / (1.0 - v);
    dzdnu = 0.0;
    if (withNu) {
      const double pn = h00 * psin[k] + h10 * h * dpsin[k] + h01 * psin[k + 1] + h11 * h * dpsin[k + 1];
      dzdnu = (1.0 + z) * pn;
    }
  }
};

// t quantile of a base-scale value, with its derivatives for dual numbers
inline double ast_z(const AstZTab &t, double v, double) { double z, a, b; t.eval(v, z, a, b); return z; }
template <int P> inline Dual<P> ast_z(const AstZTab &t, const Dual<P> &v, const Dual<P> &nu) {
  double z, dv, dn;
  t.eval(v.v, z, dv, dn);
  Dual<P> r;
  r.v = z;
  for (int k = 0; k < P; k++) r.d[k] = dv * v.d[k] + dn * nu.d[k];
  return r;
}

// log density of the ast copula from the t quantiles (zx, zy) of the two arguments
template <class T> inline T ast_ld_z(const T &zx, const T &zy, const T &nu, const T &c0) {
  using std::log1p;
  const T x2 = zx * zx / nu, y2 = zy * zy / nu;
  return c0 + (-(nu + 2.0) / 2.0) * log1p(x2 + y2) + ((nu + 1.0) / 2.0) * (log1p(x2) + log1p(y2));
}

// u-scale value 1 - 2 pt_upper(z, nu + 1), clamped (the h-function of the ast copula from the scaled quantile z)
inline double nl_ast_exit(double z, double nu) { return clampu(1.0 - 2.0 * R::pt(z, nu + 1.0, 0, 0)); }
template <int P> inline Dual<P> nl_ast_exit(const Dual<P> &z, const Dual<P> &nu) {
  const double u = 1.0 - 2.0 * R::pt(z.v, nu.v + 1.0, 0, 0);
  Dual<P> r;
  if (u < EPS) { r.v = EPS; return r; }
  if (u > 1.0 - EPS) { r.v = 1.0 - EPS; return r; }
  const double dz = 2.0 * R::dt(z.v, nu.v + 1.0, 0), dn = -2.0 * nl_pt_dk(z.v, nu.v + 1.0);
  r.v = u;
  for (int k = 0; k < P; k++) r.d[k] = dz * z.d[k] + dn * nu.d[k];
  return r;
}

// Quadrature policy of an ast base: log density and h-functions on the base scale through the quantile table
template <class T> struct AstPolicy {
  const AstZTab *tab;
  T nu, c0;
  inline T ld(const T &a, const T &b) const { return ast_ld_z(ast_z(*tab, a, nu), ast_z(*tab, b, nu), nu, c0); }
  inline void hb(const T &a, const T &b, T &h1, T &h2) const {
    using std::sqrt;
    const T za = ast_z(*tab, a, nu), zb = ast_z(*tab, b, nu);
    const T k1 = nu + 1.0;
    h1 = nl_ast_exit(zb * sqrt(k1 / (nu + za * za)), nu);   // P(V2 <= b | V1 = a)
    h2 = nl_ast_exit(za * sqrt(k1 / (nu + zb * zb)), nu);   // P(V1 <= a | V2 = b)
  }
};

}  // namespace
