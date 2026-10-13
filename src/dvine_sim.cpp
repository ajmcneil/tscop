// C++ kernels for the d-vine simulation (exact stationary simulation by
// inverting the forward chain) and the forward Rosenblatt transform of tied
// sdvinecop models, with the same pair copulas as the likelihood kernel
// (independence, Gauss, Clayton, Gumbel, Joe, Frank with rotations, and IVT
// pairs with these bases). They follow dv_simulate_r and dv_rosenblatt_r
// (R/dvine-engine.R) step by step.
#include "pair.h"

namespace {

struct Level {
  int fam, rot, bfam, brot;
  bool ivt;
  double th, bth, d1, d2;
};

inline double pn(double x) { return R::pnorm(x, 0, 1, 1, 0); }
inline double qn(double p) { return R::qnorm(p, 0, 1, 1, 0); }

// h-function of a (non-IVT) pair copula: cv = 1: P(B <= b | A = a), cv = 2: P(A <= a | B = b)
inline double h_plain(int fam, int rot, double th, double a, double b, int cv) {
  if (fam == GAUSS) {
    const double s = std::sqrt(1.0 - th * th);
    return cv == 1 ? pn((qn(clampin(b)) - th * qn(clampin(a))) / s)
                   : pn((qn(clampin(a)) - th * qn(clampin(b))) / s);
  }
  double lc, h1, h2;
  pair_eval(fam, rot, th, a, b, true, false, lc, h1, h2);
  return cv == 1 ? h1 : h2;
}

// inverse of h(a, .) (cv = 1) of a rotation 0 copula, by closed form or safeguarded Newton
inline double hinv0(int fam, double th, double a, double p) {
  a = clampin(a);
  p = clampin(p);
  switch (fam) {
    case INDEP: return p;
    case GAUSS: return pn(th * qn(a) + std::sqrt(1.0 - th * th) * qn(p));
    case CLAYTON: {
      const double la = std::log(a);
      const double g = -th / (1.0 + th) * (std::log(p) + (th + 1.0) * la);
      const double q = std::expm1(g) - std::expm1(-th * la);   // b^-th - 1
      return std::exp(-std::log1p(q) / th);
    }
    case FRANK: {
      const double Pa = std::expm1(-th * a), E = std::expm1(-th);
      const double q = p * E / (Pa + 1.0 - p * Pa);              // exp(-th b) - 1
      return -std::log1p(q) / th;
    }
    default: break;
  }
  // Gumbel and Joe: h1(a, b) increasing in b with derivative c(a, b)
  double lo = EPS_IN, hi = 1.0 - EPS_IN, b = p;
  for (int it = 0; it < 100; it++) {
    double lc, h1, h2;
    if (fam == GUMBEL) gumbel0(a, b, th, lc, h1, h2); else joe0(a, b, th, lc, h1, h2);
    const double f = h1 - p;
    if (std::fabs(f) < 1e-15) break;                       // converged: keep b
    if (f > 0) hi = b; else lo = b;
    double bn = b - f / std::exp(lc);
    if (!std::isfinite(bn) || bn <= lo || bn >= hi) bn = 0.5 * (lo + hi);
    if (std::fabs(bn - b) < 1e-15 * std::max(1.0, std::fabs(bn))) { b = bn; break; }
    b = bn;
  }
  return b;
}

inline double hinv_plain(int fam, int rot, double th, double a, double p) {
  switch (rot) {
    case 180: return 1.0 - hinv0(fam, th, 1.0 - a, 1.0 - p);
    case 90: return hinv0(fam, th, 1.0 - a, p);
    case 270: return 1.0 - hinv0(fam, th, a, 1.0 - p);
    default: return hinv0(fam, th, a, p);
  }
}

inline double vl(double x, double d) { return vlin(x, d); }
inline double al(double x, double d) { return alin(x, d); }

inline double lv_h(const Level &L, double a, double b, int cv) {
  if (!L.ivt) return clampu(h_plain(L.fam, L.rot, L.th, a, b, cv));
  const double x = vl(a, L.d1), y = vl(b, L.d2);
  const double hb = clampu(h_plain(L.bfam, L.brot, L.bth, x, y, cv));
  return clampu(cv == 1 ? L.d2 - al(b, L.d2) * hb : L.d1 - al(a, L.d1) * hb);
}

inline double lv_hinv(const Level &L, double a, double p) {
  if (!L.ivt) return clampu(hinv_plain(L.fam, L.rot, L.th, a, p));
  const double hb = clampu(hinv_plain(L.bfam, L.brot, L.bth, vl(a, L.d1), vl(p, L.d2)));
  return clampu(L.d2 - al(p, L.d2) * hb);
}

std::vector<Level> make_levels(const IntegerVector &fam, const IntegerVector &rot, const NumericVector &par,
                               const IntegerVector &ivt, const IntegerVector &bfam, const IntegerVector &brot,
                               const NumericVector &bpar, double d1, double d2) {
  const int K = fam.size();
  std::vector<Level> lv(K);
  for (int j = 0; j < K; j++) {
    lv[j].fam = fam[j]; lv[j].rot = rot[j]; lv[j].th = par[j];
    lv[j].ivt = ivt[j] != 0; lv[j].bfam = bfam[j]; lv[j].brot = brot[j]; lv[j].bth = bpar[j];
    lv[j].d1 = d1; lv[j].d2 = d2;
  }
  return lv;
}

}  // namespace

// Exact stationary simulation from independent U(0, 1) innovations w (n x nrep).
// [[Rcpp::export]]
NumericMatrix dv_simulate_cpp(NumericMatrix w, IntegerVector fam, IntegerVector rot, NumericVector par,
                              IntegerVector ivt, IntegerVector bfam, IntegerVector brot, NumericVector bpar,
                              double d1, double d2) {
  const int n = w.nrow(), nrep = w.ncol(), K = fam.size();
  const std::vector<Level> lv = make_levels(fam, rot, par, ivt, bfam, brot, bpar, d1, d2);
  NumericMatrix u(n, nrep);
  std::vector<double> Bprev(K + 1), Bnew(K + 1), Fm(K + 1);
  for (int r = 0; r < nrep; r++) {
    for (int t = 0; t < n; t++) {
      const int k = std::min(t, K);
      double x = w(t, r);
      for (int j = k; j >= 1; j--) {
        x = lv_hinv(lv[j - 1], Bprev[j - 1], x);
        Fm[j - 1] = x;
      }
      u(t, r) = x;
      const int nb = std::min(k, K - 1);
      Bnew[0] = x;
      for (int i = 1; i <= nb; i++) Bnew[i] = lv_h(lv[i - 1], Bprev[i - 1], Fm[i - 1], 2);
      for (int i = 0; i <= nb; i++) Bprev[i] = Bnew[i];
    }
    if ((r & 63) == 0) Rcpp::checkUserInterrupt();
  }
  return u;
}

// Forward Rosenblatt transform r_t = P(U_t <= u_t | previous min(t - 1, K) values).
// [[Rcpp::export]]
NumericVector dv_rosenblatt_cpp(NumericVector u, IntegerVector fam, IntegerVector rot, NumericVector par,
                                IntegerVector ivt, IntegerVector bfam, IntegerVector brot, NumericVector bpar,
                                double d1, double d2) {
  const int n = u.size();
  const int K = std::min<int>(fam.size(), n - 1);
  NumericVector r = clone(u);
  if (K < 1) return r;
  const std::vector<Level> lv = make_levels(fam, rot, par, ivt, bfam, brot, bpar, d1, d2);
  std::vector<double> A(u.begin(), u.end() - 1), B(u.begin() + 1, u.end()), An(n), Bn(n);
  int m = n - 1;
  for (int j = 1; j <= K; j++) {
    const Level &L = lv[j - 1];
    r[j] = lv_h(L, A[0], B[0], 1);
    if (j < K) {
      for (int i = 0; i < m - 1; i++) { An[i] = lv_h(L, A[i], B[i], 2); Bn[i] = lv_h(L, A[i + 1], B[i + 1], 1); }
      for (int i = 0; i < m - 1; i++) { A[i] = An[i]; B[i] = Bn[i]; }
      m -= 1;
    } else if (m > 1) {
      for (int i = 0; i < m - 1; i++) r[K + 1 + i] = lv_h(L, A[i + 1], B[i + 1], 1);
    }
  }
  return r;
}
