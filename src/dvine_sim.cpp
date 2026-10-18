// C++ kernels for the d-vine simulation (exact stationary simulation by
// inverting the forward chain) and the forward Rosenblatt transform of tied
// sdvinecop models, with the same pair copulas as the likelihood kernel
// (independence, Gauss, Clayton, Gumbel, Joe, Frank with rotations, and IVT
// pairs with these bases, and nonlinear IVT levels with the Clayton, Gumbel, Joe and Frank bases, see nl.h). They follow dv_simulate_r and dv_rosenblatt_r
// (R/dvine-engine.R) step by step.
#include "pair.h"
#include "ast.h"
#include "nl.h"
#include <memory>

namespace {

struct Level {
  int fam, rot, bfam, brot;
  bool ivt;
  bool nl = false;               // nonlinear IVT level (base in bfam, brot, bth; shapes k1, k2 of the v2p transforms)
  double th, th2, bth, d1, d2;   // th2: the second parameter of a plain BB1 level
  double k1 = 1.0, k2 = 1.0;
  std::shared_ptr<AstZTab> tab;  // t-quantile table of an ast base of a nonlinear level
  double c0 = 0.0;               // its log density constant
};

inline double pn(double x) { return R::pnorm(x, 0, 1, 1, 0); }
inline double qn(double p) { return R::qnorm(p, 0, 1, 1, 0); }

// ---- absolute spherical t copula, exact (th = nu; rotation 0) ----
// h(given u | x): 1 - 2 pt_upper(b sqrt((nu + 1) / (nu + a^2)), nu + 1) with a, b the quantiles of the
// given and the other variable
inline double ast_h_ab_d(double a, double b, double nu) {
  if (std::isinf(b)) return 1.0;
  if (std::isinf(a)) return 0.0;
  return 1.0 - 2.0 * R::pt(b * std::sqrt((nu + 1.0) / (nu + a * a)), nu + 1.0, 0, 0);
}
inline double ast_dens_d(double a, double b, double nu) {
  const double x = ast_q(a, nu), y = ast_q(b, nu);
  const double c0 = R::lgammafn((nu + 2.0) / 2.0) - 2.0 * R::lgammafn((nu + 1.0) / 2.0) + R::lgammafn(nu / 2.0);
  const double x2 = x * x / nu, y2 = y * y / nu;
  return std::exp(c0 - (nu + 2.0) / 2.0 * std::log1p(x2 + y2) + (nu + 1.0) / 2.0 * (std::log1p(x2) + std::log1p(y2)));
}
inline double ast_hinv_d(double a, double p, double nu) {
  const double qa = ast_q(a, nu);
  const double b = std::sqrt((nu + qa * qa) / (nu + 1.0)) * R::qt((1.0 - p) / 2.0, nu + 1.0, 0, 0);
  return 1.0 - 2.0 * R::pt(b, nu, 0, 0);
}

// ---- t copula (th = rho, th2 = nu) ----
inline double tc_h_d(double a, double b, double rho, double nu, int cv) {
  const double x = R::qt(clampin(a), nu, 1, 0), y = R::qt(clampin(b), nu, 1, 0);
  const double gv = cv == 1 ? x : y, ot = cv == 1 ? y : x;     // given and other
  return R::pt((ot - rho * gv) / std::sqrt((nu + gv * gv) * (1.0 - rho * rho) / (nu + 1.0)), nu + 1.0, 1, 0);
}
inline double tc_dens_d(double a, double b, double rho, double nu) {
  const double x = R::qt(clampin(a), nu, 1, 0), y = R::qt(clampin(b), nu, 1, 0);
  const double c0 = R::lgammafn((nu + 2.0) / 2.0) + R::lgammafn(nu / 2.0) - 2.0 * R::lgammafn((nu + 1.0) / 2.0) - 0.5 * std::log1p(-rho * rho);
  const double Q = (x * x - 2.0 * rho * x * y + y * y) / (1.0 - rho * rho);
  return std::exp(c0 - (nu + 2.0) / 2.0 * std::log1p(Q / nu) + (nu + 1.0) / 2.0 * (std::log1p(x * x / nu) + std::log1p(y * y / nu)));
}
inline double tc_hinv_d(double a, double p, double rho, double nu) {
  const double x = R::qt(clampin(a), nu, 1, 0);
  const double y = rho * x + std::sqrt((nu + x * x) * (1.0 - rho * rho) / (nu + 1.0)) * R::qt(p, nu + 1.0, 1, 0);
  return R::pt(y, nu, 1, 0);
}

// h-function of a (non-IVT) pair copula: cv = 1: P(B <= b | A = a), cv = 2: P(A <= a | B = b)
inline double h_plain(int fam, int rot, double th, double th2, double a, double b, int cv) {
  if (fam == TCOP) return tc_h_d(a, b, th, th2, cv);
  if (fam == AST) return cv == 1 ? ast_h_ab_d(ast_q(a, th), ast_q(b, th), th) : ast_h_ab_d(ast_q(b, th), ast_q(a, th), th);
  if (fam == GAUSS) {
    const double s = std::sqrt(1.0 - th * th);
    return cv == 1 ? pn((qn(clampin(b)) - th * qn(clampin(a))) / s)
                   : pn((qn(clampin(a)) - th * qn(clampin(b))) / s);
  }
  double lc, h1, h2;
  pair_eval(fam, rot, th, th2, a, b, true, false, lc, h1, h2);
  return cv == 1 ? h1 : h2;
}

// inverse of h(a, .) (cv = 1) of a rotation 0 copula, by closed form or safeguarded Newton
inline double hinv0(int fam, double th, double th2, double a, double p) {
  a = clampin(a);
  p = clampin(p);
  switch (fam) {
    case AST: return ast_hinv_d(a, p, th);
    case TCOP: return tc_hinv_d(a, p, th, th2);
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
    if (fam == GUMBEL) gumbel0(a, b, th, lc, h1, h2); else if (fam == JOE) joe0(a, b, th, lc, h1, h2); else bb10(a, b, th, th2, lc, h1, h2);
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

inline double hinv_plain(int fam, int rot, double th, double th2, double a, double p) {
  switch (rot) {
    case 180: return 1.0 - hinv0(fam, th, th2, 1.0 - a, 1.0 - p);
    case 90: return hinv0(fam, th, th2, 1.0 - a, p);
    case 270: return 1.0 - hinv0(fam, th, th2, a, 1.0 - p);
    default: return hinv0(fam, th, th2, a, p);
  }
}

inline double vl(double x, double d) { return vlin(x, d); }
inline double al(double x, double d) { return alin(x, d); }

// run f(ld, hb) with the log density and h-functions of the base copula of a nonlinear level
template <class F> inline double nl_run(const Level &L, F &&f) {
  if (L.bfam == AST) {
    const AstPolicy<double> ap{L.tab.get(), L.bth, L.c0};
    return f([&](double a, double b) { return ap.ld(a, b); }, [&](double a, double b, double &h1, double &h2) { ap.hb(a, b, h1, h2); });
  }
  return f([&](double a, double b) { return nl_base_ld(L.bfam, L.brot, L.bth, L.bth, a, b); },
           [&](double a, double b, double &h1, double &h2) { double lc; pair_eval(L.bfam, L.brot, L.bth, L.bth, a, b, true, true, lc, h1, h2); });
}

inline double nl_eps(const Level &L) { return L.bfam == AST ? EPS : EPS_IN; }

inline double lv_h(const Level &L, double a, double b, int cv) {
  if (L.nl) {
    const double eps = nl_eps(L);
    return nl_run(L, [&](const auto &ld, const auto &hb) {
      if (cv == 1) return nl_cdf_g(L.d2, L.k2, b, v2p_T(a, L.d1, L.k1, eps), false, ld, hb, eps);
      return nl_cdf_g(L.d1, L.k1, a, v2p_T(b, L.d2, L.k2, eps), true, ld, hb, eps);
    });
  }
  if (!L.ivt) return clampu(h_plain(L.fam, L.rot, L.th, L.th2, a, b, cv));
  const double x = vl(a, L.d1), y = vl(b, L.d2);
  const double hb = clampu(h_plain(L.bfam, L.brot, L.bth, 0.0, x, y, cv));
  return clampu(cv == 1 ? L.d2 - al(b, L.d2) * hb : L.d1 - al(a, L.d1) * hb);
}

// inverse of h(a, .) (cv = 1) of a nonlinear IVT level: safeguarded Newton iteration on the monotone distribution
// function, with the density c_V(T1(a), T2(b)) as derivative and bisection whenever a step leaves the bracket
inline double nl_hinv(const Level &L, double a, double p) {
  a = clampu(a); p = clampu(p);
  const double eps = nl_eps(L);
  const double ta = v2p_T(a, L.d1, L.k1, eps);
  return nl_run(L, [&](const auto &ld, const auto &hb) {
    double lo = 0.0, hi = 1.0, b = p;
    for (int it = 0; it < 80; it++) {
      const double f = nl_cdf_g(L.d2, L.k2, b, ta, false, ld, hb, eps) - p;
      if (f > 0) hi = b; else lo = b;
      if (std::fabs(f) < 1e-13 || hi - lo < 1e-13) break;
      const double dens = std::exp(ld(ta, v2p_T(b, L.d2, L.k2, eps)));
      double bn = b - f / dens;
      if (!std::isfinite(bn) || bn <= lo || bn >= hi) bn = 0.5 * (lo + hi);
      b = bn;
    }
    return b;
  });
}

inline double lv_hinv(const Level &L, double a, double p) {
  if (L.nl) return clampu(nl_hinv(L, a, p));
  if (!L.ivt) return clampu(hinv_plain(L.fam, L.rot, L.th, L.th2, a, p));
  const double hb = clampu(hinv_plain(L.bfam, L.brot, L.bth, 0.0, vl(a, L.d1), vl(p, L.d2)));
  return clampu(L.d2 - al(p, L.d2) * hb);
}

// copula density of a (non-IVT) pair copula at (a, b)
inline double dens_plain(int fam, int rot, double th, double th2, double a, double b) {
  if (fam == AST) return ast_dens_d(a, b, th);
  if (fam == TCOP) return tc_dens_d(a, b, th, th2);
  if (fam == GAUSS) {
    const double x = qn(a), y = qn(b), r2 = th * th;
    return std::exp(-0.5 * std::log1p(-r2) - (r2 * (x * x + y * y) - 2.0 * th * x * y) / (2.0 * (1.0 - r2)));
  }
  double lc, h1, h2;
  pair_eval(fam, rot, th, th2, a, b, false, false, lc, h1, h2);
  return std::exp(lc);
}

inline double lv_dens(const Level &L, double a, double b) {
  if (L.nl) {
    const double eps = nl_eps(L);
    return nl_run(L, [&](const auto &ld, const auto &) { return std::exp(ld(v2p_T(a, L.d1, L.k1, eps), v2p_T(b, L.d2, L.k2, eps))); });
  }
  if (!L.ivt) return dens_plain(L.fam, L.rot, L.th, L.th2, a, b);
  return dens_plain(L.bfam, L.brot, L.bth, 0.0, vl(a, L.d1), vl(b, L.d2));
}

std::vector<Level> make_levels(const IntegerVector &fam, const IntegerVector &rot, const NumericVector &par,
                               const IntegerVector &ivt, const IntegerVector &bfam, const IntegerVector &brot,
                               const NumericVector &bpar, double d1, double d2, const IntegerVector &nl, double k1,
                               double k2) {
  const int K = fam.size();
  std::vector<Level> lv(K);
  for (int j = 0; j < K; j++) {
    lv[j].nl = (j < nl.size() && nl[j] != 0); lv[j].k1 = k1; lv[j].k2 = k2;
    if (lv[j].nl && bfam[j] == AST) {
      lv[j].tab = std::make_shared<AstZTab>();
      lv[j].tab->build(bpar[j], false);
      const double nu = bpar[j];
      lv[j].c0 = R::lgammafn((nu + 2.0) / 2.0) - 2.0 * R::lgammafn((nu + 1.0) / 2.0) + R::lgammafn(nu / 2.0);
    }
    lv[j].fam = fam[j]; lv[j].rot = rot[j]; lv[j].th = par[j];
    lv[j].ivt = ivt[j] != 0; lv[j].bfam = bfam[j]; lv[j].brot = brot[j]; lv[j].bth = bpar[j]; lv[j].th2 = bpar[j];
    lv[j].d1 = d1; lv[j].d2 = d2;
  }
  return lv;
}

}  // namespace

// Exact stationary simulation from independent U(0, 1) innovations w (n x nrep).
// [[Rcpp::export]]
NumericMatrix dv_simulate_cpp(NumericMatrix w, IntegerVector fam, IntegerVector rot, NumericVector par,
                              IntegerVector ivt, IntegerVector bfam, IntegerVector brot, NumericVector bpar,
                              double d1, double d2, Nullable<IntegerVector> nl = R_NilValue, double k1 = 1.0,
                              double k2 = 1.0) {
  const int n = w.nrow(), nrep = w.ncol(), K = fam.size();
  const IntegerVector nlv = nl.isNotNull() ? IntegerVector(nl) : IntegerVector(0);
  const std::vector<Level> lv = make_levels(fam, rot, par, ivt, bfam, brot, bpar, d1, d2, nlv, k1, k2);
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
                                double d1, double d2, Nullable<IntegerVector> nl = R_NilValue, double k1 = 1.0,
                                double k2 = 1.0) {
  const int n = u.size();
  const int K = std::min<int>(fam.size(), n - 1);
  NumericVector r = clone(u);
  if (K < 1) return r;
  const IntegerVector nlv = nl.isNotNull() ? IntegerVector(nl) : IntegerVector(0);
  const std::vector<Level> lv = make_levels(fam, rot, par, ivt, bfam, brot, bpar, d1, d2, nlv, k1, k2);
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

// One-step predictive distribution of the next value given the data u: type 1
// (distribution function), 2 (quantile function) or 3 (density), evaluated at x.
// The backward variables of the last k = min(K, n) values are computed by the
// level recursion on that window (as dv_backward_tail), then the conditional
// pair copulas are applied in turn.
// [[Rcpp::export]]
NumericVector dv_predictive_cpp(NumericVector u, NumericVector x, int type, IntegerVector fam, IntegerVector rot,
                                NumericVector par, IntegerVector ivt, IntegerVector bfam, IntegerVector brot,
                                NumericVector bpar, double d1, double d2, Nullable<IntegerVector> nl = R_NilValue,
                                double k1 = 1.0, double k2 = 1.0) {
  const int n = u.size(), K = fam.size();
  const int k = std::min(K, n);
  const int nx = x.size();
  NumericVector out(nx);
  if (k < 1) {
    for (int i = 0; i < nx; i++) out[i] = type == 3 ? 1.0 : x[i];
    return out;
  }
  const IntegerVector nlv = nl.isNotNull() ? IntegerVector(nl) : IntegerVector(0);
  const std::vector<Level> lv = make_levels(fam, rot, par, ivt, bfam, brot, bpar, d1, d2, nlv, k1, k2);
  // backward variables B_i(n - i), i = 0, ..., k - 1, of the last k values
  std::vector<double> Bt(k), A(k), B(k), An(k), Bn(k);
  for (int i = 0; i < k; i++) A[i] = u[n - k + i];
  Bt[0] = A[k - 1];
  if (k >= 2) {
    for (int i = 0; i < k - 1; i++) { B[i] = A[i + 1]; }
    int m = k - 1;   // length of the current A (first m window values) and B
    for (int i = 1; i <= k - 1; i++) {
      const Level &L = lv[i - 1];
      Bt[i] = lv_h(L, A[m - 1], B[m - 1], 2);
      if (i < k - 1) {
        for (int j = 0; j < m - 1; j++) { An[j] = lv_h(L, A[j], B[j], 2); Bn[j] = lv_h(L, A[j + 1], B[j + 1], 1); }
        for (int j = 0; j < m - 1; j++) { A[j] = An[j]; B[j] = Bn[j]; }
        m -= 1;
      }
    }
  }
  for (int i = 0; i < nx; i++) {
    double f = x[i];
    if (type == 1) {
      for (int j = 1; j <= k; j++) f = lv_h(lv[j - 1], Bt[j - 1], f, 1);
      out[i] = f;
    } else if (type == 2) {
      for (int j = k; j >= 1; j--) f = lv_hinv(lv[j - 1], Bt[j - 1], f);
      out[i] = f;
    } else {
      double d = 1.0;
      for (int j = 1; j <= k; j++) {
        d *= lv_dens(lv[j - 1], Bt[j - 1], f);
        f = lv_h(lv[j - 1], Bt[j - 1], f, 1);
      }
      out[i] = d;
    }
  }
  return out;
}
