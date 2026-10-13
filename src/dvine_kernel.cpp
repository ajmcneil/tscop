// C++ kernel for the d-vine likelihood of tied sdvinecop models: plain pair
// copulas (independence, Gauss, Clayton, Gumbel, Joe, Frank, all four
// rotations) and inverse-v-transformed (IVT) pair copulas with linear
// v-transforms and these bases, evaluated by the same level recursion as the R
// reference implementation dv_loglik_r (R/dvine-engine.R). The recursion is a
// template over the scalar type: double gives the log-likelihood, and a
// forward-mode dual number with P derivative directions gives the
// log-likelihood and its gradient with respect to the pair copula parameters
// and the v-transform fulcrums (the R side chains these to the free
// parameters). Portable scalar code: no OpenMP, no compiler-specific flags.
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace {

const double EPS = 1e-12;       // clamp of the R reference (dv_clamp)
const double EPS_IN = 1e-10;    // rvinecopulib clamps the arguments of its copula functions at this value
const int MAXP = 12;            // largest number of derivative directions

enum Family { INDEP = 0, GAUSS = 1, CLAYTON = 2, GUMBEL = 3, JOE = 4, FRANK = 5 };

// ---- forward-mode dual numbers with P derivative directions ----------------

template <int P>
struct Dual {
  double v;
  double d[P];
  Dual() : v(0.0) { for (int k = 0; k < P; k++) d[k] = 0.0; }
  Dual(double x) : v(x) { for (int k = 0; k < P; k++) d[k] = 0.0; }
};

template <int P> inline Dual<P> operator+(const Dual<P> &a, const Dual<P> &b) {
  Dual<P> r; r.v = a.v + b.v; for (int k = 0; k < P; k++) r.d[k] = a.d[k] + b.d[k]; return r;
}
template <int P> inline Dual<P> operator+(const Dual<P> &a, double b) { Dual<P> r = a; r.v += b; return r; }
template <int P> inline Dual<P> operator+(double a, const Dual<P> &b) { return b + a; }
template <int P> inline Dual<P> operator-(const Dual<P> &a, const Dual<P> &b) {
  Dual<P> r; r.v = a.v - b.v; for (int k = 0; k < P; k++) r.d[k] = a.d[k] - b.d[k]; return r;
}
template <int P> inline Dual<P> operator-(const Dual<P> &a, double b) { Dual<P> r = a; r.v -= b; return r; }
template <int P> inline Dual<P> operator-(double a, const Dual<P> &b) {
  Dual<P> r; r.v = a - b.v; for (int k = 0; k < P; k++) r.d[k] = -b.d[k]; return r;
}
template <int P> inline Dual<P> operator-(const Dual<P> &a) {
  Dual<P> r; r.v = -a.v; for (int k = 0; k < P; k++) r.d[k] = -a.d[k]; return r;
}
template <int P> inline Dual<P> operator*(const Dual<P> &a, const Dual<P> &b) {
  Dual<P> r; r.v = a.v * b.v; for (int k = 0; k < P; k++) r.d[k] = a.d[k] * b.v + a.v * b.d[k]; return r;
}
template <int P> inline Dual<P> operator*(const Dual<P> &a, double b) {
  Dual<P> r; r.v = a.v * b; for (int k = 0; k < P; k++) r.d[k] = a.d[k] * b; return r;
}
template <int P> inline Dual<P> operator*(double a, const Dual<P> &b) { return b * a; }
template <int P> inline Dual<P> operator/(const Dual<P> &a, const Dual<P> &b) {
  Dual<P> r; const double ib = 1.0 / b.v; r.v = a.v * ib;
  for (int k = 0; k < P; k++) r.d[k] = (a.d[k] - r.v * b.d[k]) * ib;
  return r;
}
template <int P> inline Dual<P> operator/(const Dual<P> &a, double b) { return a * (1.0 / b); }
template <int P> inline Dual<P> operator/(double a, const Dual<P> &b) {
  Dual<P> r; r.v = a / b.v; const double f = -r.v / b.v;
  for (int k = 0; k < P; k++) r.d[k] = f * b.d[k];
  return r;
}
template <int P> inline Dual<P> &operator+=(Dual<P> &a, const Dual<P> &b) { a = a + b; return a; }

// derivative of a unary function: value fv and derivative multiplier fd at a.v
template <int P> inline Dual<P> chain(const Dual<P> &a, double fv, double fd) {
  Dual<P> r; r.v = fv; for (int k = 0; k < P; k++) r.d[k] = fd * a.d[k]; return r;
}
template <int P> inline Dual<P> log(const Dual<P> &a) { return chain(a, std::log(a.v), 1.0 / a.v); }
template <int P> inline Dual<P> exp(const Dual<P> &a) { const double e = std::exp(a.v); return chain(a, e, e); }
template <int P> inline Dual<P> expm1(const Dual<P> &a) { return chain(a, std::expm1(a.v), std::exp(a.v)); }
template <int P> inline Dual<P> log1p(const Dual<P> &a) { return chain(a, std::log1p(a.v), 1.0 / (1.0 + a.v)); }
template <int P> inline Dual<P> sqrt(const Dual<P> &a) { const double s = std::sqrt(a.v); return chain(a, s, 0.5 / s); }
template <int P> inline Dual<P> fabs(const Dual<P> &a) { return a.v < 0.0 ? -a : a; }

inline double val(double x) { return x; }
template <int P> inline double val(const Dual<P> &x) { return x.v; }

// standard normal distribution function and quantile function
inline double pnorm_(double x) { return R::pnorm(x, 0, 1, 1, 0); }
template <int P> inline Dual<P> pnorm_(const Dual<P> &x) {
  return chain(x, R::pnorm(x.v, 0, 1, 1, 0), R::dnorm(x.v, 0, 1, 0));
}
inline double qnorm_(double p) { return R::qnorm(p, 0, 1, 1, 0); }
template <int P> inline Dual<P> qnorm_(const Dual<P> &p) {
  const double q = R::qnorm(p.v, 0, 1, 1, 0);
  return chain(p, q, 1.0 / R::dnorm(q, 0, 1, 0));
}

// clamps: zero derivative where the clamp is active
template <class T> inline T clampu(const T &x) {
  const double v = val(x);
  return v < EPS ? T(EPS) : (v > 1.0 - EPS ? T(1.0 - EPS) : x);
}
template <class T> inline T clampin(const T &x) {
  const double v = val(x);
  return v < EPS_IN ? T(EPS_IN) : (v > 1.0 - EPS_IN ? T(1.0 - EPS_IN) : x);
}

// ---- rotation 0 pair copulas: log density and h-functions -------------------
// h1 = P(B <= b | A = a), h2 = P(A <= a | B = b), for a, b in (0, 1).

template <class T>
inline void clayton0(const T &a, const T &b, const T &th, T &lc, T &h1, T &h2) {
  using std::log; using std::exp; using std::log1p; using std::expm1;
  const T la = log(a), lb = log(b);
  // ta = a^-th, tb = b^-th; s - 1 = (ta - 1) + (tb - 1) is formed from expm1 so that
  // log(s) = log1p(...) stays accurate for small th (nearly independent lags)
  const T ea = expm1(-th * la), eb = expm1(-th * lb);
  const T ta = ea + 1.0, tb = eb + 1.0;
  const T ls = log1p(ea + eb);
  lc = log1p(th) - (th + 1.0) * (la + lb) - (2.0 + 1.0 / th) * ls;
  const T f = exp(-(1.0 + 1.0 / th) * ls);
  h1 = (ta / a) * f;
  h2 = (tb / b) * f;
}

template <class T>
inline void gumbel0(const T &a, const T &b, const T &th, T &lc, T &h1, T &h2) {
  using std::log; using std::exp; using std::log1p;
  const T x = -log(a), y = -log(b);
  const T lx = log(x), ly = log(y);
  const T tx = exp(th * lx), ty = exp(th * ly);
  const T t = tx + ty;
  const T lt = log(t);
  const T Tt = exp(lt / th);
  lc = -Tt - log(a) - log(b) + (th - 1.0) * (lx + ly) + (2.0 / th - 2.0) * lt + log1p((th - 1.0) / Tt);
  const T g = exp(-Tt + (1.0 / th - 1.0) * lt);
  h1 = g * exp((th - 1.0) * lx) / a;
  h2 = g * exp((th - 1.0) * ly) / b;
}

template <class T>
inline void joe0(const T &a, const T &b, const T &th, T &lc, T &h1, T &h2) {
  using std::log; using std::exp;
  const T ub = 1.0 - a, vb = 1.0 - b;
  const T la = log(ub), lb = log(vb);
  const T A = exp(th * la), B = exp(th * lb);
  const T s = A + B - A * B;
  const T ls = log(s);
  lc = (1.0 / th - 2.0) * ls + (th - 1.0) * (la + lb) + log(th - 1.0 + s);
  const T f = exp((1.0 / th - 1.0) * ls);
  h1 = (A / ub) * (1.0 - B) * f;
  h2 = (B / vb) * (1.0 - A) * f;
}

template <class T>
inline void frank0(const T &a, const T &b, const T &th, T &lc, T &h1, T &h2) {
  using std::log; using std::expm1;
  const T Pu = expm1(-th * a), Pv = expm1(-th * b), E = expm1(-th);
  const T D = E + Pu * Pv;
  const T eu = Pu + 1.0, ev = Pv + 1.0;
  lc = log(th * (-E) / (D * D)) - th * (a + b);
  h1 = eu * Pv / D;
  h2 = ev * Pu / D;
}

// Pair copula of the given family and rotation (rvinecopulib convention) at
// (a, b); h-values only if needh.
template <class T>
inline void pair_eval(int fam, int rot, const T &th, const T &a, const T &b, bool needh, bool clampout,
                      T &lc, T &h1, T &h2) {
  if (fam == INDEP) { lc = T(0.0); h1 = b; h2 = a; return; }
  const bool fa = (rot == 90 || rot == 180), fb = (rot == 270 || rot == 180);
  const T aa = clampin(fa ? 1.0 - a : a), bb = clampin(fb ? 1.0 - b : b);
  T g1, g2;
  switch (fam) {
    case CLAYTON: clayton0(aa, bb, th, lc, g1, g2); break;
    case GUMBEL: gumbel0(aa, bb, th, lc, g1, g2); break;
    case JOE: joe0(aa, bb, th, lc, g1, g2); break;
    case FRANK: frank0(aa, bb, th, lc, g1, g2); break;
    default: stop("unsupported family in the kernel");
  }
  if (needh) {
    h1 = fb ? 1.0 - g1 : g1;
    h2 = fa ? 1.0 - g2 : g2;
    if (clampout) { h1 = clampu(h1); h2 = clampu(h2); }
  }
}

template <class T> inline T vlin(const T &x, const T &d) {
  return val(x) > val(d) ? (x - d) / (1.0 - d) : (d - x) / d;
}
template <class T> inline T alin(const T &x, const T &d) { return val(x) > val(d) ? d - 1.0 : d; }

template <class T> inline T gauss_ll(const T &x, const T &y, const T &rho) {
  using std::log1p;
  const T r2 = rho * rho;
  return -0.5 * log1p(-r2) - (r2 * (x * x + y * y) - 2.0 * rho * x * y) / (2.0 * (1.0 - r2));
}

// The level recursion. Level j (0-based) has family code fam[j], rotation
// rot[j], parameter par[j]; an IVT level has ivt[j] = 1 and its base copula in
// bfam[j], brot[j], bpar[j]. The v-transform fulcrums are d1, d2.
template <class T>
T loglik_impl(const NumericVector &u, const IntegerVector &fam, const IntegerVector &rot,
              const std::vector<T> &par, const IntegerVector &ivt, const IntegerVector &bfam,
              const IntegerVector &brot, const std::vector<T> &bpar, const T &d1, const T &d2) {
  using std::sqrt; using std::log1p;
  const int n = u.size();
  const int K = std::min<int>(fam.size(), n - 1);
  if (K < 1) return T(0.0);
  int m = n - 1;
  std::vector<T> A(u.begin(), u.end() - 1), B(u.begin() + 1, u.end());
  std::vector<T> sA, sB, h1(m), h2(m);
  bool zs = false, bs = false, zb = false;
  T ll(0.0);
  for (int j = 0; j < K; j++) {
    const bool last = (j == K - 1);
    if (bs && !ivt[j]) {
      if (zb) {
        for (int i = 0; i < m; i++) { A[i] = pnorm_(A[i]); B[i] = pnorm_(B[i]); }
        zb = false;
      }
      for (int i = 0; i < m; i++) { A[i] = clampu(d1 - sA[i] * A[i]); B[i] = clampu(d2 - sB[i] * B[i]); }
      bs = false;
    }
    const bool plain_gauss = !ivt[j] && fam[j] == GAUSS;
    if (plain_gauss) {
      if (!zs) {
        for (int i = 0; i < m; i++) { A[i] = qnorm_(A[i]); B[i] = qnorm_(B[i]); }
        zs = true;
      }
      const T rho = par[j];
      const T s = sqrt(1.0 - rho * rho);
      for (int i = 0; i < m; i++) ll = ll + gauss_ll(A[i], B[i], rho);
      if (!last) {
        for (int i = 0; i < m - 1; i++) {
          h1[i] = (A[i] - rho * B[i]) / s;
          h2[i] = (B[i + 1] - rho * A[i + 1]) / s;
        }
        for (int i = 0; i < m - 1; i++) { A[i] = h1[i]; B[i] = h2[i]; }
        m -= 1;
      }
    } else if (ivt[j]) {
      if (!bs) {
        if (zs) {
          for (int i = 0; i < m; i++) { A[i] = pnorm_(A[i]); B[i] = pnorm_(B[i]); }
          zs = false;
        }
        sA.assign(m, T(0.0)); sB.assign(m, T(0.0));
        for (int i = 0; i < m; i++) {
          sA[i] = alin(A[i], d1); sB[i] = alin(B[i], d2);
          A[i] = vlin(A[i], d1); B[i] = vlin(B[i], d2);
        }
        bs = true; zb = false;
      }
      const bool gb = (bfam[j] == GAUSS);
      if (gb && !zb) {
        for (int i = 0; i < m; i++) { A[i] = qnorm_(clampu(A[i])); B[i] = qnorm_(clampu(B[i])); }
        zb = true;
      } else if (!gb && zb) {
        for (int i = 0; i < m; i++) { A[i] = pnorm_(A[i]); B[i] = pnorm_(B[i]); }
        zb = false;
      }
      if (gb) {
        const T rho = bpar[j];
        const T s = sqrt(1.0 - rho * rho);
        for (int i = 0; i < m; i++) ll = ll + gauss_ll(A[i], B[i], rho);
        if (!last) {
          for (int i = 0; i < m - 1; i++) {
            h1[i] = (A[i] - rho * B[i]) / s;
            h2[i] = (B[i + 1] - rho * A[i + 1]) / s;
          }
          for (int i = 0; i < m - 1; i++) { A[i] = h1[i]; B[i] = h2[i]; }
        }
      } else {
        for (int i = 0; i < m; i++) {
          T lc;
          pair_eval(bfam[j], brot[j], bpar[j], A[i], B[i], !last, true, lc, h1[i], h2[i]);
          ll = ll + lc;
        }
        if (!last) {
          for (int i = 0; i < m - 1; i++) { A[i] = h2[i]; B[i] = h1[i + 1]; }
        }
      }
      if (!last) {
        for (int i = 0; i < m - 1; i++) sB[i] = sB[i + 1];
        m -= 1;
      }
    } else {
      if (zs) {
        for (int i = 0; i < m; i++) { A[i] = pnorm_(A[i]); B[i] = pnorm_(B[i]); }
        zs = false;
      }
      for (int i = 0; i < m; i++) {
        T lc;
        pair_eval(fam[j], rot[j], par[j], A[i], B[i], !last, false, lc, h1[i], h2[i]);
        ll = ll + lc;
      }
      if (!last) {
        for (int i = 0; i < m - 1; i++) { A[i] = h2[i]; B[i] = h1[i + 1]; }
        m -= 1;
      }
    }
  }
  return ll;
}

// Gradient with P derivative directions: dpar and dbpar are K x P matrices of
// the derivatives of the level parameters, dd1 and dd2 the derivatives of the
// fulcrums, all with respect to the P directions.
template <int P>
List grad_impl(const NumericVector &u, const IntegerVector &fam, const IntegerVector &rot,
               const NumericVector &par, const NumericMatrix &dpar, const IntegerVector &ivt,
               const IntegerVector &bfam, const IntegerVector &brot, const NumericVector &bpar,
               const NumericMatrix &dbpar, double d1, double d2, const NumericVector &dd1,
               const NumericVector &dd2) {
  const int K = fam.size();
  std::vector<Dual<P> > p(K), bp(K);
  for (int j = 0; j < K; j++) {
    p[j].v = par[j]; bp[j].v = bpar[j];
    for (int k = 0; k < P; k++) { p[j].d[k] = dpar(j, k); bp[j].d[k] = dbpar(j, k); }
  }
  Dual<P> D1(d1), D2(d2);
  for (int k = 0; k < P; k++) { D1.d[k] = dd1[k]; D2.d[k] = dd2[k]; }
  Dual<P> ll = loglik_impl<Dual<P> >(u, fam, rot, p, ivt, bfam, brot, bp, D1, D2);
  NumericVector g(P);
  for (int k = 0; k < P; k++) g[k] = ll.d[k];
  if (!std::isfinite(ll.v)) return List::create(_["value"] = R_NegInf, _["gradient"] = g);
  return List::create(_["value"] = ll.v, _["gradient"] = g);
}

}  // namespace

// Log-likelihood of the d-vine.
// [[Rcpp::export]]
double dv_loglik_cpp(NumericVector u, IntegerVector fam, IntegerVector rot, NumericVector par,
                     IntegerVector ivt, IntegerVector bfam, IntegerVector brot, NumericVector bpar,
                     double d1, double d2) {
  const int K = fam.size();
  std::vector<double> p(par.begin(), par.end()), bp(bpar.begin(), bpar.end());
  (void)K;
  const double ll = loglik_impl<double>(u, fam, rot, p, ivt, bfam, brot, bp, d1, d2);
  return std::isfinite(ll) ? ll : R_NegInf;
}

// Log-likelihood and gradient with respect to P (1 to 12) directions along which the
// level parameters (K x P matrices dpar, dbpar) and the fulcrums (dd1, dd2) vary.
// [[Rcpp::export]]
List dv_loglik_grad_cpp(NumericVector u, IntegerVector fam, IntegerVector rot, NumericVector par,
                        NumericMatrix dpar, IntegerVector ivt, IntegerVector bfam, IntegerVector brot,
                        NumericVector bpar, NumericMatrix dbpar, double d1, double d2,
                        NumericVector dd1, NumericVector dd2) {
  const int P = dpar.ncol();
  switch (P) {
#define CASE(PP) case PP: return grad_impl<PP>(u, fam, rot, par, dpar, ivt, bfam, brot, bpar, dbpar, d1, d2, dd1, dd2);
    CASE(1) CASE(2) CASE(3) CASE(4) CASE(5) CASE(6) CASE(7) CASE(8) CASE(9) CASE(10) CASE(11) CASE(12)
#undef CASE
    default: stop("the kernel gradient supports 1 to 12 directions");
  }
}
