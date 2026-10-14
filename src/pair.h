// Shared building blocks of the d-vine kernels (likelihood, simulation,
// Rosenblatt transform): dual numbers, clamps and the pair copula functions of
// the Clayton, Gumbel, Joe and Frank copulas (rotation 0) with rotations. The
// functions are templates over the scalar type (double or a dual number).
#pragma once
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;
namespace {

const double EPS = 1e-12;       // clamp of the R reference (dv_clamp)
const double EPS_IN = 1e-10;    // rvinecopulib clamps the arguments of its copula functions at this value
const int MAXP = 12;            // largest number of derivative directions

enum Family { INDEP = 0, GAUSS = 1, CLAYTON = 2, GUMBEL = 3, JOE = 4, FRANK = 5, BB1 = 6 };

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

// BB1 (Joe's two-parameter Clayton-Gumbel family; th > 0, dl >= 1):
// C(u, v) = (1 + t)^(-1/th), t = (x^dl + y^dl)^(1/dl), x = u^-th - 1, y = v^-th - 1.
template <class T>
inline void bb10(const T &a, const T &b, const T &th, const T &dl, T &lc, T &h1, T &h2) {
  using std::log; using std::exp; using std::log1p; using std::expm1;
  const T la = log(a), lb = log(b);
  const T x = expm1(-th * la), y = expm1(-th * lb);
  const T lx = log(x), ly = log(y);
  const T s = exp(dl * lx) + exp(dl * ly);
  const T ls = log(s);
  const T t = exp(ls / dl);
  const T l1t = log1p(t);
  lc = (dl - 1.0) * (lx + ly) - (th + 1.0) * (la + lb) - (1.0 / th + 2.0) * l1t + (1.0 / dl - 2.0) * ls +
       log((1.0 + th) * t + th * (dl - 1.0) * (1.0 + t));
  const T AB = exp(-(1.0 / th + 1.0) * l1t + (1.0 / dl - 1.0) * ls);
  h1 = AB * exp((dl - 1.0) * lx - (th + 1.0) * la);
  h2 = AB * exp((dl - 1.0) * ly - (th + 1.0) * lb);
}

// Pair copula of the given family and rotation (rvinecopulib convention) at
// (a, b); h-values only if needh.
template <class T>
inline void pair_eval(int fam, int rot, const T &th, const T &th2, const T &a, const T &b, bool needh, bool clampout,
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
    case BB1: bb10(aa, bb, th, th2, lc, g1, g2); break;
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

}  // namespace
