// C++ kernel for the d-vine likelihood of tied sdvinecop models: plain pair
// copulas (independence, Gauss, Clayton, Gumbel, Joe, Frank, all four
// rotations) and inverse-v-transformed (IVT) pair copulas with linear
// v-transforms and these bases, evaluated by the same level recursion as the R
// reference implementation dv_loglik_r (R/dvine-engine.R). Portable scalar
// code: no OpenMP, no compiler-specific flags.
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace {

const double EPS = 1e-12;       // clamp of the R reference (dv_clamp)
const double EPS_IN = 1e-10;    // rvinecopulib clamps the arguments of its copula functions at this value

enum Family { INDEP = 0, GAUSS = 1, CLAYTON = 2, GUMBEL = 3, JOE = 4, FRANK = 5 };

inline double clampu(double x) { return x < EPS ? EPS : (x > 1.0 - EPS ? 1.0 - EPS : x); }
inline double clampin(double x) { return x < EPS_IN ? EPS_IN : (x > 1.0 - EPS_IN ? 1.0 - EPS_IN : x); }

// Rotation 0 pair copulas: log density and the h-functions
// h1 = P(B <= b | A = a), h2 = P(A <= a | B = b), for a, b in (0, 1).
inline void clayton0(double a, double b, double th, double &lc, double &h1, double &h2) {
  const double la = std::log(a), lb = std::log(b);
  const double ta = std::exp(-th * la), tb = std::exp(-th * lb);
  const double s = ta + tb - 1.0;
  const double ls = std::log(s);
  lc = std::log1p(th) - (th + 1.0) * (la + lb) - (2.0 + 1.0 / th) * ls;
  const double f = std::exp(-(1.0 + 1.0 / th) * ls);
  h1 = (ta / a) * f;
  h2 = (tb / b) * f;
}

inline void gumbel0(double a, double b, double th, double &lc, double &h1, double &h2) {
  const double x = -std::log(a), y = -std::log(b);
  const double lx = std::log(x), ly = std::log(y);
  const double tx = std::exp(th * lx), ty = std::exp(th * ly);
  const double t = tx + ty;
  const double lt = std::log(t);
  const double T = std::exp(lt / th);
  // log C = -T
  lc = -T - std::log(a) - std::log(b) + (th - 1.0) * (lx + ly) + (2.0 / th - 2.0) * lt +
       std::log1p((th - 1.0) / T);
  const double g = std::exp(-T + (1.0 / th - 1.0) * lt);
  h1 = g * std::exp((th - 1.0) * lx) / a;
  h2 = g * std::exp((th - 1.0) * ly) / b;
}

inline void joe0(double a, double b, double th, double &lc, double &h1, double &h2) {
  const double ub = 1.0 - a, vb = 1.0 - b;
  const double la = std::log(ub), lb = std::log(vb);
  const double A = std::exp(th * la), B = std::exp(th * lb);
  const double s = A + B - A * B;
  const double ls = std::log(s);
  lc = (1.0 / th - 2.0) * ls + (th - 1.0) * (la + lb) + std::log(th - 1.0 + s);
  const double f = std::exp((1.0 / th - 1.0) * ls);
  h1 = (A / ub) * (1.0 - B) * f;
  h2 = (B / vb) * (1.0 - A) * f;
}

inline void frank0(double a, double b, double th, double &lc, double &h1, double &h2) {
  const double Pu = std::expm1(-th * a), Pv = std::expm1(-th * b), E = std::expm1(-th);
  const double D = E + Pu * Pv;
  const double eu = Pu + 1.0, ev = Pv + 1.0;
  lc = std::log(th * (-E) / (D * D)) - th * (a + b);
  h1 = eu * Pv / D;
  h2 = ev * Pu / D;
}

// Pair copula of the given family and rotation (rvinecopulib convention) at
// (a, b); h-values only if needh.
inline void pair_eval(int fam, int rot, double th, double a, double b, bool needh, bool clampout,
                      double &lc, double &h1, double &h2) {
  if (fam == INDEP) { lc = 0.0; h1 = b; h2 = a; return; }
  const bool fa = (rot == 90 || rot == 180), fb = (rot == 270 || rot == 180);
  double aa = clampin(fa ? 1.0 - a : a), bb = clampin(fb ? 1.0 - b : b);
  double g1, g2;
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

inline double dv_vlin(double x, double d) { return std::fabs(x - d) / (x > d ? 1.0 - d : d); }
inline double dv_alin(double x, double d) { return x > d ? d - 1.0 : d; }

inline double gauss_ll(double x, double y, double rho) {
  const double r2 = rho * rho;
  return -0.5 * std::log1p(-r2) - (r2 * (x * x + y * y) - 2.0 * rho * x * y) / (2.0 * (1.0 - r2));
}

}  // namespace

// Log-likelihood of the d-vine. Level j (0-based) has family code fam[j],
// rotation rot[j], parameter par[j]; an IVT level has ivt[j] = 1 and its base
// copula in bfam[j], brot[j], bpar[j]. The v-transform fulcrums are d1, d2.
// [[Rcpp::export]]
double dv_loglik_cpp(NumericVector u, IntegerVector fam, IntegerVector rot, NumericVector par,
                     IntegerVector ivt, IntegerVector bfam, IntegerVector brot, NumericVector bpar,
                     double d1, double d2) {
  const int n = u.size();
  const int K = std::min<int>(fam.size(), n - 1);
  if (K < 1) return 0.0;
  int m = n - 1;
  std::vector<double> A(u.begin(), u.end() - 1), B(u.begin() + 1, u.end());
  std::vector<double> sA, sB, h1(m), h2(m);
  bool zs = false, bs = false, zb = false;
  double ll = 0.0;
  for (int j = 0; j < K; j++) {
    const bool last = (j == K - 1);
    if (bs && !ivt[j]) {
      if (zb) {
        for (int i = 0; i < m; i++) { A[i] = R::pnorm(A[i], 0, 1, 1, 0); B[i] = R::pnorm(B[i], 0, 1, 1, 0); }
        zb = false;
      }
      for (int i = 0; i < m; i++) { A[i] = clampu(d1 - sA[i] * A[i]); B[i] = clampu(d2 - sB[i] * B[i]); }
      bs = false;
    }
    const bool plain_gauss = !ivt[j] && fam[j] == GAUSS;
    if (plain_gauss) {
      if (!zs) {
        for (int i = 0; i < m; i++) { A[i] = R::qnorm(A[i], 0, 1, 1, 0); B[i] = R::qnorm(B[i], 0, 1, 1, 0); }
        zs = true;
      }
      const double rho = par[j], s = std::sqrt(1.0 - rho * rho);
      for (int i = 0; i < m; i++) ll += gauss_ll(A[i], B[i], rho);
      if (!last) {
        for (int i = 0; i < m - 1; i++) {
          h1[i] = (A[i] - rho * B[i]) / s;          // new A
          h2[i] = (B[i + 1] - rho * A[i + 1]) / s;  // new B
        }
        for (int i = 0; i < m - 1; i++) { A[i] = h1[i]; B[i] = h2[i]; }
        m -= 1;
      }
    } else if (ivt[j]) {
      if (!bs) {
        if (zs) {
          for (int i = 0; i < m; i++) { A[i] = R::pnorm(A[i], 0, 1, 1, 0); B[i] = R::pnorm(B[i], 0, 1, 1, 0); }
          zs = false;
        }
        sA.assign(m, 0.0); sB.assign(m, 0.0);
        for (int i = 0; i < m; i++) {
          sA[i] = dv_alin(A[i], d1); sB[i] = dv_alin(B[i], d2);
          A[i] = dv_vlin(A[i], d1); B[i] = dv_vlin(B[i], d2);
        }
        bs = true; zb = false;
      }
      const bool gb = (bfam[j] == GAUSS);
      if (gb && !zb) {
        for (int i = 0; i < m; i++) { A[i] = R::qnorm(clampu(A[i]), 0, 1, 1, 0); B[i] = R::qnorm(clampu(B[i]), 0, 1, 1, 0); }
        zb = true;
      } else if (!gb && zb) {
        for (int i = 0; i < m; i++) { A[i] = R::pnorm(A[i], 0, 1, 1, 0); B[i] = R::pnorm(B[i], 0, 1, 1, 0); }
        zb = false;
      }
      if (gb) {
        const double rho = bpar[j], s = std::sqrt(1.0 - rho * rho);
        for (int i = 0; i < m; i++) ll += gauss_ll(A[i], B[i], rho);
        if (!last) {
          for (int i = 0; i < m - 1; i++) {
            h1[i] = (A[i] - rho * B[i]) / s;
            h2[i] = (B[i + 1] - rho * A[i + 1]) / s;
          }
          for (int i = 0; i < m - 1; i++) { A[i] = h1[i]; B[i] = h2[i]; }
        }
      } else {
        for (int i = 0; i < m; i++) {
          double lc;
          pair_eval(bfam[j], brot[j], bpar[j], A[i], B[i], !last, true, lc, h1[i], h2[i]);
          ll += lc;
        }
        if (!last) {
          for (int i = 0; i < m - 1; i++) { A[i] = h2[i]; B[i] = h1[i + 1]; }
        }
      }
      if (!last) {
        for (int i = 0; i < m - 1; i++) { sA[i] = sA[i]; sB[i] = sB[i + 1]; }
        m -= 1;
      }
    } else {
      if (zs) {
        for (int i = 0; i < m; i++) { A[i] = R::pnorm(A[i], 0, 1, 1, 0); B[i] = R::pnorm(B[i], 0, 1, 1, 0); }
        zs = false;
      }
      for (int i = 0; i < m; i++) {
        double lc;
        pair_eval(fam[j], rot[j], par[j], A[i], B[i], !last, false, lc, h1[i], h2[i]);
        ll += lc;
      }
      if (!last) {
        for (int i = 0; i < m - 1; i++) { A[i] = h2[i]; B[i] = h1[i + 1]; }
        m -= 1;
      }
    }
  }
  return std::isfinite(ll) ? ll : R_NegInf;
}
