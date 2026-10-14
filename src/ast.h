// Absolute spherical t (ast) copula: level maps for the t-scale recursion of the
// d-vine likelihood. With x = t_nu^-1((1 + u) / 2) (upper-tail quantile of (1 - u) / 2),
// the h-function of a pair is h = 1 - 2 pt_upper(z, nu + 1) with
// z = b sqrt((nu + 1) / (nu + a^2)) (a the quantile of the given variable, b that of the
// other), so the next level's t-scale value is x' = qt_upper(pt_upper(z, nu + 1), nu2)
// =: G(z), with nu2 the degrees of freedom of the next level. G is evaluated from a table of
// piecewise Chebyshev polynomials of ln G in ln z (13 panels of width 2, degree 14, 195 nodes;
// linear below e^-12 and exact above e^14), accurate to the 2e-10 floor of the exact
// composition (3e-7 for nu < 1, where R's qt is itself that accurate) and nine times cheaper than
// the two special-function calls. Only for double (no derivatives with respect to nu).
#pragma once
#include <Rcpp.h>
#include <cmath>
#include <vector>

namespace {

inline double ast_q(double u, double nu) { return R::qt((1.0 - u) / 2.0, nu, 0, 0); }   // upper-tail quantile
inline double ast_Gexact(double z, double nu, double nu2) { return R::qt(R::pt(z, nu + 1.0, 0, 0), nu2, 0, 0); }

struct AstMap {
#ifndef AST_D
#define AST_D 14
#endif
#ifndef AST_W
#define AST_W 2.0
#endif
#ifndef AST_NP
#define AST_NP 13
#endif
  static constexpr int D = AST_D;              // polynomial degree
  static constexpr double W = AST_W;           // panel width in ln z
  static constexpr double SMIN = -12.0;
  static constexpr int NP = AST_NP;
  static constexpr double SMAX = SMIN + NP * W;
  double nu = 0, nu2 = 0, c1 = 0, zlo = 0, zhi = 0, xlo = 0, xhi = 0;   // xlo, xhi: images of the clamps of the u scale
  double xt10 = 0;                             // t copula: upper-tail quantile of 1e-10 (the clamp of rvinecopulib's arguments)
  std::vector<double> C;                       // NP x (D + 1) Chebyshev coefficients
  // optional tables of d ln G / d nu and d ln G / d nu2 (finite differences of the exact composition)
  std::vector<double> C1, C2;
  bool withD = false;
  double lc1_nu = 0, lc1_nu2 = 0;               // d ln c1 / d nu and / d nu2

  void build(double nu_, double nu2_, bool withD_ = false) {
    nu = nu_; nu2 = nu2_; withD = withD_;
    c1 = R::dt(0.0, nu + 1.0, 0) / R::dt(0.0, nu2, 0);    // G(z) = c1 z for small z
    zlo = std::exp(SMIN); zhi = std::exp(SMAX);
    xlo = ast_q(1e-12, nu2); xhi = ast_q(1.0 - 1e-12, nu2);
    xt10 = R::qt(1e-10, nu2, 0, 0);
    C.assign(NP * (D + 1), 0.0);
    if (withD) { C1.assign(NP * (D + 1), 0.0); C2.assign(NP * (D + 1), 0.0); }
    double node[D + 1], cs[D + 1][D + 1];
    for (int k = 0; k <= D; k++) {
      node[k] = std::cos(M_PI * (k + 0.5) / (D + 1));
      for (int j = 0; j <= D; j++) cs[j][k] = std::cos(j * M_PI * (k + 0.5) / (D + 1));
    }
    for (int p = 0; p < NP; p++) {
      const double a = SMIN + p * W;
      double y[D + 1];
      for (int k = 0; k <= D; k++) y[k] = std::log(ast_Gexact(std::exp(a + (node[k] + 1.0) / 2.0 * W), nu, nu2));
      for (int j = 0; j <= D; j++) {
        double s = 0.0;
        for (int k = 0; k <= D; k++) s += y[k] * cs[j][k];
        C[p * (D + 1) + j] = 2.0 / (D + 1) * s;
      }
      C[p * (D + 1)] *= 0.5;
      if (withD) {
        const double h1 = 1e-4 * nu, h2 = 1e-4 * nu2;
        double y1[D + 1], y2[D + 1];
        for (int k = 0; k <= D; k++) {
          const double z = std::exp(a + (node[k] + 1.0) / 2.0 * W);
          y1[k] = (std::log(ast_Gexact(z, nu + h1, nu2)) - std::log(ast_Gexact(z, nu - h1, nu2))) / (2.0 * h1);
          y2[k] = (std::log(ast_Gexact(z, nu, nu2 + h2)) - std::log(ast_Gexact(z, nu, nu2 - h2))) / (2.0 * h2);
        }
        for (int j = 0; j <= D; j++) {
          double s1 = 0.0, s2 = 0.0;
          for (int k = 0; k <= D; k++) { s1 += y1[k] * cs[j][k]; s2 += y2[k] * cs[j][k]; }
          C1[p * (D + 1) + j] = 2.0 / (D + 1) * s1;
          C2[p * (D + 1) + j] = 2.0 / (D + 1) * s2;
        }
        C1[p * (D + 1)] *= 0.5;
        C2[p * (D + 1)] *= 0.5;
      }
    }
    if (withD) {
      lc1_nu = 0.5 * R::digamma((nu + 2.0) / 2.0) - 0.5 * R::digamma((nu + 1.0) / 2.0) - 0.5 / (nu + 1.0);
      lc1_nu2 = -(0.5 * R::digamma((nu2 + 1.0) / 2.0) - 0.5 * R::digamma(nu2 / 2.0) - 0.5 / nu2);
    }
  }

  // G(z) and its derivatives with respect to nu and nu2 (build with withD = true)
  inline double evalD(double z, double &gnu, double &gnu2) const {
    if (z < zlo) { const double g = c1 * z; gnu = g * lc1_nu; gnu2 = g * lc1_nu2; return g; }
    if (z > zhi) {
      const double g = ast_Gexact(z, nu, nu2), h1 = 1e-4 * nu, h2 = 1e-4 * nu2;
      gnu = g * (std::log(ast_Gexact(z, nu + h1, nu2)) - std::log(ast_Gexact(z, nu - h1, nu2))) / (2.0 * h1);
      gnu2 = g * (std::log(ast_Gexact(z, nu, nu2 + h2)) - std::log(ast_Gexact(z, nu, nu2 - h2))) / (2.0 * h2);
      return g;
    }
    const double s = std::log(z);
    int p = (int)((s - SMIN) / W);
    if (p > NP - 1) p = NP - 1;
    const double t = 2.0 * (s - SMIN - p * W) / W - 1.0;
    const double *c = &C[p * (D + 1)], *d1 = &C1[p * (D + 1)], *d2 = &C2[p * (D + 1)];
    double b1 = 0.0, b2 = 0.0, e1 = 0.0, e2 = 0.0, f1 = 0.0, f2 = 0.0;
    for (int j = D; j >= 1; j--) {
      const double b0 = 2.0 * t * b1 - b2 + c[j]; b2 = b1; b1 = b0;
      const double e0 = 2.0 * t * e1 - e2 + d1[j]; e2 = e1; e1 = e0;
      const double f0 = 2.0 * t * f1 - f2 + d2[j]; f2 = f1; f1 = f0;
    }
    const double g = std::exp(t * b1 - b2 + c[0]);
    gnu = g * (t * e1 - e2 + d1[0]);
    gnu2 = g * (t * f1 - f2 + d2[0]);
    return g;
  }

  inline double eval(double z) const {
    if (z < zlo) return c1 * z;
    if (z > zhi) return ast_Gexact(z, nu, nu2);
    const double s = std::log(z);
    int p = (int)((s - SMIN) / W);
    if (p > NP - 1) p = NP - 1;
    const double t = 2.0 * (s - SMIN - p * W) / W - 1.0;
    const double *c = &C[p * (D + 1)];
    double b1 = 0.0, b2 = 0.0;
    for (int j = D; j >= 1; j--) { const double b0 = 2.0 * t * b1 - b2 + c[j]; b2 = b1; b1 = b0; }
    return std::exp(t * b1 - b2 + c[0]);
  }
};

}  // namespace
