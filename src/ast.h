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
  std::vector<double> C;                       // NP x (D + 1) Chebyshev coefficients

  void build(double nu_, double nu2_) {
    nu = nu_; nu2 = nu2_;
    c1 = R::dt(0.0, nu + 1.0, 0) / R::dt(0.0, nu2, 0);    // G(z) = c1 z for small z
    zlo = std::exp(SMIN); zhi = std::exp(SMAX);
    xlo = ast_q(1e-12, nu2); xhi = ast_q(1.0 - 1e-12, nu2);
    C.assign(NP * (D + 1), 0.0);
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
    }
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
