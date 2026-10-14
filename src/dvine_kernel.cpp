// C++ kernel for the d-vine likelihood of tied sdvinecop models: plain pair
// copulas (independence, Gauss, Clayton, Gumbel, Joe, Frank, all four
// rotations) and inverse-v-transformed (IVT) pair copulas with linear
// v-transforms and these bases, evaluated by the same level recursion as the R
// reference implementation dv_loglik_r (R/dvine-engine.R). The recursion is a
// template over the scalar type: double gives the log-likelihood, and a
// forward-mode dual number with P derivative directions gives the
// log-likelihood and its gradient (see pair.h).
#include "pair.h"
#include "ast.h"
#include <type_traits>

namespace {

// The level recursion. Level j (0-based) has family code fam[j], rotation
// rot[j], parameter par[j]; an IVT level has ivt[j] = 1 and its base copula in
// bfam[j], brot[j], bpar[j]. The v-transform fulcrums are d1, d2.
template <class T>
T loglik_impl(const std::vector<T> &udat, const IntegerVector &fam, const IntegerVector &rot,
              const std::vector<T> &par, const IntegerVector &ivt, const IntegerVector &bfam,
              const IntegerVector &brot, const std::vector<T> &bpar, const T &d1, const T &d2) {
  using std::sqrt; using std::log1p;
  const int n = udat.size();
  const int K = std::min<int>(fam.size(), n - 1);
  if (K < 1) return T(0.0);
  int m = n - 1;
  std::vector<T> A(udat.begin(), udat.end() - 1), B(udat.begin() + 1, udat.end());
  std::vector<T> sA, sB, h1(m), h2(m);
  bool zs = false, bs = false, zb = false;
  bool ts = false;   // the arrays hold t-scale values of the current ast level (double only)
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
    if constexpr (std::is_same<T, double>::value) {
      // ---- absolute spherical t (ast) levels: plain (rotation 0) or the base of an IVT pair ----
      auto is_ast = [&](int jj) { return ivt[jj] ? bfam[jj] == AST : fam[jj] == AST; };
      if (is_ast(j)) {
        if (ivt[j]) {
          if (!bs) {
            if (zs) { for (int i = 0; i < m; i++) { A[i] = pnorm_(A[i]); B[i] = pnorm_(B[i]); } zs = false; }
            sA.assign(m, T(0.0)); sB.assign(m, T(0.0));
            for (int i = 0; i < m; i++) {
              sA[i] = alin(A[i], d1); sB[i] = alin(B[i], d2);
              A[i] = vlin(A[i], d1); B[i] = vlin(B[i], d2);
            }
            bs = true; zb = false;
          }
          if (zb) { for (int i = 0; i < m; i++) { A[i] = pnorm_(A[i]); B[i] = pnorm_(B[i]); } zb = false; }
        } else if (zs) {
          for (int i = 0; i < m; i++) { A[i] = pnorm_(A[i]); B[i] = pnorm_(B[i]); }
          zs = false;
        }
        const double nu = ivt[j] ? bpar[j] : par[j];
        if (!ts) for (int i = 0; i < m; i++) { A[i] = ast_q(A[i], nu); B[i] = ast_q(B[i], nu); }
        // log density: lgamma terms once, three log1p per pair
        const double c0 = R::lgammafn((nu + 2.0) / 2.0) - 2.0 * R::lgammafn((nu + 1.0) / 2.0) + R::lgammafn(nu / 2.0);
        double s = 0.0;
        for (int i = 0; i < m; i++) {
          const double x2 = A[i] * A[i] / nu, y2 = B[i] * B[i] / nu;
          s += c0 - (nu + 2.0) / 2.0 * std::log1p(x2 + y2) + (nu + 1.0) / 2.0 * (std::log1p(x2) + std::log1p(y2));
        }
        ll = ll + s;
        if (!last) {
          const bool carry = is_ast(j + 1) && (bool(ivt[j + 1]) == bool(ivt[j]));
          AstMap amap;
          if (carry) amap.build(nu, ivt[j + 1] ? bpar[j + 1] : par[j + 1]);
          const double k1 = nu + 1.0;
          for (int i = 0; i < m - 1; i++) {
            const double zA = A[i] * std::sqrt(k1 / (nu + B[i] * B[i]));          // h2 at pair i
            const double zB = B[i + 1] * std::sqrt(k1 / (nu + A[i + 1] * A[i + 1]));   // h1 at pair i + 1
            if (carry) {
              h1[i] = std::min(std::max(amap.eval(zA), amap.xlo), amap.xhi);
              h2[i] = std::min(std::max(amap.eval(zB), amap.xlo), amap.xhi);
            } else {
              h1[i] = clampu(1.0 - 2.0 * R::pt(zA, k1, 0, 0));
              h2[i] = clampu(1.0 - 2.0 * R::pt(zB, k1, 0, 0));
            }
          }
          for (int i = 0; i < m - 1; i++) { A[i] = h1[i]; B[i] = h2[i]; }
          ts = carry;
          if (ivt[j]) for (int i = 0; i < m - 1; i++) sB[i] = sB[i + 1];
          m -= 1;
        }
        continue;
      }
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
          pair_eval(bfam[j], brot[j], bpar[j], bpar[j], A[i], B[i], !last, true, lc, h1[i], h2[i]);
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
        pair_eval(fam[j], rot[j], par[j], bpar[j], A[i], B[i], !last, false, lc, h1[i], h2[i]);   // bpar[j]: the second parameter of BB1
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
List grad_impl(const NumericVector &u, const NumericMatrix &du, const IntegerVector &fam, const IntegerVector &rot,
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
  // the data enter with their own derivatives when du (n x P) is supplied (data that depend on the
  // parameters, for example a v-transformed series); otherwise they are constants
  const int n = u.size();
  std::vector<Dual<P> > ud(n);
  const bool dd = du.nrow() == n && du.ncol() == P;
  for (int i = 0; i < n; i++) {
    ud[i].v = u[i];
    if (dd) for (int k = 0; k < P; k++) ud[i].d[k] = du(i, k);
  }
  Dual<P> ll = loglik_impl<Dual<P> >(ud, fam, rot, p, ivt, bfam, brot, bp, D1, D2);
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
  const std::vector<double> ud(u.begin(), u.end());
  const double ll = loglik_impl<double>(ud, fam, rot, p, ivt, bfam, brot, bp, d1, d2);
  return std::isfinite(ll) ? ll : R_NegInf;
}

// Log-likelihood and gradient with respect to P (1 to 12) directions along which the
// level parameters (K x P matrices dpar, dbpar), the fulcrums (dd1, dd2) and, if du is an
// n x P matrix, the data u vary.
// [[Rcpp::export]]
List dv_loglik_grad_cpp(NumericVector u, NumericMatrix du, IntegerVector fam, IntegerVector rot, NumericVector par,
                        NumericMatrix dpar, IntegerVector ivt, IntegerVector bfam, IntegerVector brot,
                        NumericVector bpar, NumericMatrix dbpar, double d1, double d2,
                        NumericVector dd1, NumericVector dd2) {
  const int P = dpar.ncol();
  switch (P) {
#define CASE(PP) case PP: return grad_impl<PP>(u, du, fam, rot, par, dpar, ivt, bfam, brot, bpar, dbpar, d1, d2, dd1, dd2);
    CASE(1) CASE(2) CASE(3) CASE(4) CASE(5) CASE(6) CASE(7) CASE(8) CASE(9) CASE(10) CASE(11) CASE(12)
#undef CASE
    default: stop("the kernel gradient supports 1 to 12 directions");
  }
}
