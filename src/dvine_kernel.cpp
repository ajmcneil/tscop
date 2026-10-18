// C++ kernel for the d-vine likelihood of tied sdvinecop models: plain pair
// copulas (independence, Gauss, Clayton, Gumbel, Joe, Frank, all four
// rotations) and inverse-v-transformed (IVT) pair copulas with linear
// v-transforms and these bases, and IVT levels with nonlinear (v2p) v-transforms (src/nl.h; Clayton,
// Gumbel, Joe and Frank bases, rotations; the ast and Gaussian bases use the R path), evaluated by the same level recursion as the R
// reference implementation dv_loglik_r (R/dvine-engine.R). The recursion is a
// template over the scalar type: double gives the log-likelihood, and a
// forward-mode dual number with P derivative directions gives the
// log-likelihood and its gradient (see pair.h).
#include "pair.h"
#include "ast.h"
#include "kendall.h"
#include "nl.h"
#include <type_traits>

namespace {

// ---- absolute spherical t: scalar helpers with derivatives (double and dual numbers) ----

inline double lgamma_(double a) { return R::lgammafn(a); }
template <int P> inline Dual<P> lgamma_(const Dual<P> &a) { return chain(a, R::lgammafn(a.v), R::digamma(a.v)); }

// d pt_upper(z, k) / dk by a central difference (pt itself is accurate to 1e-15)
inline double pt_dk(double z, double k) {
  const double h = 1e-4 * k;
  return (R::pt(z, k + h, 0, 0) - R::pt(z, k - h, 0, 0)) / (2.0 * h);
}

// t-scale value x = qt_upper((1 - u) / 2, nu) of a u-scale value
inline double ast_entry(double u, double nu) { return ast_q(u, nu); }
template <int P> inline Dual<P> ast_entry(const Dual<P> &u, const Dual<P> &nu) {
  const double x = ast_q(u.v, nu.v), ft = R::dt(x, nu.v, 0);
  Dual<P> r; r.v = x;
  const double dxdu = 0.5 / ft, dxdnu = pt_dk(x, nu.v) / ft;   // implicit: pt_upper(x, nu) = (1 - u) / 2
  for (int k = 0; k < P; k++) r.d[k] = dxdu * u.d[k] + dxdnu * nu.d[k];
  return r;
}

// next level's t-scale value G(z) through the level map, clamped to the image of the u-scale clamps
inline double ast_carry(const AstMap &m, double z, double, double) {
  const double g = m.eval(z);
  return g < m.xlo ? m.xlo : (g > m.xhi ? m.xhi : g);
}
template <int P> inline Dual<P> ast_carry(const AstMap &m, const Dual<P> &z, const Dual<P> &nu, const Dual<P> &nu2) {
  double gnu, gnu2;
  const double g = m.evalD(z.v, gnu, gnu2);
  Dual<P> r;
  if (g < m.xlo) { r.v = m.xlo; return r; }
  if (g > m.xhi) { r.v = m.xhi; return r; }
  const double gz = R::dt(z.v, nu.v + 1.0, 0) / R::dt(g, nu2.v, 0);   // dG/dz = dt(z, nu + 1) / dt(G, nu2)
  r.v = g;
  for (int k = 0; k < P; k++) r.d[k] = gz * z.d[k] + gnu * nu.d[k] + gnu2 * nu2.d[k];
  return r;
}

// ---- t copula (rho in par, nu in bpar of a plain level): signed t-scale values x = qt(u, nu) ----
inline double tc_entry(double u, double nu) { return R::qt(clampin(u), nu, 1, 0); }
template <int P> inline Dual<P> tc_entry(const Dual<P> &u, const Dual<P> &nu) {
  const double uc = clampin(u.v);
  const double x = R::qt(uc, nu.v, 1, 0), ft = R::dt(x, nu.v, 0);
  Dual<P> r; r.v = x;
  const double dxdu = (uc == u.v) ? 1.0 / ft : 0.0, dxdnu = pt_dk(x, nu.v) / ft;
  for (int k = 0; k < P; k++) r.d[k] = dxdu * u.d[k] + dxdnu * nu.d[k];
  return r;
}
// next level's signed t-scale value: sign(w) G(|w|) through the same level map (G is odd), clamped
// to the image of rvinecopulib's clamp of its arguments
inline double tc_carry(const AstMap &m, double w, double, double) {
  double g = m.eval(std::fabs(w));
  if (g > m.xt10) g = m.xt10;
  return w < 0.0 ? -g : g;
}
template <int P> inline Dual<P> tc_carry(const AstMap &m, const Dual<P> &w, const Dual<P> &nu, const Dual<P> &nu2) {
  const double aw = std::fabs(w.v), sg = w.v < 0.0 ? -1.0 : 1.0;
  double gnu, gnu2;
  const double g = m.evalD(aw, gnu, gnu2);
  Dual<P> r;
  if (g > m.xt10) { r.v = sg * m.xt10; return r; }
  const double gz = R::dt(aw, nu.v + 1.0, 0) / R::dt(g, nu2.v, 0);
  r.v = sg * g;
  for (int k = 0; k < P; k++) r.d[k] = gz * w.d[k] + sg * (gnu * nu.d[k] + gnu2 * nu2.d[k]);
  return r;
}
// u-scale value pt(w, nu + 1) (leaving a run of t levels; not clamped, as the h-functions of rvinecopulib)
inline double tc_exit(double w, double nu) { return R::pt(w, nu + 1.0, 1, 0); }
template <int P> inline Dual<P> tc_exit(const Dual<P> &w, const Dual<P> &nu) {
  Dual<P> r; r.v = R::pt(w.v, nu.v + 1.0, 1, 0);
  const double dw = R::dt(w.v, nu.v + 1.0, 0), dn = -pt_dk(w.v, nu.v + 1.0);
  for (int k = 0; k < P; k++) r.d[k] = dw * w.d[k] + dn * nu.d[k];
  return r;
}

// u-scale value 1 - 2 pt_upper(z, nu + 1), clamped (leaving a run of ast levels)
inline double ast_exit(double z, double nu) { return clampu(1.0 - 2.0 * R::pt(z, nu + 1.0, 0, 0)); }
template <int P> inline Dual<P> ast_exit(const Dual<P> &z, const Dual<P> &nu) {
  const double u = 1.0 - 2.0 * R::pt(z.v, nu.v + 1.0, 0, 0);
  Dual<P> r;
  if (u < EPS) { r.v = EPS; return r; }
  if (u > 1.0 - EPS) { r.v = 1.0 - EPS; return r; }
  const double dz = 2.0 * R::dt(z.v, nu.v + 1.0, 0), dn = -2.0 * pt_dk(z.v, nu.v + 1.0);
  r.v = u;
  for (int k = 0; k < P; k++) r.d[k] = dz * z.d[k] + dn * nu.d[k];
  return r;
}

// The level recursion. Level j (0-based) has family code fam[j], rotation
// rot[j], parameter par[j]; an IVT level has ivt[j] = 1 and its base copula in
// bfam[j], brot[j], bpar[j]. The v-transform fulcrums are d1, d2.
template <class T>
T loglik_impl(const std::vector<T> &udat, const IntegerVector &fam, const IntegerVector &rot,
              const std::vector<T> &par, const IntegerVector &ivt, const IntegerVector &bfam,
              const IntegerVector &brot, const std::vector<T> &bpar, const T &d1, const T &d2,
              const IntegerVector &nl, const T &k1, const T &k2, std::vector<double> *tau_out = nullptr) {
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
    // Kendall's tau of the pairs fed to pair copula j, on the scale the arrays hold (normal scores, t
    // scale or base scale: tau is invariant to increasing transformations of each margin)
    auto rec = [&]() {
      if constexpr (std::is_same<T, double>::value) {
        if (tau_out) (*tau_out)[j] = kendall_tau_b(A.data(), B.data(), m);
      }
    };
    if (bs && !ivt[j]) {
      if (zb) {
        for (int i = 0; i < m; i++) { A[i] = pnorm_(A[i]); B[i] = pnorm_(B[i]); }
        zb = false;
      }
      for (int i = 0; i < m; i++) { A[i] = clampu(d1 - sA[i] * A[i]); B[i] = clampu(d2 - sB[i] * B[i]); }
      bs = false;
    }
    {
      // ---- absolute spherical t (ast) levels: plain (rotation 0) or the base of an IVT pair ----
      // t-scale values of a run of ast levels are carried through the level maps G (src/ast.h); for
      // dual numbers the maps also carry the derivatives with respect to the degrees of freedom.
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
        const T nu = ivt[j] ? bpar[j] : par[j];
        if (!ts) for (int i = 0; i < m; i++) { A[i] = ast_entry(A[i], nu); B[i] = ast_entry(B[i], nu); }
        rec();
        // log density: lgamma terms once, three log1p per pair
        const T c0 = lgamma_((nu + 2.0) / 2.0) - 2.0 * lgamma_((nu + 1.0) / 2.0) + lgamma_(nu / 2.0);
        const T e1 = -(nu + 2.0) / 2.0, e2 = (nu + 1.0) / 2.0;
        T s = c0 * double(m);
        for (int i = 0; i < m; i++) {
          const T x2 = A[i] * A[i] / nu, y2 = B[i] * B[i] / nu;
          s = s + e1 * log1p(x2 + y2) + e2 * (log1p(x2) + log1p(y2));
        }
        ll = ll + s;
        if (!last) {
          const bool carry = is_ast(j + 1) && (bool(ivt[j + 1]) == bool(ivt[j]));
          const T nu2 = ivt[j + 1] ? bpar[j + 1] : par[j + 1];
          AstMap amap;
          if (carry) amap.build(val(nu), val(nu2), !std::is_same<T, double>::value);
          const T k1 = nu + 1.0;
          for (int i = 0; i < m - 1; i++) {
            const T zA = A[i] * sqrt(k1 / (nu + B[i] * B[i]));                    // h2 at pair i
            const T zB = B[i + 1] * sqrt(k1 / (nu + A[i + 1] * A[i + 1]));        // h1 at pair i + 1
            if (carry) { h1[i] = ast_carry(amap, zA, nu, nu2); h2[i] = ast_carry(amap, zB, nu, nu2); }
            else { h1[i] = ast_exit(zA, nu); h2[i] = ast_exit(zB, nu); }
          }
          for (int i = 0; i < m - 1; i++) { A[i] = h1[i]; B[i] = h2[i]; }
          ts = carry;
          if (ivt[j]) for (int i = 0; i < m - 1; i++) sB[i] = sB[i + 1];
          m -= 1;
        }
        continue;
      }
    }
    {
      // ---- t copula levels (plain; rho = par, nu = bpar): signed t-scale values carried through level maps ----
      auto is_tc = [&](int jj) { return !ivt[jj] && fam[jj] == TCOP; };
      if (is_tc(j)) {
        if (zs) { for (int i = 0; i < m; i++) { A[i] = pnorm_(A[i]); B[i] = pnorm_(B[i]); } zs = false; }
        const T rho = par[j], nu = bpar[j];
        if (!ts) for (int i = 0; i < m; i++) { A[i] = tc_entry(A[i], nu); B[i] = tc_entry(B[i], nu); }
        rec();
        const T omr = 1.0 - rho * rho;
        const T c0 = lgamma_((nu + 2.0) / 2.0) + lgamma_(nu / 2.0) - 2.0 * lgamma_((nu + 1.0) / 2.0) - 0.5 * log1p(-rho * rho);
        const T e1 = -(nu + 2.0) / 2.0, e2 = (nu + 1.0) / 2.0;
        T s = c0 * double(m);
        for (int i = 0; i < m; i++) {
          const T Q = (A[i] * A[i] - 2.0 * rho * A[i] * B[i] + B[i] * B[i]) / omr;
          s = s + e1 * log1p(Q / nu) + e2 * (log1p(A[i] * A[i] / nu) + log1p(B[i] * B[i] / nu));
        }
        ll = ll + s;
        if (!last) {
          const bool carry = is_tc(j + 1);
          AstMap amap;
          if (carry) amap.build(val(nu), val(bpar[j + 1]), !std::is_same<T, double>::value);
          const T k1 = nu + 1.0;
          for (int i = 0; i < m - 1; i++) {
            const T wA = (A[i] - rho * B[i]) / sqrt((nu + B[i] * B[i]) * omr / k1);                 // h2 at pair i
            const T wB = (B[i + 1] - rho * A[i + 1]) / sqrt((nu + A[i + 1] * A[i + 1]) * omr / k1);   // h1 at pair i + 1
            if (carry) { h1[i] = tc_carry(amap, wA, nu, bpar[j + 1]); h2[i] = tc_carry(amap, wB, nu, bpar[j + 1]); }
            else { h1[i] = tc_exit(wA, nu); h2[i] = tc_exit(wB, nu); }
          }
          for (int i = 0; i < m - 1; i++) { A[i] = h1[i]; B[i] = h2[i]; }
          ts = carry;
          m -= 1;
        }
        continue;
      }
    }
    if (j < nl.size() && nl[j]) {
      // ---- nonlinear IVT level (v2p transforms, shared fulcrums d1, d2 and shapes k1, k2): evaluated on the u scale ----
      if (zs) {
        for (int i = 0; i < m; i++) { A[i] = pnorm_(A[i]); B[i] = pnorm_(B[i]); }
        zs = false;
      }
      rec();
      std::vector<T> tA(m), tB(m);
      for (int i = 0; i < m; i++) { tA[i] = v2p_T(A[i], d1, k1); tB[i] = v2p_T(B[i], d2, k2); }
      for (int i = 0; i < m; i++) ll = ll + nl_base_ld(bfam[j], brot[j], bpar[j], bpar[j], tA[i], tB[i]);
      if (!last) {
        for (int i = 0; i < m - 1; i++) {
          h2[i] = nl_cdf(bfam[j], brot[j], bpar[j], bpar[j], d1, k1, A[i], tB[i], true);          // P(U1 <= a_i | U2 = b_i)
          h1[i + 1] = nl_cdf(bfam[j], brot[j], bpar[j], bpar[j], d2, k2, B[i + 1], tA[i + 1], false);  // P(U2 <= b | U1 = a), pair i + 1
        }
        for (int i = 0; i < m - 1; i++) { A[i] = h2[i]; B[i] = h1[i + 1]; }
        m -= 1;
      }
      continue;
    }
    const bool plain_gauss = !ivt[j] && fam[j] == GAUSS;
    if (plain_gauss) {
      if (!zs) {
        for (int i = 0; i < m; i++) { A[i] = qnorm_(A[i]); B[i] = qnorm_(B[i]); }
        zs = true;
      }
      rec();
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
      rec();
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
      rec();
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
               const NumericVector &dd2, const IntegerVector &nl, double k1, double k2,
               const NumericVector &dk1, const NumericVector &dk2) {
  const int K = fam.size();
  std::vector<Dual<P> > p(K), bp(K);
  for (int j = 0; j < K; j++) {
    p[j].v = par[j]; bp[j].v = bpar[j];
    for (int k = 0; k < P; k++) { p[j].d[k] = dpar(j, k); bp[j].d[k] = dbpar(j, k); }
  }
  Dual<P> D1(d1), D2(d2), K1(k1), K2(k2);
  for (int k = 0; k < P; k++) {
    D1.d[k] = dd1[k]; D2.d[k] = dd2[k];
    if (dk1.size() == P) K1.d[k] = dk1[k];
    if (dk2.size() == P) K2.d[k] = dk2[k];
  }
  // the data enter with their own derivatives when du (n x P) is supplied (data that depend on the
  // parameters, for example a v-transformed series); otherwise they are constants
  const int n = u.size();
  std::vector<Dual<P> > ud(n);
  const bool dd = du.nrow() == n && du.ncol() == P;
  for (int i = 0; i < n; i++) {
    ud[i].v = u[i];
    if (dd) for (int k = 0; k < P; k++) ud[i].d[k] = du(i, k);
  }
  Dual<P> ll = loglik_impl<Dual<P> >(ud, fam, rot, p, ivt, bfam, brot, bp, D1, D2, nl, K1, K2);
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
                     double d1, double d2, Nullable<IntegerVector> nl = R_NilValue, double k1 = 1.0, double k2 = 1.0) {
  const int K = fam.size();
  std::vector<double> p(par.begin(), par.end()), bp(bpar.begin(), bpar.end());
  (void)K;
  const std::vector<double> ud(u.begin(), u.end());
  const IntegerVector nlv = nl.isNotNull() ? IntegerVector(nl) : IntegerVector(0);
  const double ll = loglik_impl<double>(ud, fam, rot, p, ivt, bfam, brot, bp, d1, d2, nlv, k1, k2);
  return std::isfinite(ll) ? ll : R_NegInf;
}

// Log-likelihood and gradient with respect to P (1 to 12) directions along which the
// level parameters (K x P matrices dpar, dbpar), the fulcrums (dd1, dd2) and, if du is an
// n x P matrix, the data u vary.
// [[Rcpp::export]]
List dv_loglik_grad_cpp(NumericVector u, NumericMatrix du, IntegerVector fam, IntegerVector rot, NumericVector par,
                        NumericMatrix dpar, IntegerVector ivt, IntegerVector bfam, IntegerVector brot,
                        NumericVector bpar, NumericMatrix dbpar, double d1, double d2,
                        NumericVector dd1, NumericVector dd2, Nullable<IntegerVector> nl = R_NilValue,
                        double k1 = 1.0, double k2 = 1.0, Nullable<NumericVector> dk1 = R_NilValue,
                        Nullable<NumericVector> dk2 = R_NilValue) {
  const int P = dpar.ncol();
  const IntegerVector nlv = nl.isNotNull() ? IntegerVector(nl) : IntegerVector(0);
  const NumericVector dk1v = dk1.isNotNull() ? NumericVector(dk1) : NumericVector(0);
  const NumericVector dk2v = dk2.isNotNull() ? NumericVector(dk2) : NumericVector(0);
  switch (P) {
#define CASE(PP) case PP: return grad_impl<PP>(u, du, fam, rot, par, dpar, ivt, bfam, brot, bpar, dbpar, d1, d2, dd1, dd2, nlv, k1, k2, dk1v, dk2v);
    CASE(1) CASE(2) CASE(3) CASE(4) CASE(5) CASE(6) CASE(7) CASE(8) CASE(9) CASE(10) CASE(11) CASE(12)
#undef CASE
    default: stop("the kernel gradient supports 1 to 12 directions");
  }
}

// Kendall's tau of the level pairs of the d-vine recursion (the empirical Kendall partial
// autocorrelations given the pair copulas): element j is the tau of the pairs fed to pair copula
// j + 1, NA beyond the last level the data length allows.
// [[Rcpp::export]]
NumericVector dv_level_tau_cpp(NumericVector u, IntegerVector fam, IntegerVector rot, NumericVector par,
                               IntegerVector ivt, IntegerVector bfam, IntegerVector brot, NumericVector bpar,
                               double d1, double d2) {
  std::vector<double> p(par.begin(), par.end()), bp(bpar.begin(), bpar.end());
  const std::vector<double> ud(u.begin(), u.end());
  std::vector<double> tau(fam.size(), NA_REAL);
  loglik_impl<double>(ud, fam, rot, p, ivt, bfam, brot, bp, d1, d2, IntegerVector(0), 1.0, 1.0, &tau);
  return wrap(tau);
}
