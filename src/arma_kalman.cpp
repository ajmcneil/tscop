// Kalman filter for the zero-mean, unit-variance Gaussian ARMA process behind an
// armacop model (Harvey's state-space form, innovation variance 1 / gamma_0), as a
// template over the scalar type: double gives the log-likelihood (R reference:
// arma_filter in R/arma-filter.R), a forward-mode dual number gives its gradient with
// respect to the (expanded) AR and MA coefficients and to the data. The same steps as
// the R code: gamma_0 from the (p + 1) x (p + 1) system, the stationary state
// covariance by doubling, then the filter with the steady-state shortcut.
#include "pair.h"

namespace {

template <class T>
T arma_ll_impl(const std::vector<T> &x, const std::vector<char> &miss, const std::vector<T> &ar,
               const std::vector<T> &ma, double tol) {
  using std::sqrt; using std::log;
  const int n = x.size(), p = ar.size(), q = ma.size();
  const int m = std::max(p, q + 1);
  // gamma_0: psi weights, c_k, and the linear system for gamma_0, ..., gamma_p
  std::vector<T> psi(q + 1), th(q + 1);
  psi[0] = T(1.0);
  for (int j = 1; j <= q; j++) {
    T s = ma[j - 1];
    for (int i = 1; i <= std::min(j, p); i++) s = s + ar[i - 1] * psi[j - i];
    psi[j] = s;
  }
  th[0] = T(1.0);
  for (int j = 1; j <= q; j++) th[j] = ma[j - 1];
  auto ck = [&](int k) {
    T s(0.0);
    if (k > q) return s;
    for (int i = 0; i <= q - k; i++) s = s + th[k + i] * psi[i];
    return s;
  };
  const int np = p + 1;
  std::vector<T> A(np * np), rhs(np), gam(np);
  for (int i = 0; i < np; i++) for (int j = 0; j < np; j++) A[i * np + j] = T(i == j ? 1.0 : 0.0);
  for (int k = 0; k <= p; k++) {
    for (int i = 1; i <= p; i++) { const int col = std::abs(k - i); A[k * np + col] = A[k * np + col] - ar[i - 1]; }
    rhs[k] = ck(k);
  }
  for (int c = 0; c < np; c++) {   // Gaussian elimination with partial pivoting (on the values)
    int piv = c;
    for (int r = c + 1; r < np; r++) if (std::fabs(val(A[r * np + c])) > std::fabs(val(A[piv * np + c]))) piv = r;
    if (piv != c) { for (int j = 0; j < np; j++) std::swap(A[c * np + j], A[piv * np + j]); std::swap(rhs[c], rhs[piv]); }
    for (int r = c + 1; r < np; r++) {
      const T f = A[r * np + c] / A[c * np + c];
      for (int j = c; j < np; j++) A[r * np + j] = A[r * np + j] - f * A[c * np + j];
      rhs[r] = rhs[r] - f * rhs[c];
    }
  }
  for (int r = np - 1; r >= 0; r--) {
    T s = rhs[r];
    for (int j = r + 1; j < np; j++) s = s - A[r * np + j] * gam[j];
    gam[r] = s / A[r * np + r];
  }
  const T g0 = gam[0];
  // state space: phi, R, Tm (first column phi, superdiagonal ones), RR
  std::vector<T> phi(m, T(0.0)), R(m, T(0.0));
  for (int i = 0; i < p; i++) phi[i] = ar[i];
  const T sg = sqrt(g0);
  R[0] = T(1.0) / sg;
  for (int j = 1; j <= q; j++) R[j] = ma[j - 1] / sg;
  std::vector<T> RR(m * m), Pm(m * m), Am(m * m), tmp(m * m), dP(m * m), A2(m * m);
  for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) RR[i * m + j] = R[i] * R[j];
  // stationary covariance by doubling: P = RR, A = Tm; dP = A P A'; P += dP; A = A A
  for (int i = 0; i < m * m; i++) { Pm[i] = RR[i]; Am[i] = T(0.0); }
  for (int i = 0; i < m; i++) { Am[i * m] = phi[i]; if (i + 1 < m) Am[i * m + i + 1] = T(1.0); }
  for (int it = 0; it < 200; it++) {
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) {      // tmp = A P
      T s(0.0); for (int k = 0; k < m; k++) s = s + Am[i * m + k] * Pm[k * m + j]; tmp[i * m + j] = s; }
    double dmax = 0.0, pmax = 0.0;
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) {      // dP = tmp A'
      T s(0.0); for (int k = 0; k < m; k++) s = s + tmp[i * m + k] * Am[j * m + k]; dP[i * m + j] = s;
      dmax = std::max(dmax, std::fabs(val(s))); }
    for (int i = 0; i < m * m; i++) { Pm[i] = Pm[i] + dP[i]; pmax = std::max(pmax, std::fabs(val(Pm[i]))); }
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) {      // A = A A
      T s(0.0); for (int k = 0; k < m; k++) s = s + Am[i * m + k] * Am[k * m + j]; A2[i * m + j] = s; }
    for (int i = 0; i < m * m; i++) Am[i] = A2[i];
    if (dmax <= 1e-16 * pmax) break;
  }
  // filter
  std::vector<T> a(m, T(0.0)), an(m), Kg(m, T(0.0)), Pc(m), TP(m * m), Pn(m * m);
  T ll(0.0);
  const double lt2pi = std::log(2.0 * M_PI);
  bool steady = false;
  // T P T' for a symmetric matrix P using the structure of Tm
  auto TPT = [&](const std::vector<T> &Q, std::vector<T> &out) {
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++)
      TP[i * m + j] = phi[i] * Q[j] + (i + 1 < m ? Q[(i + 1) * m + j] : T(0.0));
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++)
      out[i * m + j] = TP[i * m] * phi[j] + (j + 1 < m ? TP[i * m + j + 1] : T(0.0));
  };
  auto Tstep = [&](const std::vector<T> &s, std::vector<T> &out) {
    for (int i = 0; i < m; i++) out[i] = phi[i] * s[0] + (i + 1 < m ? s[i + 1] : T(0.0));
  };
  for (int t = 0; t < n; t++) {
    const T Ft = Pm[0];
    if (miss[t]) {
      Tstep(a, an); a = an;
      TPT(Pm, Pn);
      for (int k = 0; k < m * m; k++) Pm[k] = Pn[k] + RR[k];
      steady = false;
    } else {
      const T vt = x[t] - a[0];
      ll = ll - 0.5 * (lt2pi + log(Ft) + vt * vt / Ft);
      if (!steady) for (int i = 0; i < m; i++) Kg[i] = Pm[i * m] / Ft;
      std::vector<T> s(m);
      for (int i = 0; i < m; i++) s[i] = a[i] + Kg[i] * vt;
      Tstep(s, an); a = an;
      if (!steady) {
        std::vector<T> Q(m * m);
        for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) Q[i * m + j] = Pm[i * m + j] - Pm[i * m] * Pm[j * m] / Ft;
        TPT(Q, Pn);
        double dmax = 0.0;
        for (int k = 0; k < m * m; k++) { Pn[k] = Pn[k] + RR[k]; dmax = std::max(dmax, std::fabs(val(Pn[k]) - val(Pm[k]))); }
        steady = dmax < tol;
        Pm = Pn;
      }
    }
  }
  return ll;
}

template <int P>
List arma_grad_impl(const NumericVector &x, const NumericMatrix &dx, const NumericVector &ar, const NumericMatrix &dar,
                    const NumericVector &ma, const NumericMatrix &dma, double tol) {
  const int n = x.size(), p = ar.size(), q = ma.size();
  std::vector<Dual<P> > xd(n), a(p), b(q);
  std::vector<char> miss(n);
  const bool dd = dx.nrow() == n && dx.ncol() == P;
  for (int t = 0; t < n; t++) {
    miss[t] = std::isnan(x[t]);
    xd[t].v = miss[t] ? 0.0 : x[t];
    if (dd && !miss[t]) for (int k = 0; k < P; k++) xd[t].d[k] = dx(t, k);
  }
  for (int i = 0; i < p; i++) { a[i].v = ar[i]; for (int k = 0; k < P; k++) a[i].d[k] = dar(i, k); }
  for (int i = 0; i < q; i++) { b[i].v = ma[i]; for (int k = 0; k < P; k++) b[i].d[k] = dma(i, k); }
  Dual<P> ll = arma_ll_impl<Dual<P> >(xd, miss, a, b, tol);
  NumericVector g(P);
  for (int k = 0; k < P; k++) g[k] = ll.d[k];
  return List::create(_["value"] = ll.v, _["gradient"] = g);
}

}  // namespace

// Gaussian log-likelihood of the normal scores x (NaN for missing values) of the ARMA
// process with expanded coefficients ar and ma.
// [[Rcpp::export]]
double arma_loglik_cpp(NumericVector x, NumericVector ar, NumericVector ma, double tol) {
  const int n = x.size();
  std::vector<double> xd(n), a(ar.begin(), ar.end()), b(ma.begin(), ma.end());
  std::vector<char> miss(n);
  for (int t = 0; t < n; t++) { miss[t] = std::isnan(x[t]); xd[t] = miss[t] ? 0.0 : x[t]; }
  return arma_ll_impl<double>(xd, miss, a, b, tol);
}

// Log-likelihood and gradient along P (1 to 12) directions: the derivatives of the AR
// coefficients (p x P), MA coefficients (q x P) and, if dx is an n x P matrix, of the data.
// [[Rcpp::export]]
List arma_loglik_grad_cpp(NumericVector x, NumericMatrix dx, NumericVector ar, NumericMatrix dar, NumericVector ma,
                          NumericMatrix dma, double tol) {
  const int P = dar.ncol();
  switch (P) {
#define CASE(PP) case PP: return arma_grad_impl<PP>(x, dx, ar, dar, ma, dma, tol);
    CASE(1) CASE(2) CASE(3) CASE(4) CASE(5) CASE(6) CASE(7) CASE(8) CASE(9) CASE(10) CASE(11) CASE(12)
#undef CASE
    default: stop("the gradient supports 1 to 12 directions");
  }
}
