// Kalman filter for the Gaussian AR-factor copula (arfactcop): the log-likelihood of
// the normal scores x of the model x_t = sum_i sqrt(w_i) a_i,t + sqrt(w0) e_t with
// independent unit-variance AR(1) factors a_i (rates lambda_i). A template over the
// scalar type: double gives the log-likelihood (R reference: afc_filter in
// R/arfactcop.R), a forward-mode dual number gives its gradient with respect to the
// parameters and the data. Same steady-state shortcut as the R code.
#include "pair.h"

namespace {

template <class T>
T afc_ll_impl(const std::vector<T> &x, const std::vector<char> &miss, const std::vector<T> &lam,
              const std::vector<T> &w, const T &w0, double tol) {
  using std::sqrt; using std::log;
  const int n = x.size(), m = lam.size();
  std::vector<T> zv(m), a(m), Pz(m), Kg(m);
  std::vector<T> P(m * m), Pn(m * m), L2(m * m), Q(m * m);
  for (int i = 0; i < m; i++) {
    zv[i] = sqrt(w[i]);
    a[i] = T(0.0);
    for (int j = 0; j < m; j++) {
      L2[i * m + j] = lam[i] * lam[j];
      P[i * m + j] = T(i == j ? 1.0 : 0.0);
      Q[i * m + j] = i == j ? 1.0 - lam[i] * lam[i] : T(0.0);
    }
  }
  const double lt2pi = std::log(2.0 * M_PI);
  T ll(0.0);
  bool steady = false;
  for (int t = 0; t < n; t++) {
    for (int i = 0; i < m; i++) {
      T s(0.0);
      for (int j = 0; j < m; j++) s = s + P[i * m + j] * zv[j];
      Pz[i] = s;
    }
    T Ft = w0, mu(0.0);
    for (int i = 0; i < m; i++) { Ft = Ft + zv[i] * Pz[i]; mu = mu + zv[i] * a[i]; }
    if (miss[t]) {
      for (int i = 0; i < m; i++) a[i] = lam[i] * a[i];
      for (int k = 0; k < m * m; k++) P[k] = L2[k] * P[k] + Q[k];
      steady = false;
    } else {
      const T vt = x[t] - mu;
      ll = ll - 0.5 * (lt2pi + log(Ft) + vt * vt / Ft);
      if (!steady) for (int i = 0; i < m; i++) Kg[i] = Pz[i] / Ft;
      for (int i = 0; i < m; i++) a[i] = lam[i] * (a[i] + Kg[i] * vt);
      if (!steady) {
        double dmax = 0.0;
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            const int k = i * m + j;
            Pn[k] = L2[k] * (P[k] - Pz[i] * Pz[j] / Ft) + Q[k];
            dmax = std::max(dmax, std::fabs(val(Pn[k]) - val(P[k])));
          }
        }
        steady = dmax < tol;
        for (int k = 0; k < m * m; k++) P[k] = Pn[k];
      }
    }
  }
  return ll;
}

template <int P>
List afc_grad_impl(const NumericVector &x, const NumericMatrix &dx, const NumericVector &lambda,
                   const NumericMatrix &dlam, const NumericVector &w, const NumericMatrix &dw, double w0,
                   const NumericVector &dw0, double tol) {
  const int n = x.size(), m = lambda.size();
  std::vector<Dual<P> > xd(n), lam(m), ww(m);
  std::vector<char> miss(n);
  const bool dd = dx.nrow() == n && dx.ncol() == P;
  for (int t = 0; t < n; t++) {
    miss[t] = std::isnan(x[t]);
    xd[t].v = miss[t] ? 0.0 : x[t];
    if (dd && !miss[t]) for (int k = 0; k < P; k++) xd[t].d[k] = dx(t, k);
  }
  for (int i = 0; i < m; i++) {
    lam[i].v = lambda[i]; ww[i].v = w[i];
    for (int k = 0; k < P; k++) { lam[i].d[k] = dlam(i, k); ww[i].d[k] = dw(i, k); }
  }
  Dual<P> W0(w0);
  for (int k = 0; k < P; k++) W0.d[k] = dw0[k];
  Dual<P> ll = afc_ll_impl<Dual<P> >(xd, miss, lam, ww, W0, tol);
  NumericVector g(P);
  for (int k = 0; k < P; k++) g[k] = ll.d[k];
  return List::create(_["value"] = ll.v, _["gradient"] = g);
}

}  // namespace

// Gaussian log-likelihood of the normal scores x (NaN for missing values).
// [[Rcpp::export]]
double afc_loglik_cpp(NumericVector x, NumericVector lambda, NumericVector w, double w0, double tol) {
  const int n = x.size(), m = lambda.size();
  std::vector<double> xd(n), lam(lambda.begin(), lambda.end()), ww(w.begin(), w.end());
  std::vector<char> miss(n);
  for (int t = 0; t < n; t++) { miss[t] = std::isnan(x[t]); xd[t] = miss[t] ? 0.0 : x[t]; }
  (void)m;
  return afc_ll_impl<double>(xd, miss, lam, ww, w0, tol);
}

// Log-likelihood and gradient along P (1 to 12) directions: the derivatives of the rates
// (m x P), weights (m x P), noise weight (P) and, if dx is an n x P matrix, of the data.
// [[Rcpp::export]]
List afc_loglik_grad_cpp(NumericVector x, NumericMatrix dx, NumericVector lambda, NumericMatrix dlam, NumericVector w,
                         NumericMatrix dw, double w0, NumericVector dw0, double tol) {
  const int P = dlam.ncol();
  switch (P) {
#define CASE(PP) case PP: return afc_grad_impl<PP>(x, dx, lambda, dlam, w, dw, w0, dw0, tol);
    CASE(1) CASE(2) CASE(3) CASE(4) CASE(5) CASE(6) CASE(7) CASE(8) CASE(9) CASE(10) CASE(11) CASE(12)
#undef CASE
    default: stop("the gradient supports 1 to 12 directions");
  }
}
