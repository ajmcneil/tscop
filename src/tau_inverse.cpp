// Vectorised Kendall's tau inverses for the Joe and Frank copulas (see
// R/tau-inverse.R for the formulas and the R reference versions).
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

namespace {

// Joe: tau = 1 - (2 / theta) sum_n c_n r^(n-1), r = 1 - 2 / theta
void joe_tau_d(double th, const NumericVector &c, double &tau, double &dtau) {
  const double r = 1.0 - 2.0 / th;
  double s = 0.0, ds = 0.0, rp = 1.0;   // rp = r^(n-1)
  const int N = c.size();
  for (int n = 1; n <= N; n++) {
    s += c[n - 1] * rp;
    if (n >= 2) ds += c[n - 1] * (n - 1) * (rp / r);   // c_n (n-1) r^(n-2)
    rp *= r;
    if (std::fabs(rp) * c[n - 1] < 1e-19 && n > 3) break;
  }
  if (r == 0.0) ds = c[1];   // limit: c_2 * 1 * r^0 at n = 2
  tau = 1.0 - 2.0 / th * s;
  dtau = 2.0 / (th * th) * s - 2.0 / th * ds * 2.0 / (th * th);
}

// Frank, theta > 0: series for theta <= 2, exponential sum beyond
void frank_tau_d(double th, const NumericVector &f, double &tau, double &dtau) {
  if (th <= 2.0) {
    const double t2 = th * th;
    double s = 0.0, ds = 0.0, p = 1.0;   // p = t2^(k-1)
    const int N = f.size();
    for (int k = 1; k <= N; k++) {
      s += f[k - 1] * p;
      ds += (2 * k - 1) * f[k - 1] * p;
      p *= t2;
    }
    tau = -4.0 * s * th;
    dtau = -4.0 * ds;
  } else {
    double tail = 0.0;
    for (int m = 1; m <= 200; m++) {
      const double e = std::exp(-m * th);
      tail += e * (th / m + 1.0 / (double(m) * m));
      if (e < 1e-19) break;
    }
    const double D1 = (M_PI * M_PI / 6.0 - tail) / th;
    tau = 1.0 - 4.0 / th * (1.0 - D1);
    dtau = 4.0 / (th * th) * (1.0 - 2.0 * D1) + 4.0 / (th * std::expm1(th));
  }
}

}  // namespace

// [[Rcpp::export]]
List joe_tau_cpp(NumericVector theta, NumericVector c) {
  const int n = theta.size();
  NumericVector tau(n), dtau(n);
  for (int i = 0; i < n; i++) joe_tau_d(theta[i], c, tau[i], dtau[i]);
  return List::create(_["tau"] = tau, _["dtau"] = dtau);
}

// [[Rcpp::export]]
List frank_tau_cpp(NumericVector theta, NumericVector f) {
  const int n = theta.size();
  NumericVector tau(n), dtau(n);
  for (int i = 0; i < n; i++) {
    const double a = std::fabs(theta[i]);
    frank_tau_d(a, f, tau[i], dtau[i]);
    if (theta[i] < 0) tau[i] = -tau[i];
  }
  return List::create(_["tau"] = tau, _["dtau"] = dtau);
}

// Joe parameter at Kendall's tau in [0, 1)
// [[Rcpp::export]]
NumericVector joe_par_cpp(NumericVector tau, NumericVector c) {
  const int n = tau.size();
  NumericVector out(n, 1.0);
  for (int i = 0; i < n; i++) {
    const double t = tau[i];
    if (!(t > 1e-14)) continue;
    double lo = 1.0, hi = std::max(3.0, 4.0 / (1.0 - t));
    double tv, dv;
    for (;;) { joe_tau_d(hi, c, tv, dv); if (tv >= t) break; hi *= 2.0; }
    double x = 1.0 + 1.5 * t / (1.0 - 0.5 * t * t);
    if (x <= lo || x >= hi) x = 0.5 * (lo + hi);
    for (int it = 0; it < 100; it++) {
      joe_tau_d(x, c, tv, dv);
      const double fx = tv - t;
      if (fx > 0) hi = x; else lo = x;
      double xn = x - fx / dv;
      if (!std::isfinite(xn) || xn <= lo || xn >= hi) xn = 0.5 * (lo + hi);
      const bool done = std::fabs(xn - x) <= 1e-15 * (xn - 1.0);
      x = xn;
      if (done) break;
    }
    out[i] = x;
  }
  return out;
}

// Frank parameter at Kendall's tau in (-1, 1)
// [[Rcpp::export]]
NumericVector frank_par_cpp(NumericVector tau, NumericVector f) {
  const int n = tau.size();
  NumericVector out(n, 0.0);
  for (int i = 0; i < n; i++) {
    const double t = std::fabs(tau[i]);
    // |theta| is kept at 1e-10 or more: the Frank density is evaluated through
    // expm1 and degenerates at exactly 0 (tau of 1e-11 or less is independence)
    if (t < 1.1e-11) { out[i] = tau[i] < 0 ? -1e-10 : 1e-10; continue; }
    double lo = 0.0, hi = std::max(10.0, 8.0 / (1.0 - t));
    double tv, dv;
    for (;;) { frank_tau_d(hi, f, tv, dv); if (tv >= t) break; hi *= 2.0; }
    double x = std::min(9.0 * t, 0.5 * hi);
    for (int it = 0; it < 100; it++) {
      frank_tau_d(x, f, tv, dv);
      const double fx = tv - t;
      if (fx > 0) hi = x; else lo = x;
      double xn = x - fx / dv;
      if (!std::isfinite(xn) || xn <= lo || xn >= hi) xn = 0.5 * (lo + hi);
      const bool done = std::fabs(xn - x) <= 1e-15 * std::fabs(xn);
      x = xn;
      if (done) break;
    }
    out[i] = tau[i] < 0 ? -x : x;
  }
  return out;
}
