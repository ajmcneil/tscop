// Kendall's tau-b in O(m log m) (Knight 1966): sort by x, then count the inversions of y with a
// merge sort; ties in x, in y and in both are handled. Used for the empirical Kendall partial
// autocorrelations (ktauhat) and exported to R for testing.
#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace {

inline double kendall_tau_b(const double *x, const double *y, int m) {
  if (m < 2) return NAN;
  std::vector<int> idx(m);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int a, int b) { return x[a] < x[b] || (x[a] == x[b] && y[a] < y[b]); });
  const double n0 = 0.5 * double(m) * (m - 1);
  // ties in x and in both
  double n1 = 0.0, n3 = 0.0;
  {
    int i = 0;
    while (i < m) {
      int j = i;
      while (j + 1 < m && x[idx[j + 1]] == x[idx[i]]) j++;
      const double t = j - i + 1;
      n1 += 0.5 * t * (t - 1);
      int a = i;                                  // within the x tie, y is sorted: count joint ties
      while (a <= j) {
        int b = a;
        while (b + 1 <= j && y[idx[b + 1]] == y[idx[a]]) b++;
        const double u = b - a + 1;
        n3 += 0.5 * u * (u - 1);
        a = b + 1;
      }
      i = j + 1;
    }
  }
  // inversions of y in the order given by idx (merge sort, bottom up)
  std::vector<double> v(m), w(m);
  for (int i = 0; i < m; i++) v[i] = y[idx[i]];
  double dis = 0.0;
  for (int width = 1; width < m; width *= 2) {
    for (int lo = 0; lo < m; lo += 2 * width) {
      const int mid = std::min(lo + width, m), hi = std::min(lo + 2 * width, m);
      int a = lo, b = mid, k = lo;
      while (a < mid && b < hi) {
        if (v[b] < v[a]) { dis += mid - a; w[k++] = v[b++]; } else { w[k++] = v[a++]; }
      }
      while (a < mid) w[k++] = v[a++];
      while (b < hi) w[k++] = v[b++];
    }
    v.swap(w);
  }
  // ties in y (v is now sorted)
  double n2 = 0.0;
  {
    int i = 0;
    while (i < m) {
      int j = i;
      while (j + 1 < m && v[j + 1] == v[i]) j++;
      const double t = j - i + 1;
      n2 += 0.5 * t * (t - 1);
      i = j + 1;
    }
  }
  const double denom = std::sqrt((n0 - n1) * (n0 - n2));
  if (!(denom > 0.0)) return NAN;
  return (n0 - n1 - n2 + n3 - 2.0 * dis) / denom;
}

}  // namespace
