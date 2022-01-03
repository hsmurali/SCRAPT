#ifndef TERNARY_SORT_HPP
#define TERNARY_SORT_HPP

#include <algorithm>

namespace ternary_sort
{

  void vecswap(int i, int j, int n, const char *x[])
  {   
    while (n-- > 0) {
      std::swap(x[i], x[j]);
      i++;
      j++;
    }
  } // void vecswap(int i, int j, int n, const char *x[])

  void ssort1(const char *x[], int n, int depth)
  {
    int    a, b, c, d, r, v;
    if (n <= 1)
      return;
    // Randomness here means that outputs will not always be the same. Maybe.
    a = rand() % n;
    std::swap(x[0], x[a]);
    v = x[0][depth];
    a = b = 1;
    c = d = n-1;
    while (true) {
      while (b <= c && (r = x[b][depth]-v) <= 0) {
	if (r == 0) { std::swap(x[a], x[b]); a++; }
	b++;
      }
      while (b <= c && (r = x[c][depth]-v) >= 0) {
	if (r == 0) { std::swap(x[c], x[d]); d--; }
	c--;
      }
      if (b > c) break;
      std::swap(x[b], x[c]);
      b++;
      c--;
    }
    r = std::min(a, b-a);     vecswap(0, b-r, r, x);
    r = std::min(d-c, n-d-1); vecswap(b, n-r, r, x);
    r = b-a; ssort1(x, r, depth);
    if (x[r][depth] != 0)
      ssort1(x + r, a + n-d-1, depth+1);
    r = d-c; ssort1(x + n-r, r, depth);
  } // void ssort1(const char *x[], int n, int depth)

  void ternarySort(const char *x[], int n)
  { 
    ssort1(x, n, 0); 
  } // void ternarySort(const char *x[], int n)

} // namespace ternary_sort

#endif // TERNARY_SORT_HPP
