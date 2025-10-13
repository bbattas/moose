#pragma once
#include <cstddef>      //size::t
#include "MooseError.h" // mooseAssert / mooseError (or use <cassert> assert())
#include <utility>      // std::pair
#include <cmath>        // std::sqrt

namespace GBPairPacking
{
/**
Pack (i,j) with 0 <= i <= j < N into a single index in the
upper-triangular (including diagonal) storage of size N*(N+1)/2.

Layout is column-by-column (j is the column):
  [ (0,0)
    (0,1) (1,1)
    (0,2) (1,2) (2,2)
    ... ]

Example:
  const unsigned T = PairPacking::count_upper_incl_diag(N);
  std::vector<RankTwoTensor> H(T);
  H[PairPacking::pack_upper_incl_diag(i, j)] = Bij;
*/
static inline constexpr unsigned
pack_upper_incl_diag(unsigned i, unsigned j) noexcept
{
  // caller is responsible for ensuring i <= j and j < N
  return j * (j + 1) / 2 + i;
}

// Optional helper: number of stored (i,j) with i<=j<N
static inline constexpr unsigned
count_upper_incl_diag(unsigned N) noexcept
{
  return N * (N + 1) / 2;
}

/**
 * Count of (i,j) pairs with 0 <= i < j < N
 * = N*(N-1)/2
 * N = _op_num
 */
static inline constexpr std::size_t
count_upper(std::size_t N) noexcept
{
  return N * (N - 1) / 2;
}

/**
 * Pack (i,j) with 0 <= i < j < N into a single index in the
 * strict upper-triangular storage of size count_upper(N).
 * (N = _op_num)
 *
 * Layout is column-by-column (by j):
 *   j=1: (0,1)
 *   j=2: (0,2) (1,2)
 *   j=3: (0,3) (1,3) (2,3)
 *   ...
 *
 * Index = base_of_column(j) + i
 * base_of_column(j) = j*(j-1)/2
 *
 * Example:
 *   const auto T = GBPairPacking::count_upper(N);
 *   std::vector<RankTwoTensor> H(T);
 *   H[GBPairPacking::pack_upper(i, j)] = Bij; // with i<j
 */
static inline constexpr std::size_t
pack_upper(std::size_t i, std::size_t j) noexcept
{
  // mooseAssert(i < j, "pack_upper expects i < j");
  if (!(i < j))
    mooseError("pack_upper: need i < j, got i=", i, " j=", j);
  return j * (j - 1) / 2 + i;
}

// Optional: fully-checked versions (use in non-hot code)
inline std::size_t
pack_upper_checked(std::size_t i, std::size_t j, std::size_t N)
{
  if (!(i < j && j < N))
    mooseError("pack_upper_checked: need 0 <= i < j < N, got i=", i, " j=", j, " N=", N);
  const auto k = j * (j - 1) / 2 + i;
  // ensure k is within the allocated range
  if (k >= count_upper(N))
    mooseError("pack_upper_checked: computed index out of range: k=", k);
  return k;
}

// ---------- UNpackers (k -> (i,j)) ----------

// Strict upper (no diagonal), columns j = 1..N-1, rows i = 0..j-1
// Use: auto [i, j] = GBPairPacking::unpack_upper(k);
inline std::pair<std::size_t, std::size_t>
unpack_upper(std::size_t k)
{
  // find column j s.t. j*(j-1)/2 <= k < (j+1)*j/2
  const double kd = static_cast<double>(k);
  const std::size_t j =
      static_cast<std::size_t>(std::floor((1.0 + std::sqrt(1.0 + 8.0 * kd)) / 2.0));
  const std::size_t base = j * (j - 1) / 2;
  const std::size_t i = k - base; // 0 <= i < j
  return {i, j};
}

// Upper including diagonal, columns j = 0..N-1, rows i = 0..j
inline std::pair<std::size_t, std::size_t>
unpack_upper_incl_diag(std::size_t k)
{
  // find column j s.t. j*(j+1)/2 <= k < (j+1)*(j+2)/2
  const double kd = static_cast<double>(k);
  const std::size_t j =
      static_cast<std::size_t>(std::floor((std::sqrt(8.0 * kd + 1.0) - 1.0) / 2.0));
  const std::size_t base = j * (j + 1) / 2;
  const std::size_t i = k - base; // 0 <= i <= j
  return {i, j};
}

}
