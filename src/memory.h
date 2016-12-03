#ifndef MEMORY_H
#define MEMORY_H

#include <vector>

/** \brief Resize 2D vector.
 *
 * Dimension and number of points should have been
 * assigned a value before reaching this function.
 *
 * @param v        - vector of points
 * @param N        - number of points
 * @param D        - dimension of points
 */
template<typename T>
void resize_2D_vector(std::vector<T>& v, int N, int D)
{
  v.resize(N * D);
}

#endif /*MEMORY_H*/