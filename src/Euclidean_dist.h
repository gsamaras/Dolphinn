#ifndef EUCLIDEAN_DIST_H
#define EUCLIDEAN_DIST_H

#include <vector>

/** \brief Euclidean distance squared.
 *
 * @param it1       - first point
 * @param it1_end   - end of first point
 * @param it2       - second point
 * @return          - the Euclidean distance of p1-p2
 */
template<typename iterator>
float squared_Eucl_distance(iterator it1, iterator it1_end, iterator it2)
{
  float squared_distance = 0.;
  float diff;
  for (; it1 < it1_end; ++it1, ++it2)
  {
    diff = *it1 - *it2;
    squared_distance += diff * diff;
  }
  return squared_distance;
}

/** \brief Report a point's index (if any) that has Euclidean distance
 * less or equal than a given radius.
 *
 * @param pointset        - 1D vector of all points
 * @param points_idxs     - indices of candidate points
 * @param D               - dimension of points
 * @param query_point     - vector containing only the coordinates of the query point
 * @param squared_radius  - square value of given radius
 * @param threshold       - max number of points to check
 * @return                - the index of the point. -1 if not found.
 */
template <typename iterator>
int Euclidean_distance_within_radius(iterator pointset, const std::vector<int>& points_idxs,
 const int D, iterator query_point, const int squared_radius, const int threshold)
{
  int size = points_idxs.size();
  for(unsigned int i = 0; i < threshold && i < size; ++i)
  {
    if(squared_Eucl_distance(query_point, query_point + D, pointset + points_idxs[i] * D) <= squared_radius)
      return points_idxs[i];
  }
  return -1;
}


#endif /*EUCLIDEAN_DIST_H*/