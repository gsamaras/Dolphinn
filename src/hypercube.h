#ifndef HYPERCUBE_H
#define HYPERCUBE_H

#include <vector>
#include "hash.h"

#include <thread>
#include <iterator>
#include <utility>

namespace Dolphinn
{
  template <typename T, typename bitT>
  class Hypercube
  {
    // The 'K' hash-functions that we are going to use. Only the last one will be used to query,
    // but we need all of them to map the query on arrival, first.
    std::vector<StableHashFunction<T>> H;
    // original dimension of points
    const int D;
    // mapped dimension of points (dimension of the Hypercube)
    const int K;
    // Reference of an 1D vector of points, emulating a 2D, with N rows and D columns per row.
    const std::vector<T>& pointset;
    public:
    /** \brief Constructor that creates in parallel a 
      * vector from a stable distribution.
      *
      * Assign 'D' random values from a normal distribution
      * N(0,1/sqrt(D)).
      *
      * @param pointset    - 1D vector of points, emulating a 2D, with N rows and D columns per row.
      * @param N           - number of points
      * @param D           - dimension of points
      * @param K           - dimension of Hypercube (and of the mapped points)
      * @param threads_no  - number of threads to be created. Default value is 'std::thread::hardware_concurrency()'.
      * @param r           - parameter of Stable Distribution. Default value is 4.
   */
    Hypercube(const std::vector<T>& pointset, const int N, const int D, const int K, const int threads_no = std::thread::hardware_concurrency(), const int r = 4/*3 or 8*/)
      : D(D), K(K), pointset(pointset)
    {
      if(threads_no >= K || ((K - 1) % threads_no) != 0)
      {
        std::cout << "Threads number is greater or equal to K (dimension of Hypercube). Or  (threads_no MOD (K - 1)) != 0. Construction aborted..." << std::endl;
        return;
      }
      std::vector<bitT> mapped_pointset(N * K);

      if(threads_no == 1)
      {
        for(int k = 0; k < K - 1; ++k)
        {
          H.emplace_back(D, r);
          //H[k].print_a();
          H[k].hash(pointset, N, D);
          //H[k].print_stats();

          H[k].assign_random_bit(mapped_pointset, k, K);
        }
        H.emplace_back(D, r);
        H[K - 1].hash(pointset, N, D);
        H[K - 1].assign_random_bit_and_fill_hashtable_cube(mapped_pointset, K);

        //H[K - 1].print_hashtable_cube();
      }
      else
      {
        std::vector<std::thread> threads;

        std::vector<std::vector<StableHashFunction<T>>> subvectors(threads_no);
        const int subvector_size = (K - 1)/threads_no;
        //std::cout << "subvector_size = " << subvector_size << std::endl;
        for (int i = 0; i < threads_no; ++i)
          threads.push_back(std::thread(populate_vector_of_hash_functions, std::ref(subvectors[i]), subvector_size, D, r, std::ref(pointset), N, std::ref(mapped_pointset), i * subvector_size, K));

        for (auto& th : threads)
          th.join();

        for(auto& subv: subvectors)
        {
          //for(auto& s: subv)
          //{
            H.insert(
              H.end(),
              std::make_move_iterator(subv.begin()),
              std::make_move_iterator(subv.end())
            );
          //}
        }
        H.emplace_back(D, r);
        H[K - 1].hash(pointset, N, D);
        H[K - 1].assign_random_bit_and_fill_hashtable_cube(mapped_pointset, K);

        //H[K - 1].print_hashtable_cube();
      }
    } 

    /** \brief Populate the vector of hash functions.
      * Helper function for the Constructor in a parallel environment.
      *
      * @param H                 - vector of Hash Functions
      * @param n_vec             - nubmer of hash function to be inserted
      * @param D                 - dimension of the original points
      * @param r                 - Stable Distirbution parameter
      * @param pointset          - original points
      * @param N                 - number of origial points
      * @param mapped_pointset   - vector of mapped points (to be poppulated)
      * @param k_start           - starting index of mapped_pointset to be poppulated in parallel
      * @param K                 - dimension of Hypercube
    */
    static void populate_vector_of_hash_functions(std::vector<StableHashFunction<T>>& H, const int n_vec, const int D, const int r, const std::vector<T>& pointset, const int N, std::vector<bitT>& mapped_pointset, const int k_start, const int K)
    {
      for (int i = 0; i < n_vec; ++i)
      {
        H.emplace_back(D, r, k_start + i);
        H[i].hash(pointset, N, D);
        //std::cout << k_start << " " << i << std::endl;
        H[i].assign_random_bit(mapped_pointset, k_start + i, K);
      }
    }

    /** \brief Radius query the Hamming cube.
      *
      * @param query               - vector of queries
      * @param Q                   - number of queries
      * @param radius              - find a point within r with query
      * @param MAX_PNTS_TO_SEARCH  - threshold
      * @param results_idxs        - indices of Q points, where Eucl(point[i], query[i]) <= r
      * @param threads_no          - number of threads to be created. Default value is 'std::thread::hardware_concurrency()'.
    */
    void radius_query(const std::vector<T>& query, const int Q, const int radius, const int MAX_PNTS_TO_SEARCH, std::vector<int>& results_idxs, const int threads_no = std::thread::hardware_concurrency())
    {
      std::vector<bitT> mapped_query(Q * K);
      if(threads_no == 1)
      {
        for(int q = 0; q < Q; ++q)
        {
          for(int k = 0; k < K; ++k)
          {
            H[k].assign_random_bit_query((std::begin(query) + q * D), (std::begin(mapped_query) + q * K), k);
          }
          results_idxs[q] = H[K - 1].radius_query(std::string(mapped_query.begin() + q * K, mapped_query.begin() + (q + 1) * K), radius, K, MAX_PNTS_TO_SEARCH, pointset.begin(), query.begin() + q * D);
        }
      }
      else
      {
        std::vector<std::thread> threads;

        const int batch = Q/threads_no;
        //std::cout << "subvector_size = " << subvector_size << std::endl;
        for (int i = 0; i < threads_no - 1; ++i)
          threads.push_back(std::thread(execute_radius_queries, std::ref(H), std::ref(query), std::ref(mapped_query), i * batch, (i + 1) * batch, K, D, std::ref(pointset), radius, MAX_PNTS_TO_SEARCH, std::ref(results_idxs)));
        threads.push_back(std::thread(execute_radius_queries, std::ref(H), std::ref(query), std::ref(mapped_query), (threads_no - 1) * batch, Q, K, D, std::ref(pointset), radius, MAX_PNTS_TO_SEARCH, std::ref(results_idxs)));
    
        for (auto& th : threads)
          th.join();
      }
      // per query, find candidates
      /*for(int q = 0; q < Q; ++q)
      {
        //std::cout << "Query no. " << q << std::endl;
        results_idxs[q] = H[K - 1].radius_query(std::string(mapped_query.begin() + q * K, mapped_query.begin() + (q + 1) * K), radius, K, MAX_PNTS_TO_SEARCH, pointset.begin(), query.begin() + q * D);
        //std::cout << "Query no. " << q << " completed" << std::endl;
      }*/
    }

    /** \brief Execute specified portion of Radius Queries.
      * Helper function for the Constructor in a parallel environment.
      *
      * @param H                    - vector of Hash Functions
      * @param query                - vector of all queries
      * @param mapped query         - vector of all (to be) mapped queries
      * @param q_start              - starting index of query to execute
      * @param q_end                - ending index of query to execute
      * @param K                    - dimension of Hypercube
      * @param D                    - dimension of original points and queries
      * @param pointset             - original points
      * @param radius               - radius to query with
      * @param MAX_PNTS_TO_SEARCH   - threshold when searching
      * @param results_idxs         - The index of the point-answer in i-th posistion, for i-th query, -1 if not found.
    */
    static void execute_radius_queries(std::vector<StableHashFunction<T>>& H, const std::vector<T>& query, std::vector<bitT>& mapped_query, const int q_start, const int q_end, const int K, const int D, const std::vector<T>& pointset, const int radius, const int MAX_PNTS_TO_SEARCH, std::vector<int>& results_idxs)
    {
      for(int q = q_start; q < q_end; ++q)
      {
        for(int k = 0; k < K; ++k)
        {
          H[k].assign_random_bit_query((std::begin(query) + q * D), (std::begin(mapped_query) + q * K), k);
        }
        results_idxs[q] = H[K - 1].radius_query(std::string(mapped_query.begin() + q * K, mapped_query.begin() + (q + 1) * K), radius, K, MAX_PNTS_TO_SEARCH, pointset.begin(), query.begin() + q * D);
      }
    }

    /** \brief Radius query the Hamming cube.
      *
      * @param query               - vector of queries
      * @param Q                   - number of queries
      * @param MAX_PNTS_TO_SEARCH  - threshold
      * @param results_idxs_dists  - indices and distances of Q points, where the (Approximate) Nearest Neighbors are stored.
      * @param threads_no          - number of threads to be created. Default value is 'std::thread::hardware_concurrency()'.
    */
    void nearest_neighbor_query(const std::vector<T>& query, const int Q, const int MAX_PNTS_TO_SEARCH, std::vector<std::pair<int, float>>& results_idxs_dists, const int threads_no = std::thread::hardware_concurrency())
    {
      std::vector<bitT> mapped_query(Q * K);
      if(threads_no == 1)
      {
        for(int q = 0; q < Q; ++q)
        {
          for(int k = 0; k < K; ++k)
          {
            H[k].assign_random_bit_query((std::begin(query) + q * D), (std::begin(mapped_query) + q * K), k);
          }
          results_idxs_dists[q] = H[K - 1].nearest_neighbor_query(std::string(mapped_query.begin() + q * K, mapped_query.begin() + (q + 1) * K), K, MAX_PNTS_TO_SEARCH, pointset.begin(), query.begin() + q * D);
        }
      }
      else
      {
        std::vector<std::thread> threads;

        const int batch = Q/threads_no;
        //std::cout << "subvector_size = " << subvector_size << std::endl;
        for (int i = 0; i < threads_no - 1; ++i)
          threads.push_back(std::thread(execute_nearest_neighbor_queries, std::ref(H), std::ref(query), std::ref(mapped_query), i * batch, (i + 1) * batch, K, D, std::ref(pointset), MAX_PNTS_TO_SEARCH, std::ref(results_idxs_dists)));
        threads.push_back(std::thread(execute_nearest_neighbor_queries, std::ref(H), std::ref(query), std::ref(mapped_query), (threads_no - 1) * batch, Q, K, D, std::ref(pointset), MAX_PNTS_TO_SEARCH, std::ref(results_idxs_dists)));
    
        for (auto& th : threads)
          th.join();
      }
    }

    /** \brief Execute specified portion of Nearest Neighbor Queries.
      * Helper function for the Constructor in a parallel environment.
      *
      * @param H                    - vector of Hash Functions
      * @param query                - vector of all queries
      * @param mapped query         - vector of all (to be) mapped queries
      * @param q_start              - starting index of query to execute
      * @param q_end                - ending index of query to execute
      * @param K                    - dimension of Hypercube
      * @param D                    - dimension of original points and queries
      * @param pointset             - original points
      * @param MAX_PNTS_TO_SEARCH   - threshold when searching
      * @param results_idxs_dists  - indices and distances of Q points, where the (Approximate) Nearest Neighbors are stored.
    */
    static void execute_nearest_neighbor_queries(std::vector<StableHashFunction<T>>& H, const std::vector<T>& query, std::vector<bitT>& mapped_query, const int q_start, const int q_end, const int K, const int D, const std::vector<T>& pointset, const int MAX_PNTS_TO_SEARCH, std::vector<std::pair<int, float>>& results_idxs_dists)
    {
      for(int q = q_start; q < q_end; ++q)
      {
        for(int k = 0; k < K; ++k)
        {
          H[k].assign_random_bit_query((std::begin(query) + q * D), (std::begin(mapped_query) + q * K), k);
        }
        results_idxs_dists[q] = H[K - 1].nearest_neighbor_query(std::string(mapped_query.begin() + q * K, mapped_query.begin() + (q + 1) * K), K, MAX_PNTS_TO_SEARCH, pointset.begin(), query.begin() + q * D);
      }
    }

    /** \brief Print how many points are assigned to every vertex.
      * Empty vertices (if any) are not printed (because we do not store them).
      *
    */
    void print_no_of_assigned_points_per_vertex()
    {
      H[K - 1].print_hashtable_cube();
    }

  };
}

#endif /* HYPERCUBE_H */
