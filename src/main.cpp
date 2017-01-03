#include <iostream>
#include <vector>

#include "IO.h"
#include "memory.h"
#include "hypercube.h"
//#include "hashCompareFalconn.h"

#include <ctime>
#include <ratio>
#include <chrono>

#define N 60000
#define D 784
#define K 3 //floor(log2(N)/2)
#define Q 10000
#define T int
#define bitT char
#define MAX_PNTS_TO_SEARCH 60000 //N * 1 / 100
#define RADIUS 1

//     /usr/bin/time -l ./fish
//     ssh gsamaras@195.134.71.24
int main()
{
	if(MAX_PNTS_TO_SEARCH > N)
	{
    	std::cerr << "MAX_PNTS_TO_SEARCH > N" << std::endl;
        return -1;
  	}

  	// vector is actually 1D, emulating a 2D vector
  	std::vector<T> pointset(N * D);

  	const std::string data_filename = "/home/gsamaras/Desktop/Data/Sphere/dim512/Sphere_nb_P" + std::to_string(N) + "_dim" + std::to_string(D) + ".txt";
  	const std::string query_filename = "/home/gsamaras/Desktop/Data/Sphere/dim512/queries_Sphere_nb_P" + std::to_string(Q)  + "_dim" + std::to_string(D) + ".txt";

  	std::cout << "N = " << N << ", D = " << D << ", K = " << K << ", MAX_PNTS_TO_SEARCH = " << MAX_PNTS_TO_SEARCH << std::endl;

  	//readfvecs<T>(pointset, N, D, "../create_pointset/siftsmall/siftsmall_base.fvecs");
  	//read_points<T>(pointset, N, D, data_filename.c_str());
  	read_points_IDX_format<T>(pointset, N, D, "/Users/gsamaras/Code/C++/create_pointset/MNIST/train-images-idx3-ubyte");

	//print_2D_vector<T>(pointset, N, D);

  	using namespace std::chrono;
  	high_resolution_clock::time_point t1 = high_resolution_clock::now();
 
  	Dolphinn::Hypercube<T, bitT> hypercube(pointset, N, D, K, 2);

  	high_resolution_clock::time_point t2 = high_resolution_clock::now();
 
  	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
 
  	std::cout << "Build: " << time_span.count() << " seconds.\n";


  	hypercube.print_no_of_assigned_points_per_vertex();
  
  	//print_2D_vector<bitT>(mapped_pointset, N, K, true);


  	// QUERY
  	std::vector<T> query(Q * D);
  	//readfvecs<T>(query, Q, D, "../create_pointset/siftsmall/siftsmall_base.fvecs");
  	//read_points<T>(query, Q, D, query_filename.c_str());
  	read_points_IDX_format<T>(query, Q, D, "/Users/gsamaras/Code/C++/create_pointset/MNIST/t10k-images-idx3-ubyte");

 
  	std::vector<int> results_idxs(Q);

  	t1 = high_resolution_clock::now();

    hypercube.radius_query(query, Q, RADIUS, MAX_PNTS_TO_SEARCH, results_idxs, 2);

  	t2 = high_resolution_clock::now();
  	time_span = duration_cast<duration<double>>(t2 - t1);
 
  	std::cout << "Search: " << time_span.count()/(double)Q << " seconds.\n";

  	t1 = high_resolution_clock::now();

  	int squared_radius = RADIUS * RADIUS;
  	std::vector<int> brute_results_idxs(Q);
  	for(int q = 0; q < Q; ++q)
  	{
    	for(int n = 0; n < N; ++n)
    	{
      		/*if(q == 8 && (n == 8 || n == 0))
      		{
        		int dist = squared_Eucl_distance(query.begin() + q * D, (query.begin() + q * D) + D, pointset.begin() + n * D);

        		std::cout << dist << " " << squared_radius << std::endl;
        		std::cout << (dist <= squared_radius) << std::endl;

      		}*/
      		if(squared_Eucl_distance(query.begin() + q * D, (query.begin() + q * D) + D, pointset.begin() + n * D) <= squared_radius)
      		{
        		brute_results_idxs[q] = q;
        		break;
      		}
      		else
      		{
        		brute_results_idxs[q] = -1;
    		}
    	}
  	}

  	t2 = high_resolution_clock::now();
  	time_span = duration_cast<duration<double>>(t2 - t1);
 
  	std::cout << "Brute force: " << time_span.count()/(double)Q << " seconds.\n";


  	print_1D_vector(results_idxs);
  	print_1D_vector(brute_results_idxs);
  	int correct = 0;
  	for(int q = 0; q < Q; ++q)
    	if( (results_idxs[q] != -1 && brute_results_idxs[q] != -1) || (results_idxs[q] == -1 && brute_results_idxs[q] == -1))
      		correct++;
  	std::cout << "Correct = " << (correct * 100)/(double)Q << "% , correct = " << correct << ", Q = " << Q << std::endl;

  	return 0;
}
