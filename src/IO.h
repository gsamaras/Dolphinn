#ifndef IO_H
#define IO_H

#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <utility>

/** \brief Read a collection of points from file.
 *
 * Dimension and number of points should have been
 * assigned a value before reaching this function.
 *
 * @param v        - vector of points
 * @param N        - number of points
 * @param D        - dimension of points
 * @param filename - input file
 */
template<typename T>
void read_points(std::vector<T>& v, int N, int D, const char* filename)
{
  std::ifstream infile;

  infile.open(filename);
  if (!infile)
    std::cout << "File not found!" << std::endl;

  int hRead = 0;
  for (int i = 0; i < N && infile; ++i) {
    for (int j = 0; j < D; ++j)
      infile >> v[i * D + j];
    hRead++;
  }
  if (hRead != N)
    std::cout << "ERROR, read less than " << N << " points!!\n\n";
}

/** \brief Read a collection of points from file. Every line of the file has this format: D x_1 ... x_D
 *
 * Dimension and number of points should have been
 * assigned a value before reaching this function.
 *
 * @param v        - vector of points
 * @param N        - number of points
 * @param D        - dimension of points
 * @param filename - input file
 */
template<typename T>
void readfvecs(std::vector<T>& v, int N, int D, const char* filename) {
  if(!std::is_same<T, float>::value && !std::is_same<T, double>::value)
    std::cout << "WARNING!!!!!!!!! You are reasding fvecs, but T is neither float nor double!" << std::endl;
  FILE* fid;
  fid = fopen(filename, "rb");

  if (!fid) {
    printf("I/O error : Unable to open the file %s\n", filename);
    //std::cerr << "Error: " << strerror(errno) << std::endl;
  }

  // we assign the return value of fread() to 'sz' just to suppress a warning
  int foundD;
  size_t sz = fread(&foundD, sizeof(foundD), 1, fid);
  fseek(fid, 0L, SEEK_END);
  sz = ftell(fid);
  int foundN = sz / (1 * 4 + D * 4);
  printf("foundN = %d, N = %d, foundD = %d, D = %d, |%s|\n", foundN, N, foundD, D, filename);
  if(N != foundN || D != foundD)
    std::cout << "WARNING, see above!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

  fseek(fid, 0L, SEEK_SET);
  //std::cout << ds.dim() << " " << ds.size() << "\n";
  int c = 0;
  float value;
  int i, j;
  for(i = 0; i < N; ++i) {
    sz = fread(&D, sizeof(D), 1, fid);
    //printf("%d\n", D);
    for (j = 0; j < D; ++j) {
      sz = fread(&value, sizeof(value), 1, fid);
      //if(c >= 279619)
      //printf("j = %d, value = %f, read up to point %d\n", j, value, c);
      v[i * D + j] = value;
    }
    ++c;
    //printf("read up to %d\n", c);
  }
  if(c != N)
    printf("WARNING! Read less points than expected.\n");
}

/** \brief Helper function to read a file in IDX format.
 *
 * @param i - integer to reversed
 * @return  - reversed integer
 */
int reverseInt (int i) 
{
    unsigned char c1, c2, c3, c4;

    c1 = i & 255;
    c2 = (i >> 8) & 255;
    c3 = (i >> 16) & 255;
    c4 = (i >> 24) & 255;

    return ((int)c1 << 24) + ((int)c2 << 16) + ((int)c3 << 8) + c4;
}

/** \brief Read a file in IDX format.
 *
 * Dimension and number of points should have been
 * assigned a value before reaching this function.
 *
 * @param v        - vector of points
 * @param N        - number of points
 * @param D        - dimension of points
 * @param filename - input file
 */
template<typename T>
void read_points_IDX_format(std::vector<T>& v, int N, int D, const char* filename)
{
  if(!std::is_same<T, int>::value)
    std::cout << "WARNING!!!!!!!!! You are reasding IDX format (usually MNIST), but T is not int!" << std::endl;
  std::ifstream file (filename);

  int magic_number=0;
  int number_of_images=0;
  int n_rows=0;
  int n_cols=0;
  int coord_idx = 0;
  if (file.is_open())
  {
    file.read((char*)&magic_number,sizeof(magic_number)); 
    magic_number= reverseInt(magic_number);
    file.read((char*)&number_of_images,sizeof(number_of_images));
    number_of_images= reverseInt(number_of_images);
    file.read((char*)&n_rows,sizeof(n_rows));
    n_rows= reverseInt(n_rows);
    file.read((char*)&n_cols,sizeof(n_cols));
    n_cols= reverseInt(n_cols);
    //std::cout << "Read IDX formatted file, with n_rows = " << n_rows <<" and n_cols = " << n_cols << std::endl;
    for(int i = 0; i < number_of_images && i < N; ++i)
    {
      for(int r=0;r<n_rows;++r)
      {
        for(int c=0;c<n_cols;++c)
        {
          unsigned char temp = -1;
          file.read((char*)&temp,sizeof(temp));
          //std::cout << (int)temp << std::endl;
          v[coord_idx++] = (int)temp;
        }
      }
    }
  }
  std::cout << "Read IDX formatted file, with number_of_images = " << number_of_images << " and dimension (n_rows x n_cols) = " << (n_rows * n_cols) << std::endl;
  if (number_of_images != N)
    std::cout << "ERROR, read less than " << N << " points!!\n\n";
  if (n_rows * n_cols != D)
    std::cout << "ERROR, dimension less than " << D << " points!!\n\n";
}

/** \brief Read a custom format of Crow features,
 * based on the Oxford dataset.
 *
 * Dimension and number of points should have been
 * assigned a value before reaching this function.
 *
 * @param data     - vector of points
 * @param N        - number of points. Usually 5063.
 * @param D        - dimension of points. Usually 512.
 * @param filename - input file
 */
template<typename T>
void read_crow_features_oxford(std::vector<T>& data, int N, int D, const char* filename) {
  if(!std::is_same<T, float>::value && !std::is_same<T, double>::value)
    std::cout << "WARNING!!!!!!!!! You are reasding CROW features for Oxford dataset, but T is neither float nor double!" << std::endl;

    std::string line;
    std::ifstream infile(filename);
    int i = 0, j = 0;
    while(std::getline(infile, line)) {
        bool exists = line.find("[") != std::string::npos;
        if(exists) {
            std::istringstream iss(line);
            char bracket;
            iss >> bracket;
            float v;
            while (iss >> v) {
                data[i * D + j++] = v;
            }
            continue;
        }
        exists = line.find("]") != std::string::npos;
        std::istringstream iss(line);
        float v;
        while (iss >> v) {
            data[i * D + j++] = v;
        }
        if(exists) {
            j = 0;
            ++i;
        }
        //break;
    }
    if(i != N)
        std::cout << "WARNING!!!!!!!! READ less/more points (" << i << ") than N =" << N << std::endl;
    infile.close();
}

/** \brief Read a custom format of Crow features,
 * based on the Oxford dataset's queries.
 *
 * Dimension and number of points should have been
 * assigned a value before reaching this function.
 *
 * @param query     - vector of points
 * @param N        - number of points. Usually 55.
 * @param D        - dimension of points. Usually 512.
 * @param filename - input file
 */
template<typename T>
void read_crow_features_oxford_queries(std::vector<T>& query, int Q, int D, const char* filename) {
    std::ifstream query_file;
    float v;
    query_file.open (filename);
    for(int i = 0; i < Q; ++i) {
        for(int j = 0; j < D; ++j) {
            query_file >> v;
            query[i * D + j] = v;
        }
    }
    query_file.close();
}

/** \brief Print 2D vector.
 *
 * @param v           - 1D vector of points
 * @param N           - number of points
 * @param D           - dimension of points
 */
template<typename T>
void print_2D_vector(const std::vector<T>& v, const int N, const int D)
{
  for(unsigned int i = 0; i < N; ++i)
  {
  	for(unsigned int j = 0; j < D; ++j)
  	{
  		(std::is_same<T, char>::value) ? (std::cout << (int)v[i * D + j] << " ") : (std::cout << v[i * D + j] << " ");
  	}
  	std::cout << "\n";
  }
}

/** \brief Print 1D vector.
 *
 * @param v           - vector to be printed
 */
template<typename T>
void print_1D_vector(const std::vector<T>& v)
{
  for (auto& v_i: v)
  {
    std::cout << v_i << " ";
  }
  std::cout << "\n";
}

/** \brief Print 1D vector of pairs.
 *
 * @param v           - vector to be printed
 */
void print_1D_vector(const std::vector<std::pair<int, float>>& v)
{
  for (auto& v_i: v)
  {
    std::cout << "(" << v_i.first << ", " << v_i.second << ") ";
  }
  std::cout << "\n";
}

/** \brief Print std::string with every element casted to int.
 *
 * @param str - string
 */
void print_string_cast_int(const std::string& str)
{
  for(auto& character: str)
    std::cout << (int)character;
}

#endif /*IO_H*/