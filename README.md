# DOLPHINN
A C++ library for: Dimension reductiOn and LookuPs on a Hypercube for effIcient Near Neighbor.

Part of the Data Science Master Thesis of George Samaras, National Kapodistrian University of Athens, 2016.

src/main.cpp contains a representative example. In general, a data-structure of class Hypercube should be constructed and then one can execute Near (radius) or Nearest Neighbor queries on that Hypercube.

Note: If you are interested in Nearest Neighbor, use [DolphinnPy](https://github.com/ipsarros/DolphinnPy).

---

DOLPHINN provides with a simple, yet efficient method for the problem of computing an (approximate) nearest neighbor in high dimensions. The algorithm is based on our paper: [Practical linear-space Approximate Near Neighbors in high dimension](https://arxiv.org/pdf/1612.07405.pdf)[Avarikioti, Prof. Emiris, Psarros (original idea) and Samaras], where we show linear space and sublinear query for a specific setting of parameters.

First, N points are randomly mapped to keys in {0,1}^K, for K<=logN, by making use of the Hypeplane LSH family. Then, for a given query, candidate nearest neighbors are the ones within a small hamming radius with respect to their keys. Our approach resembles the multi-probe LSH approach but it differs on how the list of candidates is computed.
