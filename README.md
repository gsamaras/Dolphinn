# DOLPHINN
*DOLPHINN* is a **C++ [header-only](http://en.wikipedia.org/wiki/Header-only) library** for: Dimension reductiOn and LookuPs on a Hypercube for effIcient Near Neighbor.


## How to use DOLPHINN?

Just include DOLPHINN's header file. src/main.cpp contains a representative example.

Note: If you are interested in Nearest Neighbor, use [DolphinnPy](https://github.com/ipsarros/DolphinnPy).

## How fast is DOLPHINN?

On data sets with more than 1 million points in around 128 dimensions, DOLPHINN typically requires only some milliseconds per query.

---

DOLPHINN provides with a simple, yet efficient method for the problem of computing an (approximate) nearest neighbor in high dimensions. The algorithm is based on our paper: [Practical linear-space Approximate Near Neighbors in high dimension](https://arxiv.org/pdf/1612.07405.pdf)[Avarikioti, Prof. Emiris, Psarros (original idea) and Samaras], where we show linear space and sublinear query for a specific setting of parameters. Part of the Data Science Master Thesis of George Samaras, National Kapodistrian University of Athens, 2016.

First, N points are randomly mapped to keys in {0,1}^K, for K<=logN, by making use of the Hypeplane LSH family. Then, for a given query, candidate nearest neighbors are the ones within a small hamming radius with respect to their keys. Our approach resembles the multi-probe LSH approach but it differs on how the list of candidates is computed.
