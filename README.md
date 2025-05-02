[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/HlQKP7Zu)

# Sparse matrix 

Views:
- implement all products
- test transpose and diagonal products
- move semantic for 

```cpp
template <AddMulType T, StorageOrder S>
Matrix<T, S> operator*(const TransposeView<T, S> &m1, const TransposeView<T, S> &m2)
```

Optimization
- optimize compression in modified format (remove unnecessary ifs and avoid push_back)
- review normal algorithms and add parallel execution policies where possibile
  
README.md (
    need of TBB, 
    we didn't implement COO since COOmap is clearly more efficient, 
    explain why parallel implementation of compress is slower, put a test about this [keep compress_parallel in order to see the parallelization overhead]
    please run multiple times,
    results interpretation,
    structure of the code and the output)

//Final version of readme

# Matrix class

## Set up

Clone the repository with the command: 
```bash
git clone --recurse-submodules git@github.com:PACS-24-25/challenge2-male.git
```

Be aware of the fact that _TBB library_ is required to compile and execute the code.

## Implementation
We implemented the following class-structure:
- **Abstract matrix**: base abstract template class that allows greater polymorphism among derived classes
    - **Matrix**: template class than encodes two data structures for a matrix, in particular a sparse matrix → **COO Map** and **CSR/CSC** formats
        - **SquareMatrix**: template class specialized for square matrices, that encodes the **MSR/MSC** format
    - **Views** → **TransposeView** and **Diagonal view**: template classes that stand for specific views of a matrix

### Main features section?

## Test case
The code has been checked out on three different matrices:
- 5x5 Matrix: 'read_test_5x5.mtx', user-defined
- 131x131 Matrix: ['lnsp_131.mtx'](https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/lns/lnsp_131.html)
- 1182x1182 Matrix: ['e20r0000.mtx'](https://math.nist.gov/MatrixMarket/data/SPARSKIT/drivcav/e20r0000.html)

For each of those, the following processes were tested, for both storage orders:
- compress(), uncompress(), compress_mod() (only for SquareMatrix)
--------- aggiungiamo compression_test ad execute_test?
- norm<<N>>()
- multiplication with:
    - std::vector<<T>>
    - Matrix<<T, S>>
    - TransposeView<<T, S>>
    - DiagonalView<<T, S>>
    Note:
        - times of execution are obtained through std::chrono functions and printed on screen
        - operations are implemented for matrices of the same format
        - each format is examined and they are compared, as the compressed format determines a speedup
        - the execution data is saved in 'execution_time.json'
    

