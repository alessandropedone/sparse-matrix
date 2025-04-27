[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/HlQKP7Zu)

# Sparse matrix 

Modified compressed techniques (MSR, MSC)
- add all methods to views
- SquareMatrix has a transpose and diagonal view

how do i manage to compute products (views and SquareMatrix)
- for diagonalview, what matrix-matrix product shall we implement? (diag-diag, diag-matrix, matrix-diag)

test also the product with modified compressed format for SquareMatrix

Parallelization with std algorithms
- keep compress_parallel in order to see the parallelization overhead
- review normal algorithms and add parallel execution policies where possibile
  
README.md (
    need of TBB, 
    we didn't implement COO since COOmap is clearly more efficient, 
    explain why parallel implementation of compress is slower, put a test about this)