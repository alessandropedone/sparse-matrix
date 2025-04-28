[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/HlQKP7Zu)

# Sparse matrix 

Optimize compression in modified format (remove unnecessary ifs)
Test
- product with modified compressed format for SquareMatrix
- transpose and diagonal views

how do i manage to compute products (views and SquareMatrix)
- implement 3 cases for each view
- see use of private member of matrix in the view
    - possibility of public getter of const & to data members of matrix
- discuss about return types
    - create matrix constructor with transposed view
    - create squareMatrix constructor with transposed view and diagonal view

Parallelization with std algorithms
- keep compress_parallel in order to see the parallelization overhead
- review normal algorithms and add parallel execution policies where possibile
  
README.md (
    need of TBB, 
    we didn't implement COO since COOmap is clearly more efficient, 
    explain why parallel implementation of compress is slower, put a test about this)