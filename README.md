[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/HlQKP7Zu)

# Sparse matrix 

Test:
- reader() for matrix market format
- with a given matrix
- time with Chrono.hpp

Parallelization with std algorithms 
(or blas library: ask other teams and check pacs-examples)

Modified compressed techniques (MSR, MSC)
- abstract class AbstractMatrix
- 2 classes: Matrix, SquareMatrix
- Matrix has classical compressed formats, while SquareMatrix class has modified compressed formats
- Both will have transpose and diagonal views
  
README.md (need of TBB, we didn't implement COO since COOmap is clearly more efficient)

