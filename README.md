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