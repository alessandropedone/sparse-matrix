# Sparse matrix 

# Matrix class

## Set up
Clone the repository with the command: 
```bash
git clone --recurse-submodules git@github.com:PACS-24-25/challenge2-male.git
```
Be aware of the fact that _TBB_ library is required to compile and execute the code.

## Implementation
We implemented all the methods requested and we built the following class-structure:
- **AbstractMatrix**: base abstract template class that allows greater polymorphism among derived classes
    - **Matrix**: template class that encodes **COO Map** and **CSR/CSC** formats, as data structures for a matrix
        - **SquareMatrix**: template class specialized for square matrices, that encodes also the **MSR/MSC** format
    - Views
        - **TransposeView**
        - **DiagonalView**

**Note**: in `abstract_matrix.hpp` you can see how we have conceptualized a matrix.

### Design choises
1) Data formats and associated components are defined in `Storage.hpp`.
2) We used concepts to require some specific properties (only the necessary ones) about the type of the elements of matrices, as illustrated by the following piece of code.
    ```cpp
    template <typename T>
    concept AddMulType = requires(T a, T b) {
        { a + b } -> std::convertible_to<T>;
        { a * b } -> std::convertible_to<T>;
        { std::abs(a) } -> std::convertible_to<AbsReturnType_t<T>>;
    };
    ```
3) For the dynamic storage techniques, among COO format and COOmap, we opted for the latter, since it provides access to random elements with $O(log(N))$ average complexity and it yields an easy and fast way to insert the elements in order.
4) `Proxy.hpp` was conceived to provide restricted access to private data, while still allowing operations such as:
    ```cpp
    m(0, 0) = m(1, 1) + m(2, 2);
    ```
    This approach ensures encapsulation while enabling controlled manipulation of matrix elements.

### Products
Below there are the declarations of all the matrix products we have implemented as friend functions of _Matrix_ and _SquareMatrix_ classes.
```cpp
// Matrix products
template <AddMulType U, StorageOrder V>
std::vector<U> operator*(const Matrix<U, V> &m, const std::vector<U> &v);

template <AddMulType U, StorageOrder V>
Matrix<U, V> operator*(const Matrix<U, V> &m1, const Matrix<U, V> &m2);

// SquareMatrix products
template <AddMulType U, StorageOrder V>
std::vector<U> operator*(const SquareMatrix<U, V> &m, const std::vector<U> &v);

template <AddMulType U, StorageOrder V>
SquareMatrix<U, V> operator*(const SquareMatrix<U, V> &m1, const SquareMatrix<U, V> &m2);

// TransposeView products 
template <AddMulType U, StorageOrder V>
std::vector<U> operator*(const TransposeView<U, V> &m, const std::vector<U> &v);

template <AddMulType U, StorageOrder V>
Matrix<U, V> operator*(const TransposeView<U, V> &m1, const TransposeView<U, V> &m2);

// DiagonalView products
template <AddMulType U, StorageOrder V>
std::vector<U> operator*(const DiagonalView<U, V> &m, const std::vector<U> &v);

template <AddMulType U, StorageOrder V>
SquareMatrix<U, V> operator*(const DiagonalView<U, V> &m1, const DiagonalView<U, V> &m2);

template <AddMulType U, StorageOrder V>
Matrix<U, V> operator*(const Matrix<U, V> &m1, const DiagonalView<U, V> &m2);

template <AddMulType U, StorageOrder V>
Matrix<U, V> operator*(const DiagonalView<U, V> &m1, const Matrix<U, V> &m2);
```
### Parallelization 
Many methods include **parallel execution policies**, to accelerate certain procedures, as maximum search or vector filling.\
We retained the method `compress_parallel()`, available only for the _Matrix_ class, designed to perform the transition from the uncompressed format to the compressed format using a parallel approach with `std::atomic`. However, we did not further develop this idea because the overhead caused by creating an index vector is too significant, primarily due to the use of the `std::iota` function.

## Test
The code has been tested in two different ways, always starting from matrices in _Matrix Market format_.

1) Initially, with matrices of size $(5 \times 5)$:
    - `real_test_5x5.mtx`, containing real entries
    - `complex_test_5x5.mtx`, containing complex entries  
    
    These matrices allow for an easy visualization of the results in the terminal.

2) Subsequently, with larger matrices:
    - [`lnsp_131.mtx`](https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/lns/lnsp_131.html) $(131 \times 131)$
    - [`e20r0000.mtx`](https://math.nist.gov/MatrixMarket/data/SPARSKIT/drivcav/e20r0000.html) $(1182 \times 1182)$

These matrices, in particular, can be selected using the dedicated file `data/data.json`, which contains a list of the names of the matrices that will be tested by the code. To test different matrices, it will be sufficient to modify this list with the desired names after optionally adding the desired matrix to the `data` folder.

**Important**: the aforementioned information could be useful in scenarios where one wishes to avoid testing the second matrix, which has significantly larger dimensions and might require more execution time for the code.

The tests that have been performed are as follows:
- validation of the **compression** operation, corresponding to the methods `compress()`, `uncompress()`, and `compress_mod()` (the latter available only for the _SquareMatrix_ class)
- calculation of the three types of **norms**, using the method `norm<N>()`, where $N$ is the type of the norm (One, Infinity, Frobenius)
- execution time of the **matrix-vector product**, where the vector was randomly generated
- execution time of the **matrix-matrix product**, leveraging the fact that the tested matrices are square

**Speedup**: the execution times of the products are taken (exploiting `std::chrono`) for the matrix in both compressed and uncompressed formats and they are saved in `execution_time.json`, allowing them to be displayed on the screen along with the corresponding speedups, calculated as
$$\text{speedup} = \frac{\text{execution time in uncompressed format}}{\text{execution time in compressed format}},$$
so that we can appreciate the improvements in terms of speed achieved thanks to the compressed format.


