[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/HlQKP7Zu)

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
Many methods include **parallel execution policies**, to accelerate certain procedure, as maximum search or vector filling.

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

Queste ultime in particolare possono essere selezionate sfruttando l'apposito file `data/data.json`, che contiene una lista dei nomi delle che verranno testate dal codice. Per testare matrici diverse basterà quindi, dopo aver eventualmente aggiunto la matrice desiderata all'interno della cartella `data`, modificare questa lista con i nomi di interesse.

**Importante**: l'informazione appena citata potrebbe essere utile nello scenario in cui si voglia evitare di svolgere il test sulla seconda matrice, che è di dimensioni decisamente maggiori e potrebbe richiedere un maggiore tempo di esecuzione del codice.

I test che sono stati effettuati sono i seguenti:
- validità dell'operazione di **compressione**, corrispondente ai metodi `compress()`, `uncompress()` e `compress_mod()` (quest'ultimo disponibile solo per la classe _SquareMatrix_)
- calcolo delle tre tipologie di **norme**, sfruttando il metodo `norm<N>()`, where $N$ is type of the norm (One, Infinity, Frobenius)
- tempo di esecuzione del **prodotto matrice-vettore**, dove il vettore è stato generato in modo casuale
- tempo di esecuzione **prodotto matrice-matrice**, sfruttando il fatto che le matrici teste sono quadrate

**Speedup**: the execution times of the products are tested using the matrix in both compressed and uncompressed formats amd they times are saved in `execution_time.json`, allowing them to be displayed on the screen along with the corresponding speedups, calculated as
$$\text{speedup} = \frac{\text{execution time in uncompressed format}}{\text{execution time in compressed format}},$$
so that we can appreciate the improvements in terms of speed achieved thanks to the compressed format.


