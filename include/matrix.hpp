#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "storage.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <execution>
#include <cmath>

namespace algebra
{

    
    /// @brief type of norm
    enum class NormType
    {
        One,
        Infinity,
        Frobenius
    };

    /// @brief Matrix class
    /// @tparam T type of the matrix elements
    /// @tparam S storage order of the matrix (RowMajor or ColumnMajor)
    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class Matrix
    {
    public:
        // delete default constructor
        Matrix() = delete;

        /// @brief constructor with size
        /// @param rows number of rows
        /// @param cols number of columns
        Matrix(size_t rows, size_t cols)
        {
            this->rows = rows;
            this->cols = cols;
            this->compressed = false;
        };

        /// @brief set an element in the matrix (dynamic construction of the matrix)
        /// @param row row index
        /// @param col column index
        /// @param value value to set
        void set(size_t row, size_t col, const T &value);

        /// @brief check if the matrix is in a compressed format
        /// @return true if the matrix is compressed, false otherwise
        bool is_compressed() const { return compressed; };

        /// @brief compress the matrix if it is in an uncompressed format
        void compress();

        /// @brief uncompress the matrix if it is in a compressed format
        void uncompress();

        /// @brief call operator() const version
        /// @param row row index
        /// @param col column index
        /// @return element at (row, col)
        T operator()(size_t row, size_t col) const;

        /// @brief call operator() non-const version
        /// @param row row index
        /// @param col column index
        /// @return reference to the element at (row, col)
        T &operator()(size_t row, size_t col);

        /// @brief resize the matrix
        /// @param rows number of rows
        /// @param cols number of columns
        void resize_and_clear(size_t rows, size_t cols);

        /// @brief calculate the norm of the matrix
        /// @tparam N type of the norm (One, Infinity, Frobenius)
        /// @return value of the norm
        template <NormType N>
        double norm() const;

        /// @brief Function to read a matrix in Matrix Market format
        /// @param filename input file name
        void reader(const std::string &filename);

        // friend functions
        // multiply with a std::vector
        template <AddMulType U, StorageOrder V>
        friend std::vector<U> operator*(const Matrix<U, V> &m, const std::vector<U> &v);
        // multiply with another matrix
        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const Matrix<U, V> &m1, const Matrix<U, V> &m2);

    protected:
        size_t rows;             // number of rows
        size_t cols;             // number of columns
        bool compressed = false; // flag to check if the matrix is compressed

        // storage for the matrix
        // uncompressed matrix
        UncompressedStorage<T, S> uncompressed_format; // COO format
        // compressed matrix
        CompressedStorage<T> compressed_format; // CSR or CSC format
    };

}

#include "matrix.tpp"

#endif // MATRIX_HPP
