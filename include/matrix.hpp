#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "storage.hpp"
#include <vector>

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
        // default constructor
        Matrix();
        // constructor with size
        Matrix(size_t rows, size_t cols);

        // set an element in the matrix (dynamic construction of the matrix)
        void set(size_t row, size_t col, const T &value);

        // compress
        void compress();
        // uncompress
        void uncompress();
        // check if the matrix is compressed
        bool is_compressed() const;

        // call operator() const version
        T operator()(size_t i, size_t j) const;
        // call operator() non-const version
        T &operator()(size_t i, size_t j);

        // resize the matrix
        void resize(size_t rows, size_t cols);

        // calculate the norm of the matrix
        template <NormType N>
        T norm() const;

        // friend functions
        // multiply with a std::vector
        friend std::vector<T> operator*(const Matrix &m, const std::vector<T> &v);
        // multiply with another matrix
        friend Matrix operator*(const Matrix &m1, const Matrix &m2);

    private:
        size_t rows;     // number of rows
        size_t cols;     // number of columns
        bool compressed; // flag to check if the matrix is compressed

        // storage for the matrix
        // uncompressed matrix
        UncompressedStorage<T, S> uncompressed; // COO format
        // compressed matrix
        CompressedStorage<T> compressed; // CSR or CSC format
    };

    template <AddMulType T, StorageOrder S>
    std::vector<T> operator*(const Matrix<T, S> &m, const std::vector<T> &v)
    {
        // check if the matrix is compressed
        if (m.is_compressed())
        {
            // uncompress the matrix
            m.uncompress();
        }
        // multiply the matrix with the vector
        std::vector<T> result(m.rows());
        for (size_t i = 0; i < m.rows(); i++)
        {
            result[i] = 0;
            for (size_t j = 0; j < m.cols(); j++)
            {
                result[i] += m(i, j) * v[j];
            }
        }
        return result;
    }

    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const Matrix<T, S> &m1, const Matrix<T, S> &m2)
    {
        // TBD
    }
}

#endif MATRIX_HPP
