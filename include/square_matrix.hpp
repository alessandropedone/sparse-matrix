#ifndef SQUARE_MATRIX_HPP
#define SQUARE_MATRIX_HPP

#include "storage.hpp"
#include "matrix.hpp"
#include "proxy.hpp"
#include "matrix_views.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <execution>
#include <cmath>

namespace algebra
{
    /// @brief Square Matrix class, child of Matrix class
    /// @tparam T type of the matrix elements
    /// @tparam S storage order of the matrix (RowMajor or ColumnMajor)
    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class SquareMatrix : public Matrix<T, S>
    {
    public:
        // delete default constructor
        SquareMatrix() = delete;

        /// @brief constructor with size
        /// @param size number of rows and columns
        SquareMatrix(int size) : Matrix<T, S>(size, size)
        {
            this->compressed = false;
            this->modified = false;
        };

        /// @brief default copy constructor
        /// @param other matrix to copy
        SquareMatrix(const SquareMatrix &other) = default;

        /// @brief constructor from a matrix
        /// @param other matrix to copy
        SquareMatrix(const Matrix<T, S> &other) : Matrix<T, S>(other)
        {
            if (other.get_rows() != other.get_cols())
            {
                throw std::runtime_error("Matrix is not square");
            }
            this->compressed = other.is_compressed();
            this->modified = false;
        };

        /// @brief constructor from a MatrixTransposeView
        /// @note the constructed matrix is in uncompressed format
        /// @param matrix transposed view of matrix to copy
        SquareMatrix(const MatrixTransposeView<T, S> &matrixView);

        /// @brief constructor from a MatrixTransposeView
        /// @note the constructed matrix is in modified compressed format
        /// @param matrix transposed view of matrix to copy
        SquareMatrix(const MatrixDiagonalView<T, S> &matrixView);

        /// @brief default destructor
        virtual ~SquareMatrix() override = default;

        /// @brief check if the matrix is in a modified compressed format
        /// @return true if the matrix is (modified) compressed, false otherwise
        bool is_modified() const { return modified; };

        /// @brief get the size of the modified compressed matrix: it comprehends also possible zero elements in the diagonal
        /// @return size of the modified compressed matrix vectors
        const size_t get_mod_size() const;

        /// @brief compress the matrix in modified format
        void compress_mod();

        /// @brief set an element in the matrix (dynamic construction of the matrix)
        /// @param row row index
        /// @param col column index
        /// @param value value to set
        virtual void set(size_t row, size_t col, const T &value) override;

        /// @brief compress the matrix if it is in an uncompressed format
        virtual void compress() override;

        /// @brief uncompress the matrix if it is in a compressed format
        virtual void uncompress() override;

        /// @brief call operator() const version
        /// @param row row index
        /// @param col column index
        /// @return element at (row, col)
        virtual T operator()(size_t row, size_t col) const override;

        /// @brief call operator() non-const version
        /// @param row row index
        /// @param col column index
        /// @return reference to the element at (row, col) with proxy (to avoid setting zero values)
        virtual Proxy<T, S> operator()(size_t row, size_t col) override;

        /// @brief resize the matrix
        /// @param rows number of rows
        /// @param cols number of columns
        // Overload of Matrix<T, S>::resize_and_clear
        void resize_and_clear(size_t dim);

        /// @brief calculate the norm of the matrix
        /// @tparam N type of the norm (One, Infinity, Frobenius)
        /// @return value of the norm
        /// @note this function HIDES the one in the base class, since it cannot be virtual, being a template function
        template <NormType N>
        double norm() const;

        /// @brief Function to read a matrix in Matrix Market format
        /// @param filename input file name
        virtual void reader(const std::string &filename);

        /// @brief get the number of non-zero elements
        /// @return number of non-zero elements
        virtual size_t get_nnz() const override;

        /// @brief multiply with a std::vector
        /// @tparam U type of the vector elements
        /// @tparam T type of the storage order
        /// @param m matrix
        /// @param v vector
        /// @return the result of the multiplication
        /// @note this function is a friend of the Matrix class, so it can access the private members
        template <AddMulType U, StorageOrder V>
        friend std::vector<U> operator*(const SquareMatrix<U, V> &m, const std::vector<U> &v);

        /// @brief multiply with another matrix
        /// @tparam U type of the matrix elements
        /// @tparam T type of the storage order
        /// @param m1 first matrix
        /// @param m2 second matrix
        /// @return the result of the multiplication
        /// @note this function is a friend of the Matrix class, so it can access the private members
        template <AddMulType U, StorageOrder V>
        friend SquareMatrix<U, V> operator*(const SquareMatrix<U, V> &m1, const SquareMatrix<U, V> &m2);

    private:
        bool modified = false; // flag to check if the matrix is in modified compressed format

        // storage for the matrix
        ModifiedCompressedStorage<T> compressed_format_mod; // MSR or MSC format
    };

};

#include "square_matrix.tpp"

#endif // SQUARE_MATRIX_HPP