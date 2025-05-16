
/**
 * @file matrix.hpp
 * @author 
 * @brief Defines the Matrix class and related algebraic structures for sparse and dense matrix operations.
 * 
 * This header provides the main Matrix class template, supporting both row-major and column-major storage orders,
 * and a variety of algebraic operations such as multiplication, transposition, and norm calculation. The Matrix
 * class is designed for efficient manipulation of both sparse and dense matrices, with support for compressed and
 * uncompressed storage formats. The file also includes forward declarations for TransposeView and DiagonalView,
 * and declares friend functions for matrix and vector products.
 * 
 * @details
 * The Matrix class inherits from AbstractMatrix and provides a flexible interface for matrix construction,
 * manipulation, and algebraic operations. It supports dynamic and static memory management, parallel compression,
 * and reading from Matrix Market files. The file is part of the algebra namespace, which encapsulates all
 * algebraic matrix operations and representations.
 * 
 * @see storage.hpp
 * @see proxy.hpp
 * @see abstract_matrix.hpp
 * @see matrix.tpp
 * @see view_products.tpp
 */
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "storage.hpp"
#include "proxy.hpp"
#include "abstract_matrix.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <execution>
#include <cmath>

/**
 * @namespace algebra
 * @brief Contains classes and functions for algebraic operations on matrices, including
 *        various matrix representations, views, and operations such as multiplication,
 *        transposition, and norm calculation.
 *
 * The algebra namespace provides a flexible and extensible framework for working with
 * sparse and dense matrices in different storage formats (row-major, column-major),
 * supporting efficient mathematical operations and matrix manipulations.
 */
namespace algebra
{

    // forward declaration of the TransposeView class
    template <AddMulType T, StorageOrder S>
    class TransposeView;

    // forward declaration of the DiagonalView class
    template <AddMulType T, StorageOrder S>
    class DiagonalView;

    /**
     * @class Matrix
     * @brief Represents a sparse matrix with configurable storage order and element type.
     * 
     * The Matrix class provides an interface for constructing, manipulating, and performing operations
     * on sparse matrices. It supports both row-major and column-major storage orders, and can store
     * elements of any type that supports addition and multiplication. The class allows for dynamic
     * construction, compression (to CSR/CSC), and uncompression (to COO), as well as efficient
     * multiplication with vectors and other matrices. It also provides views for transposed and diagonal
     * matrices, and supports parallel compression.
     * 
     * @tparam T Type of the matrix elements. Must satisfy AddMulType concept.
     * @tparam S Storage order of the matrix (RowMajor or ColumnMajor). Defaults to RowMajor.
     * 
     * @note The Matrix class inherits from AbstractMatrix and is not default-constructible.
     * 
     * @see AbstractMatrix
     * @see TransposeView
     * @see DiagonalView
     */
    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class Matrix : public AbstractMatrix<T, S>
    {
    public:
        // delete default constructor
        Matrix() = delete;

        /// @brief constructor with size
        /// @param rows number of rows
        /// @param cols number of columns
        Matrix(size_t rows, size_t cols) : rows(rows), cols(cols) { this->compressed = false; };

        /// @brief default copy constructor
        Matrix(const Matrix &other) = default;

        /// @brief constructor from a TransposeView
        /// @note the constructed matrix is in uncompressed format
        /// @param view transposed view of matrix to copy
        Matrix(const TransposeView<T, S> &view);

        /// @brief constructor from a DiagonalView
        /// @note the constructed matrix is in uncompressed format
        /// @param view diagonal view of matrix to copy
        Matrix(const DiagonalView<T, S> &view);

        /// @brief clone method
        /// @return a pointer to the cloned object
        virtual std::unique_ptr<AbstractMatrix<T, S>> clone() const override
        {
            return std::make_unique<Matrix<T, S>>(*this);
        };

        /// @brief default destructor
        virtual ~Matrix() = default;

        /// @brief move constructor
        /// @param other matrix to move
        Matrix(Matrix &&other) noexcept;

        /// @brief move assignment operator
        /// @param other matrix to move
        /// @return reference to the moved matrix
        Matrix &operator=(Matrix &&other) noexcept;

        /// @brief set an element in the matrix (dynamic construction of the matrix)
        /// @param row row index
        /// @param col column index
        /// @param value value to set
        virtual void set(size_t row, size_t col, const T &value) override;

        /// @brief check if the matrix is in a compressed format
        /// @return true if the matrix is compressed, false otherwise
        virtual bool is_compressed() const override { return compressed; };

        /// @brief compress the matrix if it is in an uncompressed format
        virtual void compress() override;

        /// @brief compress the matrix in parallel if it is in an uncompressed format
        virtual void compress_parallel();

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
        /// @return reference to the element at (row, col) with proxy (to avoid storing zero values)
        virtual Proxy<T, S> operator()(size_t row, size_t col) override;

        /// @brief resize the matrix
        /// @param rows number of rows
        /// @param cols number of columns
        virtual void resize_and_clear(size_t rows, size_t cols);

        /// @brief calculate the norm of the matrix
        /// @tparam N type of the norm (One, Infinity, Frobenius)
        /// @return value of the norm
        template <NormType N>
        double norm() const;

        /// @brief Function to read a matrix in Matrix Market format
        /// @param filename input file name
        virtual void reader(const std::string &filename) override;

        /// @brief get the number of rows
        /// @return number of rows
        virtual size_t get_rows() const override { return rows; };

        /// @brief get the number of columns
        /// @return number of columns
        virtual size_t get_cols() const override { return cols; };

        /// @brief get the number of non-zero elements
        /// @return number of non-zero elements
        virtual size_t get_nnz() const override;

        /// @brief multiply with a std::vector
        /// @tparam U type of the vector elements
        /// @tparam V type of the storage order
        /// @param m matrix
        /// @param v vector
        /// @return the result of the multiplication
        /// @note this function is a friend of the Matrix class, so it can access the private members
        template <AddMulType U, StorageOrder V>
        friend std::vector<U> operator*(const Matrix<U, V> &m, const std::vector<U> &v);

        /// @brief multiply with another matrix
        /// @tparam U type of the matrix elements
        /// @tparam V type of the storage order
        /// @param m1 first matrix
        /// @param m2 second matrix
        /// @return the result of the multiplication
        /// @note this function is a friend of the Matrix class, so it can access the private members
        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const Matrix<U, V> &m1, const Matrix<U, V> &m2);

        /// @brief multiply a TransposeView with a std::vector
        /// @tparam U type of the vector elements
        /// @tparam V type of the storage order
        /// @param m matrix
        /// @param v vector
        /// @return the result of the multiplication
        template <AddMulType U, StorageOrder V>
        friend std::vector<U> operator*(const TransposeView<U, V> &m, const std::vector<U> &v);

        /// @brief multiply two TransposeViews
        /// @tparam U type of the matrix elements
        /// @tparam V type of the storage order
        /// @param m1 first matrix
        /// @param m2 second matrix
        /// @return the result of the multiplication
        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const TransposeView<U, V> &m1, const TransposeView<U, V> &m2);

        /// @brief multiply a DiagonalView with a std::vector
        /// @tparam U type of the vector elements
        /// @tparam V type of the storage order
        /// @param m matrix
        /// @param v vector
        /// @return the result of the multiplication
        template <AddMulType U, StorageOrder V>
        friend std::vector<U> operator*(const DiagonalView<U, V> &m, const std::vector<U> &v);

        /// @brief multiply two DiagonalViews
        /// @tparam U type of the matrix elements
        /// @tparam V type of the storage order
        /// @param m1 first matrix
        /// @param m2 second matrix
        /// @return the result of the multiplication
        template <AddMulType U, StorageOrder V>
        friend SquareMatrix<U, V> operator*(const DiagonalView<U, V> &m1, const DiagonalView<U, V> &m2);

        /// @brief multiply a DiagonalView with a Matrix
        /// @tparam U type of the matrix elements
        /// @tparam V type of the storage order
        /// @param m1 first matrix
        /// @param m2 second matrix
        /// @return the result of the multiplication
        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const Matrix<U, V> &m1, const DiagonalView<U, V> &m2);

        /// @brief multiply a DiagonalView with a Matrix
        /// @tparam U type of the matrix elements
        /// @tparam V type of the storage order
        /// @param m1 first matrix
        /// @param m2 second matrix
        /// @return the result of the multiplication
        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const DiagonalView<U, V> &m1, const Matrix<U, V> &m2);

    protected:
        size_t rows;             /// number of rows
        size_t cols;             /// number of columns
        bool compressed = false; /// flag to check if the matrix is compressed

        // storage for the matrix
        // uncompressed matrix
        UncompressedStorage<T, S> uncompressed_format; /// COO format
        // compressed matrix
        CompressedStorage<T> compressed_format; /// CSR or CSC format
    };

}

#include "matrix.tpp"
#include "view_products.tpp"

#endif // MATRIX_HPP
