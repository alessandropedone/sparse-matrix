#ifndef ABSTRACT_MATRIX_HPP
#define ABSTRACT_MATRIX_HPP

#include "storage.hpp"

namespace algebra
{
    template <typename T, StorageOrder S = StorageOrder::RowMajor>
    class AbstractMatrix
    {
    public:
        virtual ~AbstractMatrix() = default;

        /// @brief set an element in the matrix (dynamic construction of the matrix)
        /// @param row row index
        /// @param col column index
        /// @param value value to set
        virtual void set(size_t row, size_t col, const T &value) = 0;

        /// @brief check if the matrix is in a compressed format
        /// @return true if the matrix is compressed, false otherwise
        virtual bool is_compressed() = 0;

        /// @brief call operator() const version
        /// @param row row index
        /// @param col column index
        /// @return element at (row, col)
        virtual T operator()(size_t row, size_t col) const = 0;

        /// @brief call operator() non-const version
        /// @param row row index
        /// @param col column index
        /// @return reference to the element at (row, col) with proxy (to avoid storing zero values)
        virtual Proxy<T, S> operator()(size_t row, size_t col) = 0;

        /// @brief get the number of rows
        /// @return number of rows
        virtual size_t get_rows() const = 0;

        /// @brief get the number of columns
        /// @return number of columns
        virtual size_t get_cols() const = 0;

        /// @brief get the number of non-zero elements
        /// @return number of non-zero elements
        virtual size_t get_nnz() const = 0;
    };
}
#endif // ABSTRACT_MATRIX_HPP