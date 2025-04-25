#ifndef MATRIX_VIEWS_HPP
#define MATRIX_VIEWS_HPP

#include "matrix.hpp"
#include "proxy.hpp"

namespace algebra
{
    template <typename T, StorageOrder S = StorageOrder::RowMajor>
    class MatrixTransposeView
    {
    public:
        
        /// @brief delete default constructor
        MatrixTransposeView() = delete;

        /// @brief constructor
        /// @param matrix the matrix to transpose
        MatrixTransposeView(Matrix<T, S> &matrix) : matrix(matrix) {}


        /// @brief call operator non-const version
        /// @param row row index
        /// @param col column index
        /// @return a proxy to the matrix element at (row, col)
        Proxy<T,S> operator()(size_t row, size_t col)
        {
            return matrix(col, row);
        }
        T operator()(size_t row, size_t col) const
        {
            return matrix(col, row);
        }
        
        /// @brief  get the matrix
        /// @return the matrix
        Matrix<T, S> &get_matrix() const
        {
            return matrix;
        }

        // friend functions
        // multiply with a std::vector
        template <AddMulType U, StorageOrder V>
        friend std::vector<U> operator*(const MatrixTransposeView<U, V> &m, const std::vector<U> &v);
        // multiply with another matrix
        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const MatrixTransposeView<U, V> &m1, const Matrix<U, V> &m2);

    private:
        Matrix<T, S> &matrix;
    };

    template <AddMulType T, StorageOrder S>
    std::vector<T> operator*(const MatrixTransposeView<T, S> &m, const std::vector<T> &v)
    {

        if (m.rows != v.size())
        {
            throw std::invalid_argument("Matrix and vector dimensions do not match for multiplication");
        }
        std::vector<T> result(m.cols, 0);
        if (!m.is_compressed())
        {
            for (const auto &it : m.uncompressed_format)
            {
                result[it.first.col] += it.second * v[it.first.row];
            }
        }
        else
        {
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                size_t start;
                size_t end = m.compressed_format.inner[0];
                // iterate over columns of m
                for (size_t col = 0; col < m.cols; col++)
                {
                    start = end;
                    end = m.compressed_format.inner[col + 1];

                    // iterate over rows of m that are non-zero in the column "col" of m
                    for (size_t j = start; j < end; j++)
                    {
                        // row = row of m that we are currently processing
                        size_t row = m.compressed_format.outer[j];

                        // add the product of the non-zero elements to the "result" vector
                        result[col] += m.compressed_format.values[j] * v[row];
                    }
                }
            }
            else
            {
                size_t start;
                size_t end = m.compressed_format.inner[0];
                // iterate over rows of m
                for (size_t row = 0; row < m.rows; row++)
                {
                    start = end;
                    end = m.compressed_format.inner[row + 1];
                    // iterate over columns of m that are non-zero in the row "row" of m
                    for (size_t j = start; j < end; j++)
                    {
                        // col = column of m that we are currently processing
                        size_t col = m.compressed_format.outer[j];

                        // add the product of the non-zero elements to the "result" vector
                        result[col] += m.compressed_format.values[j] * v[row];
                    }
                }
            }
        }
        return result;
    }

    template <AddMulType T, StorageOrder S>
    MatrixTransposeView<T, S> operator*(const MatrixTransposeView<T, S> &m1, const MatrixTransposeView<T, S> &m2)
    {
        return MatrixTransposeView<T,S>(m2.matrix * m1.matrix);
    }
}
#endif // MATRIX_VIEWS_HPP