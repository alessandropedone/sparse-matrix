#ifndef MATRIX_PRODUCTS_TPP
#define MATRIX_PRODUCTS_TPP

#include "storage.hpp"
#include "matrix.hpp"
#include "square_matrix.hpp"
#include "matrix_views.hpp"

namespace algebra
{
    // FRIENDS: MULTIPLICATION WITH TRANSPOSE VIEW
    template <AddMulType T, StorageOrder S>
    std::vector<T> operator*(const TransposeView<T, S> &m, const std::vector<T> &v)
    {

        if (m.matrix.get_rows() != v.size())
        {
            throw std::invalid_argument("Matrix and vector dimensions do not match for multiplication");
        }
        std::vector<T> result(m.matrix.get_cols(), 0);
        if (typeid(m.matrix) == typeid(SquareMatrix<T, S>))
        {
            auto matrix = static_cast<const SquareMatrix<T, S> &>(m.matrix);
            if (matrix.is_modified())
            {
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    size_t col, start, end;

                    // off-diagonal elements
                    for (col = 0; col < matrix.get_cols() - 1; ++col)
                    {
                        // iterate over non-zero elements in the column "col" of m
                        start = matrix.compressed_format_mod.bind[col];
                        end = matrix.compressed_format_mod.bind[col + 1];
                        for (size_t j = start; j < end; ++j)
                        {
                            // row index of the non-zero element
                            size_t row = matrix.compressed_format_mod.bind[j];
                            // add off-diagonal elements
                            result[col] += matrix.compressed_format_mod.values[j] * v[row];
                        }
                    }

                    // handle last column
                    start = matrix.compressed_format_mod.bind[col];
                    end = matrix.compressed_format_mod.values.size();
                    for (size_t j = start; j < end; ++j)
                    {
                        // row index of the non-zero element
                        size_t row = matrix.compressed_format_mod.bind[j];
                        // add off-diagonal elements
                        result[col] += matrix.compressed_format_mod.values[j] * v[row];
                    }

                    // add diagonal elements
                    for (size_t i = 0; i < matrix.get_cols(); ++i)
                    {
                        result[i] += matrix.compressed_format_mod.values[i] * v[i];
                    }
                }
                else
                {
                    size_t row, start, end;

                    // off-diagonal elements
                    for (row = 0; row < matrix.get_rows() - 1; ++row)
                    {
                        // iterate over non-zero elements in the row "row" of m
                        start = matrix.compressed_format_mod.bind[row];
                        end = matrix.compressed_format_mod.bind[row + 1];
                        for (size_t j = start; j < end; ++j)
                        {
                            // col index of the non-zero element
                            size_t col = matrix.compressed_format_mod.bind[j];
                            // add off-diagonal elements
                            result[col] += matrix.compressed_format_mod.values[j] * v[row];
                        }
                    }

                    // handle last row
                    start = matrix.compressed_format_mod.bind[row];
                    end = matrix.compressed_format_mod.values.size();
                    for (size_t j = start; j < end; ++j)
                    {
                        // col index of the non-zero element
                        size_t col = matrix.compressed_format_mod.bind[j];
                        // add off-diagonal elements
                        result[col] += matrix.compressed_format_mod.values[j] * v[row];
                    }

                    // add diagonal elements
                    for (size_t i = 0; i < matrix.cols; ++i)
                    {
                        result[i] += matrix.compressed_format_mod.values[i] * v[i];
                    }
                }
            }
        }
        auto matrix = m.matrix;
        if (not matrix.is_compressed())
        {
            for (const auto &it : matrix.uncompressed_format)
            {
                result[it.first.col] += it.second * v[it.first.row];
            }
        }
        else
        {
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                size_t start;
                size_t end = matrix.compressed_format.inner[0];
                // iterate over columns of m
                for (size_t col = 0; col < matrix.get_cols(); col++)
                {
                    start = end;
                    end = matrix.compressed_format.inner[col + 1];

                    // iterate over rows of matrix that are non-zero in the column "col" of m
                    for (size_t j = start; j < end; j++)
                    {
                        // row = row of m that we are currently processing
                        size_t row = matrix.compressed_format.outer[j];

                        // add the product of the non-zero elements to the "result" vector
                        result[col] += matrix.compressed_format.values[j] * v[row];
                    }
                }
            }
            else
            {
                size_t start;
                size_t end = matrix.compressed_format.inner[0];
                // iterate over rows of m
                for (size_t row = 0; row < matrix.get_rows(); row++)
                {
                    start = end;
                    end = matrix.compressed_format.inner[row + 1];
                    // iterate over columns of m that are non-zero in the row "row" of m
                    for (size_t j = start; j < end; j++)
                    {
                        // col = column of m that we are currently processing
                        size_t col = matrix.compressed_format.outer[j];

                        // add the product of the non-zero elements to the "result" vector
                        result[col] += matrix.compressed_format.values[j] * v[row];
                    }
                }
            }
        }
        return result;
    }

    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const TransposeView<T, S> &m1, const TransposeView<T, S> &m2)
    {
        if (m1.get_cols() != m2.get_rows())
        {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }
        if (m1.is_compressed() != m2.is_compressed())
        {
            throw std::invalid_argument("Matrix compression formats do not match");
        }
        Matrix<T, S> m = m2.matrix * m1.matrix;
        TransposeView<T, S> t(m);
        Matrix<T, S> result(t);
        return result;
    }

    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const TransposeView<T, S> &m1, const Matrix<T, S> &m2)
    {
        // TO DO

        Matrix<T, S> result(m1.get_rows, m2.get_cols);

        if (m1.get_cols() != m2.get_rows())
        {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }
        if (m1.is_compressed() != m2.is_compressed())
        {
            throw std::invalid_argument("Matrix compression formats do not match");
        }

        if (not m1.is_compressed())
        {
            for (const auto &it1 : m1.uncompressed_format)
            {
                for (const auto &it2 : m2.uncompressed_format)
                {
                    if (it1.first.col == it2.first.row)
                    {
                        result(it1.first.row, it2.first.col) += it1.second * it2.second;
                    }
                }
            }
        }

        return result;
    }

    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const Matrix<T, S> &m1, const TransposeView<T, S> &m2)
    {
        return Matrix<T, S>(m1.get_rows(), m2.get_cols());
    }

    // FRIENDS: MULTIPLICATIONS WITH DIAGONAL VIEWS
    template <AddMulType T, StorageOrder S>
    std::vector<T> operator*(const DiagonalView<T, S> &m, const std::vector<T> &v)
    {
        if (m.get_rows() != v.size())
        {
            throw std::invalid_argument("Matrix and vector dimensions do not match for multiplication");
        }
        std::vector<T> result(m.get_rows(), 0);
        return result;
    };

    template <AddMulType T, StorageOrder S>
    SquareMatrix<T, S> operator*(const DiagonalView<T, S> &m1, const DiagonalView<T, S> &m2)
    {
        if (m1.get_rows() != m2.get_rows())
        {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }
        if ((m1.is_modified() != m2.is_modified()) || (m1.is_compressed() != m2.is_compressed()))
        {
            throw std::invalid_argument("Matrix compression formats do not match");
        }

        SquareMatrix<T, S> result(m1.get_rows());
        return result;
    };

    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const Matrix<T, S> &m1, const DiagonalView<T, S> &m2)
    {
        return Matrix<T, S>(m1.get_rows(), m2.get_cols());
    };

    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const DiagonalView<T, S> &m1, const Matrix<T, S> &m2)
    {
        return Matrix<T, S>(m1.get_rows, m2.get_cols());
    };
}

#endif // PRODUCTS_TPP