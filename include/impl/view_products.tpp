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

        // move semantic
        Matrix<T, S> result(m1.get_rows(), m2.get_cols());

        const auto *square_matrix1 = dynamic_cast<const SquareMatrix<T, S> *>(&m1.matrix);
        const auto *square_matrix2 = dynamic_cast<const SquareMatrix<T, S> *>(&m2.matrix);
        if (square_matrix1 && square_matrix2)
        {
            if (square_matrix1->is_modified() && square_matrix2->is_modified())
            {
                auto matrix1 = static_cast<const SquareMatrix<T, S> &>(m1.matrix);
                auto matrix2 = static_cast<const SquareMatrix<T, S> &>(m2.matrix);
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    size_t col;

                    ////// ITERATE OVER COLUMNS OF matrix1 //////

                    for (col = 0; col < matrix1.cols - 1; ++col)
                    {
                        // ADD OFF-DIAGONAL ELEMENTS OF matrix1 AND matrix2
                        size_t start = matrix1.compressed_format_mod.bind[col];
                        size_t end = matrix1.compressed_format_mod.bind[col + 1];
                        size_t k;
                        // iterate over rows of matrix1 (and columns of matrix2) that are non-zero in the column "col" of matrix1
                        for (k = start; k < end - 1; ++k)
                        {
                            // j = row of matrix1 (or column of matrix2) that we are currently processing
                            size_t j = matrix1.compressed_format_mod.bind[k];

                            // iterate over rows of matrix2 that are non-zero in the column j of matrix2
                            size_t s = matrix2.compressed_format_mod.bind[j];
                            size_t e = matrix2.compressed_format_mod.bind[j + 1];
                            for (size_t i = s; i < e; ++i)
                            {
                                // row = row of matrix2 corresponding to the index i
                                size_t row = matrix2.compressed_format_mod.bind[i];

                                // add the product of the non-zero off-diagonal elements to the "result" matrix
                                result(col, row) += matrix2.compressed_format_mod.values[i] * matrix1.compressed_format_mod.values[k];
                            }
                        }
                        // handle last column of matrix2 (last row of matrix1)
                        size_t j = matrix1.compressed_format_mod.bind[k];
                        size_t s = matrix2.compressed_format_mod.bind[j];
                        size_t e = (j + 1 == matrix2.cols) ? matrix2.compressed_format_mod.values.size() : matrix2.compressed_format_mod.bind[j + 1];
                        for (size_t i = s; i < e; ++i)
                        {
                            // row = row of matrix2 corresponding to the index i
                            size_t row = matrix2.compressed_format_mod.bind[i];

                            // add the product of the non-zero off-diagonal elements to the "result" matrix
                            result(col, row) += matrix2.compressed_format_mod.values[i] * matrix1.compressed_format_mod.values[k];
                        }

                        // ADD DIAGONAL ELEMENTS OF matrix2
                        // iterate over rows of matrix1 (and columns of matrix2) that are non-zero in the column "col" of matrix1
                        for (k = start; k < end; ++k)
                        {
                            // j = row of matrix1 (or column of matrix2) that we are currently processing
                            size_t j = matrix1.compressed_format_mod.bind[k];

                            // add the product between the diagonal element of matrix2 and current non-zero element of matrix1 to the "result" matrix
                            result(col, j) += matrix2.compressed_format_mod.values[j] * matrix1.compressed_format_mod.values[k];
                        }
                    }

                    ////// HANDLE LAST COLUMN OF matrix1 //////

                    // ADD OFF-DIAGONAL ELEMENTS OF matrix1 AND matrix2
                    size_t start = matrix1.compressed_format_mod.bind[col];
                    size_t end = matrix1.compressed_format_mod.values.size();
                    size_t k;
                    // iterate over rows of matrix1 (and columns of matrix2) that are non-zero in the column "col" of matrix1
                    for (k = start; k < end - 1; ++k)
                    {
                        // j = row of matrix1 (or column of matrix2) that we are currently processing
                        size_t j = matrix1.compressed_format_mod.bind[k];

                        // iterate over rows of matrix2 that are non-zero in the column j of matrix2
                        size_t s = matrix2.compressed_format_mod.bind[j];
                        size_t e = matrix2.compressed_format_mod.bind[j + 1];
                        for (size_t i = s; i < e; ++i)
                        {
                            // row = row of matrix2 corresponding to the index i
                            size_t row = matrix2.compressed_format_mod.bind[i];

                            // add the product of the non-zero off-diagonal elements to the "result" matrix
                            result(col, row) += matrix2.compressed_format_mod.values[i] * matrix1.compressed_format_mod.values[k];
                        }
                    }
                    // handle last column of matrix2 (last row of matrix1)
                    size_t j = matrix1.compressed_format_mod.bind[k];
                    size_t s = matrix2.compressed_format_mod.bind[j];
                    size_t e = (j + 1 == matrix2.cols) ? matrix2.compressed_format_mod.values.size() : matrix2.compressed_format_mod.bind[j + 1];
                    for (size_t i = s; i < e; ++i)
                    {
                        // row = row of matrix2 corresponding to the index i
                        size_t row = matrix2.compressed_format_mod.bind[i];

                        // add the product of the non-zero off-diagonal elements to the "result" matrix
                        result(col, row) += matrix2.compressed_format_mod.values[i] * matrix1.compressed_format_mod.values[k];
                    }

                    // ADD DIAGONAL ELEMENTS OF matrix2
                    // iterate over rows of matrix1 (and columns of matrix2) that are non-zero in the column "col" of matrix1
                    for (size_t k = start; k < end; ++k)
                    {
                        // j = row of matrix1 (or column of matrix2) that we are currently processing
                        size_t j = matrix1.compressed_format_mod.bind[k];

                        // add the product between the diagonal element of matrix2 and current non-zero element of matrix1 to the "result" matrix
                        result(col, j) += matrix2.compressed_format_mod.values[j] * matrix1.compressed_format_mod.values[k];
                    }

                    ////// ADD DIAGONAL ELEMENTS OF matrix1 //////

                    for (col = 0; col < matrix1.cols - 1; ++col)
                    {
                        // multiply each column of matrix2 with the diagonal element of matrix1
                        size_t s = matrix2.compressed_format_mod.bind[col];
                        size_t e = matrix2.compressed_format_mod.bind[col + 1];
                        // iterate over the non-zero elements of matrix2 that are in the column "col"
                        for (size_t i = s; i < e; ++i)
                        {
                            // row = row of matrix2 corresponding to the index i
                            size_t row = matrix2.compressed_format_mod.bind[i];

                            // add the product between the diagonal element of matrix1 and current non-zero element of matrix2 to the "result" matrix
                            result(col, row) += matrix2.compressed_format_mod.values[i] * matrix1.compressed_format_mod.values[col];
                        }
                    }
                    // HANDLE LAST COLUMN OF matrix2
                    // multiply each column of matrix2 with the diagonal element of matrix1
                    s = matrix2.compressed_format_mod.bind[col];
                    e = matrix2.compressed_format_mod.values.size();
                    // iterate over the non-zero elements of matrix2 that are in the column "col"
                    for (size_t i = s; i < e; ++i)
                    {
                        // row = row of matrix2 corresponding to the index i
                        size_t row = matrix2.compressed_format_mod.bind[i];

                        // add the product between the diagonal element of matrix1 and current non-zero element of matrix2 to the "result" matrix
                        result(col, row) += matrix2.compressed_format_mod.values[i] * matrix1.compressed_format_mod.values[col];
                    }

                    ////// ADD PRODUCTS OF DIAGONAL ELEMENTS //////

                    // iterate over the diagonal elements of matrix2 and matrix1
                    for (size_t i = 0; i < matrix2.rows; ++i)
                    {
                        // add the product between the diagonal elements of matrix2 and matrix1 to the "result" matrix
                        result(i, i) += matrix2.compressed_format_mod.values[i] * matrix1.compressed_format_mod.values[i];
                    }
                }
                else
                {
                    size_t row;

                    ////// ITERATE OVER COLUMNS OF matrix2 //////

                    for (row = 0; row < matrix2.rows - 1; ++row)
                    {
                        // ADD OFF-DIAGONAL ELEMENTS OF matrix2 AND matrix1
                        size_t start = matrix2.compressed_format_mod.bind[row];
                        size_t end = matrix2.compressed_format_mod.bind[row + 1];
                        size_t k;
                        // iterate over columns of matrix2 (and rows of matrix1) that are non-zero in the row "row" of matrix2
                        for (k = start; k < end - 1; ++k)
                        {
                            // j = column of matrix2 (or row of matrix1) that we are currently processing
                            size_t j = matrix2.compressed_format_mod.bind[k];

                            // iterate over columns of matrix1 that are non-zero in the row j of matrix1
                            size_t s = matrix1.compressed_format_mod.bind[j];
                            size_t e = matrix1.compressed_format_mod.bind[j + 1];
                            for (size_t i = s; i < e; ++i)
                            {
                                // col = column of matrix1 corresponding to the index i
                                size_t col = matrix1.compressed_format_mod.bind[i];

                                // add the product of the non-zero off-diagonal elements to the "result" matrix
                                result(col, row) += matrix2.compressed_format_mod.values[k] * matrix1.compressed_format_mod.values[i];
                            }
                        }
                        // handle last row of matrix1 (last column of matrix2)
                        size_t j = matrix2.compressed_format_mod.bind[k];
                        size_t s = matrix1.compressed_format_mod.bind[j];
                        size_t e = (j + 1 == matrix1.rows) ? matrix1.compressed_format_mod.values.size() : matrix1.compressed_format_mod.bind[j + 1];
                        for (size_t i = s; i < e; ++i)
                        {
                            // col = column of matrix1 corresponding to the index i
                            size_t col = matrix1.compressed_format_mod.bind[i];

                            // add the product of the non-zero off-diagonal elements to the "result" matrix
                            result(col, row) += matrix2.compressed_format_mod.values[k] * matrix1.compressed_format_mod.values[i];
                        }

                        // ADD DIAGONAL ELEMENTS OF matrix1
                        // iterate over columns of matrix2 (and rows of matrix1) that are non-zero in the row "row" of matrix2
                        for (k = start; k < end; ++k)
                        {
                            // j = column of matrix2 (or row of matrix1) that we are currently processing
                            size_t j = matrix2.compressed_format_mod.bind[k];

                            // add the product between the diagonal element of matrix1 and current non-zero element of matrix2 to the "result" matrix
                            result(j, row) += matrix2.compressed_format_mod.values[k] * matrix1.compressed_format_mod.values[j];
                        }
                    }

                    ///// HANDLE LAST ROW OF matrix2 /////

                    // ADD OFF-DIAGONAL ELEMENTS OF matrix2 AND matrix1
                    size_t start = matrix2.compressed_format_mod.bind[row];
                    size_t end = matrix2.compressed_format_mod.values.size();
                    size_t k;
                    // iterate over columns of matrix2 (and rows of matrix1) that are non-zero in the row "row" of matrix2
                    for (k = start; k < end - 1; ++k)
                    {
                        // j = column of matrix2 (or row of matrix1) that we are currently processing
                        size_t j = matrix2.compressed_format_mod.bind[k];

                        // iterate over columns of matrix1 that are non-zero in the row j of matrix1
                        size_t s = matrix1.compressed_format_mod.bind[j];
                        size_t e = matrix1.compressed_format_mod.bind[j + 1];
                        for (size_t i = s; i < e; ++i)
                        {
                            // col = column of matrix1 corresponding to the index i
                            size_t col = matrix1.compressed_format_mod.bind[i];

                            // add the product of the non-zero off-diagonal elements to the "result" matrix
                            result(col, row) += matrix2.compressed_format_mod.values[k] * matrix1.compressed_format_mod.values[i];
                        }
                    }
                    // handle last row of matrix1 (last column of matrix2)
                    size_t j = matrix2.compressed_format_mod.bind[k];
                    size_t s = matrix1.compressed_format_mod.bind[j];
                    size_t e = (j + 1 == matrix1.rows) ? matrix1.compressed_format_mod.values.size() : matrix1.compressed_format_mod.bind[j + 1];
                    for (size_t i = s; i < e; ++i)
                    {
                        // col = column of matrix1 corresponding to the index i
                        size_t col = matrix1.compressed_format_mod.bind[i];

                        // add the product of the non-zero off-diagonal elements to the "result" matrix
                        result(col, row) += matrix2.compressed_format_mod.values[k] * matrix1.compressed_format_mod.values[i];
                    }

                    // ADD DIAGONAL ELEMENTS OF matrix1
                    // iterate over columns of matrix2 (and rows of matrix1) that are non-zero in the row "row" of matrix2
                    for (size_t k = start; k < end; ++k)
                    {
                        // j = column of matrix2 (or row of matrix1) that we are currently processing
                        size_t j = matrix2.compressed_format_mod.bind[k];

                        // add the product between the diagonal element of matrix1 and current non-zero element of matrix2 to the "result" matrix
                        result(j, row) += matrix2.compressed_format_mod.values[k] * matrix1.compressed_format_mod.values[j];
                    }

                    ////// ADD DIAGONAL ELEMENTS OF matrix2 //////

                    for (row = 0; row < matrix2.rows - 1; ++row)
                    {
                        // multiply each row of matrix1 with the diagonal element of matrix2
                        size_t s = matrix1.compressed_format_mod.bind[row];
                        size_t e = matrix1.compressed_format_mod.bind[row + 1];
                        // iterate over the non-zero elements of matrix1 that are in the row "row"
                        for (size_t i = s; i < e; ++i)
                        {
                            // col = column of matrix1 corresponding to the index i
                            size_t col = matrix1.compressed_format_mod.bind[i];

                            // add the product between the diagonal element of matrix2 and current non-zero element of matrix1 to the "result" matrix
                            result(col, row) += matrix2.compressed_format_mod.values[row] * matrix1.compressed_format_mod.values[i];
                        }
                    }
                    // HANDLE LAST ROW OF matrix1
                    // multiply each row of matrix1 with the diagonal element of matrix2
                    s = matrix1.compressed_format_mod.bind[row];
                    e = matrix1.compressed_format_mod.values.size();
                    // iterate over the non-zero elements of matrix1 that are in the row "row"
                    for (size_t i = s; i < e; ++i)
                    {
                        // col = column of matrix1 corresponding to the index i
                        size_t col = matrix1.compressed_format_mod.bind[i];

                        // add the product between the diagonal element of matrix2 and current non-zero element of matrix1 to the "result" matrix
                        result(col, row) += matrix2.compressed_format_mod.values[row] * matrix1.compressed_format_mod.values[i];
                    }

                    ////// ADD PRODUCTS OF DIAGONAL ELEMENTS //////

                    // iterate over the diagonal elements of matrix2 and matrix1
                    for (size_t i = 0; i < matrix2.rows; ++i)
                    {
                        // add the product between the diagonal elements of matrix2 and matrix1 to the "result" matrix
                        result(i, i) += matrix2.compressed_format_mod.values[i] * matrix1.compressed_format_mod.values[i];
                    }
                }
            }
            else if (square_matrix1->is_modified() or square_matrix2->is_modified())
            {
                throw std::invalid_argument("Both matrices must be square and in modified format for this operation");
            }
        }
        else if (m1.is_compressed())
        {
            auto &matrix1 = m1.matrix;
            auto &matrix2 = m2.matrix;
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                // iterate over columns of matrix1
                for (size_t col = 0; col < matrix1.cols; col++)
                {
                    // iterate over rows of matrix1 (and columns of matrix2) that are non-zero in the column "col" of matrix1
                    for (size_t k = matrix1.compressed_format.inner[col]; k < matrix1.compressed_format.inner[col + 1]; k++)
                    {
                        // j = row of matrix1 (or column of matrix2) that we are currently processing
                        size_t j = matrix1.compressed_format.outer[k];

                        // iterate over rows of matrix2 that are non-zero in the column j of matrix2
                        for (size_t i = matrix2.compressed_format.inner[j]; i < matrix2.compressed_format.inner[j + 1]; i++)
                        {
                            // row = row of matrix2 corresponding to the index i
                            size_t row = matrix2.compressed_format.outer[i];

                            // add the product of the non-zero elements to the "result" matrix
                            result(col, row) += matrix2.compressed_format.values[i] * matrix1.compressed_format.values[k];
                        }
                    }
                }
            }
            else
            {
                // iterate over rows of matrix2
                for (size_t row = 0; row < matrix2.rows; row++)
                {
                    // iterate over columns of matrix2 (and rows of matrix1) that are non-zero in the row "row" of matrix2
                    for (size_t j = matrix2.compressed_format.inner[row]; j < matrix2.compressed_format.inner[row + 1]; j++)
                    {
                        // k = column of matrix2 (or row of matrix1) that we are currently processing
                        size_t k = matrix2.compressed_format.outer[j];

                        // iterate over columns of matrix1 that are non-zero in the row k of matrix1
                        for (size_t i = matrix1.compressed_format.inner[k]; i < matrix1.compressed_format.inner[k + 1]; i++)
                        {
                            // col = column of matrix1 corresponding to the index i
                            size_t col = matrix1.compressed_format.outer[i];

                            // add the product of the non-zero elements to the "result" matrix
                            result(col, row) += matrix2.compressed_format.values[j] * matrix1.compressed_format.values[i];
                        }
                    }
                }
            }
        }
        else
        {
            auto &matrix1 = m1.matrix;
            auto &matrix2 = m2.matrix;
            for (const auto &it1 : matrix1.uncompressed_format)
            {
                for (const auto &it2 : matrix2.uncompressed_format)
                {
                    if (it1.first.row == it2.first.col)
                    {
                        result(it1.first.col, it2.first.row) += it1.second * it2.second;
                    }
                }
            }
        }
        return result;
    }



    // FRIENDS: MULTIPLICATIONS WITH DIAGONAL VIEWS

    template <AddMulType T, StorageOrder S>
    std::vector<T> operator*(const DiagonalView<T, S> &m, const std::vector<T> &v)
    {
        if (m.get_cols() != v.size())
        {
            throw std::invalid_argument("Matrix and vector dimensions do not match for multiplication");
        }
        std::vector<T> result(m.get_rows(), 0);
        auto &matrix = m.matrix;
        if (matrix.is_modified())
        {
            size_t cols = matrix.get_cols();
            for (size_t i = 0; i < cols; ++i)
            {
                result[i] += matrix.compressed_format_mod.values[i] * v[i];
            }
        }
        else if (matrix.is_compressed())
        {
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                size_t start;
                size_t end = matrix.compressed_format.inner[0];
                size_t cols = matrix.get_cols();
                // iterate over columns of m
                for (size_t col = 0; col < cols; col++)
                {
                    start = end;
                    end = matrix.compressed_format.inner[col + 1];

                    // iterate over rows of matrix that are non-zero in the column "col" of m
                    for (size_t j = start; j < end; j++)
                    {
                        // row = row of m that we are currently processing
                        size_t row = matrix.compressed_format.outer[j];

                        if (col == row)
                        {
                            // add the product of the non-zero elements to the "result" vector
                            result[row] += matrix.compressed_format.values[j] * v[col];
                        }
                    }
                }
            }
            else
            {
                size_t start;
                size_t end = matrix.compressed_format.inner[0];
                size_t rows = matrix.get_cols();
                // iterate over rows of m
                for (size_t row = 0; row < rows; row++)
                {
                    start = end;
                    end = matrix.compressed_format.inner[row + 1];
                    // iterate over columns of m that are non-zero in the row "row" of m
                    for (size_t j = start; j < end; j++)
                    {
                        // col = column of m that we are currently processing
                        size_t col = matrix.compressed_format.outer[j];

                        if (col == row)
                        {
                            // add the product of the non-zero elements to the "result" vector
                            result[row] += matrix.compressed_format.values[j] * v[col];
                        }
                    }
                }
            }
        }
        else
        {
            for (const auto &it : matrix.uncompressed_format)
            {
                if (it.first.row == it.first.col)
                {
                    // add the product of the non-zero elements to the "result" vector
                    result[it.first.row] += it.second * v[it.first.col];
                }
            }
        }
        return result;
    };

    template <AddMulType T, StorageOrder S>
    SquareMatrix<T, S> operator*(const DiagonalView<T, S> &m1, const DiagonalView<T, S> &m2)
    {
        if (m1.get_cols() != m2.get_rows())
        {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }
        if ((m1.is_modified() != m2.is_modified()) || (m1.is_compressed() != m2.is_compressed()))
        {
            throw std::invalid_argument("Matrix compression formats do not match");
        }
        SquareMatrix<T, S> result(m1.get_rows());
        auto &matrix1 = m1.matrix;
        auto &matrix2 = m2.matrix;
        if (matrix1.is_modified())
        {
            size_t cols = matrix1.get_cols();
            for (size_t i = 0; i < cols; ++i)
            {
                result(i, i) = matrix1.compressed_format_mod.values[i] * matrix2.compressed_format_mod.values[i];
            }
        }
        else if (matrix1.is_compressed())
        {
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                // iterate over columns of matrix2
                for (size_t col = 0; col < matrix2.cols; col++)
                {
                    // iterate over rows of matrix2 (and columns of matrix1) that are non-zero in the column "col" of matrix2
                    for (size_t k = matrix2.compressed_format.inner[col]; k < matrix2.compressed_format.inner[col + 1]; k++)
                    {
                        // j = row of matrix2 (or column of matrix1) that we are currently processing
                        size_t j = matrix2.compressed_format.outer[k];

                        // if we are on the diagonal of matrix2
                        if (j == col)
                        {
                            // iterate over rows of matrix1 that are non-zero in the column j of matrix1
                            for (size_t i = matrix1.compressed_format.inner[j]; i < matrix1.compressed_format.inner[j + 1]; i++)
                            {
                                // row = row of matrix1 corresponding to the index i
                                size_t row = matrix1.compressed_format.outer[i];

                                // if we are on the diagonal of matrix1
                                if (row == j)
                                {
                                    // add the product of the non-zero elements to the "result" matrix
                                    result(row, col) += matrix1.compressed_format.values[i] * matrix2.compressed_format.values[k];
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                // iterate over rows of matrix1
                for (size_t row = 0; row < matrix1.rows; row++)
                {
                    // iterate over columns of matrix1 (and rows of matrix2) that are non-zero in the row "row" of matrix1
                    for (size_t j = matrix1.compressed_format.inner[row]; j < matrix1.compressed_format.inner[row + 1]; j++)
                    {
                        // k = column of matrix1 (or row of matrix2) that we are currently processing
                        size_t k = matrix1.compressed_format.outer[j];

                        // if we are the diagonal of matrix1
                        if (k == row)
                        {
                            // iterate over columns of matrix2 that are non-zero in the row k of matrix2
                            for (size_t i = matrix2.compressed_format.inner[k]; i < matrix2.compressed_format.inner[k + 1]; i++)
                            {
                                // col = column of matrix2 corresponding to the index i
                                size_t col = matrix2.compressed_format.outer[i];

                                // if we are on the diagonal of matrix2
                                if (col == k)
                                {
                                    // add the product of the non-zero elements to the "result" matrix
                                    result(row, col) += matrix1.compressed_format.values[j] * matrix2.compressed_format.values[i];
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            for (const auto &it1 : matrix1.uncompressed_format)
            {
                for (const auto &it2 : matrix2.uncompressed_format)
                {
                    if (it1.first.col == it2.first.row && it1.first.row == it1.first.col && it2.first.row == it2.first.col)
                    {
                        // add the product of the non-zero elements to the "result" vector
                        result(it1.first.row, it1.first.col) += it1.second * it2.second;
                    }
                }
            }
        }
        return result;
    };

    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const Matrix<T, S> &m1, const DiagonalView<T, S> &m2)
    {
        if (m1.get_cols() != m2.get_rows())
        {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }
        if (m1.is_compressed() != m2.is_compressed())
        {
            throw std::invalid_argument("Matrix compression formats do not match");
        }
        Matrix<T, S> result(m1.get_rows(), m2.get_cols());
        auto &matrix2 = m2.matrix;

        if (matrix2.is_modified())
        {
            if (auto square_matrix = dynamic_cast<const SquareMatrix<T, S> *>(&m1))
            {
                auto &matrix1 = *square_matrix;
                if (matrix1.is_modified())
                {
                    size_t cols = matrix1.get_cols();
                    for (size_t i = 0; i < cols; ++i)
                    {
                        result(i, i) = matrix1.compressed_format_mod.values[i] * matrix2.compressed_format_mod.values[i];
                    }
                }
                else
                {
                    throw std::invalid_argument("Matrix compression formats do not match");
                }
            }
            else
            {
                throw std::invalid_argument("Matrix compression formats do not match");
            }
        }
        else if (matrix2.is_compressed())
        {
            auto &matrix1 = m1;
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                // iterate over columns of matrix2
                for (size_t col = 0; col < matrix2.cols; col++)
                {
                    // iterate over rows of matrix2 (and columns of matrix1) that are non-zero in the column "col" of matrix2
                    for (size_t k = matrix2.compressed_format.inner[col]; k < matrix2.compressed_format.inner[col + 1]; k++)
                    {
                        // j = row of matrix2 (or column of matrix1) that we are currently processing
                        size_t j = matrix2.compressed_format.outer[k];

                        // if we are on the diagonal of matrix2
                        if (j == col)
                        {
                            // iterate over rows of matrix1 that are non-zero in the column j of matrix1
                            for (size_t i = matrix1.compressed_format.inner[j]; i < matrix1.compressed_format.inner[j + 1]; i++)
                            {
                                // row = row of matrix1 corresponding to the index i
                                size_t row = matrix1.compressed_format.outer[i];

                                // add the product of the non-zero elements to the "result" matrix
                                result(row, col) += matrix1.compressed_format.values[i] * matrix2.compressed_format.values[k];
                            }
                        }
                    }
                }
            }
            else
            {
                // iterate over rows of matrix1
                for (size_t row = 0; row < matrix1.rows; row++)
                {
                    // iterate over columns of matrix1 (and rows of matrix2) that are non-zero in the row "row" of matrix1
                    for (size_t j = matrix1.compressed_format.inner[row]; j < matrix1.compressed_format.inner[row + 1]; j++)
                    {
                        // k = column of matrix1 (or row of matrix2) that we are currently processing
                        size_t k = matrix1.compressed_format.outer[j];

                        // iterate over columns of matrix2 that are non-zero in the row k of matrix2
                        for (size_t i = matrix2.compressed_format.inner[k]; i < matrix2.compressed_format.inner[k + 1]; i++)
                        {
                            // col = column of matrix2 corresponding to the index i
                            size_t col = matrix2.compressed_format.outer[i];

                            // if we are on the diagonal of matrix2
                            if (col == k)
                            {
                                // add the product of the non-zero elements to the "result" matrix
                                result(row, col) += matrix1.compressed_format.values[j] * matrix2.compressed_format.values[i];
                            }
                        }
                    }
                }
            }
        }
        else
        {
            if (auto square_matrix = dynamic_cast<const SquareMatrix<T, S> *>(&m1))
            {
                auto &matrix1 = *square_matrix;
                if (matrix1.is_modified())
                {
                    throw std::invalid_argument("Matrix compression formats do not match");
                }
            }
            auto &matrix1 = m1;
            for (const auto &it1 : matrix1.uncompressed_format)
            {
                for (const auto &it2 : matrix2.uncompressed_format)
                {
                    if (it1.first.col == it2.first.row && it2.first.row == it2.first.col)
                    {
                        // add the product of the non-zero elements to the "result" vector
                        result(it1.first.row, it1.first.col) += it1.second * it2.second;
                    }
                }
            }
        }
        return result;
    };

    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const DiagonalView<T, S> &m1, const Matrix<T, S> &m2)
    {
        if (m1.get_cols() != m2.get_rows())
        {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }
        if (m1.is_compressed() != m2.is_compressed())
        {
            throw std::invalid_argument("Matrix compression formats do not match");
        }
        Matrix<T, S> result(m1.get_rows(), m2.get_cols());
        auto &matrix1 = m1.matrix;

        if (matrix1.is_modified())
        {
            if (auto square_matrix = dynamic_cast<const SquareMatrix<T, S> *>(&m2))
            {
                auto &matrix2 = *square_matrix;
                if (matrix2.is_modified())
                {
                    size_t cols = matrix2.get_cols();
                    for (size_t i = 0; i < cols; ++i)
                    {
                        result(i, i) = matrix1.compressed_format_mod.values[i] * matrix2.compressed_format_mod.values[i];
                    }
                }
                else
                {
                    throw std::invalid_argument("Matrix compression formats do not match");
                }
            }
            else
            {
                throw std::invalid_argument("Matrix compression formats do not match");
            }
        }
        else if (matrix1.is_compressed())
        {
            auto &matrix2 = m2;
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                // iterate over columns of matrix2
                for (size_t col = 0; col < matrix2.cols; col++)
                {
                    // iterate over rows of matrix2 (and columns of matrix1) that are non-zero in the column "col" of matrix2
                    for (size_t k = matrix2.compressed_format.inner[col]; k < matrix2.compressed_format.inner[col + 1]; k++)
                    {
                        // j = row of matrix2 (or column of matrix1) that we are currently processing
                        size_t j = matrix2.compressed_format.outer[k];

                        // iterate over rows of matrix1 that are non-zero in the column j of matrix1
                        for (size_t i = matrix1.compressed_format.inner[j]; i < matrix1.compressed_format.inner[j + 1]; i++)
                        {
                            // row = row of matrix1 corresponding to the index i
                            size_t row = matrix1.compressed_format.outer[i];

                            // if we are on the diagonal of matrix1
                            if (row == j)
                            {
                                // add the product of the non-zero elements to the "result" matrix
                                result(row, col) += matrix1.compressed_format.values[i] * matrix2.compressed_format.values[k];
                            }
                        }
                    }
                }
            }
            else
            {
                // iterate over rows of matrix1
                for (size_t row = 0; row < matrix1.rows; row++)
                {
                    // iterate over columns of matrix1 (and rows of matrix2) that are non-zero in the row "row" of matrix1
                    for (size_t j = matrix1.compressed_format.inner[row]; j < matrix1.compressed_format.inner[row + 1]; j++)
                    {
                        // k = column of matrix1 (or row of matrix2) that we are currently processing
                        size_t k = matrix1.compressed_format.outer[j];

                        // if we are the diagonal of matrix1
                        if (k == row)
                        {
                            // iterate over columns of matrix2 that are non-zero in the row k of matrix2
                            for (size_t i = matrix2.compressed_format.inner[k]; i < matrix2.compressed_format.inner[k + 1]; i++)
                            {
                                // col = column of matrix2 corresponding to the index i
                                size_t col = matrix2.compressed_format.outer[i];

                                // add the product of the non-zero elements to the "result" matrix
                                result(row, col) += matrix1.compressed_format.values[j] * matrix2.compressed_format.values[i];
                            }
                        }
                    }
                }
            }
        }
        else
        {
            if (auto square_matrix = dynamic_cast<const SquareMatrix<T, S> *>(&m2))
            {
                auto &matrix2 = *square_matrix;
                if (matrix2.is_modified())
                {
                    throw std::invalid_argument("Matrix compression formats do not match");
                }
            }
            auto &matrix2 = m2;
            for (const auto &it1 : matrix1.uncompressed_format)
            {
                for (const auto &it2 : matrix2.uncompressed_format)
                {
                    if (it1.first.col == it2.first.row && it1.first.row == it1.first.col)
                    {
                        // add the product of the non-zero elements to the "result" vector
                        result(it1.first.row, it1.first.col) += it1.second * it2.second;
                    }
                }
            }
        }
        return result;
    }
}

#endif // PRODUCTS_TPP