#ifndef MATRIX_VIEWS_HPP
#define MATRIX_VIEWS_HPP

#include "matrix.hpp"
#include "square_matrix.hpp"
#include "proxy.hpp"
#include "abstract_matrix.hpp"

#include <execution>

namespace algebra
{
    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class TransposeView final : public AbstractMatrix<T, S>
    {
    public:
        // reference to the matrix
        Matrix<T, S> &matrix;

        /// @brief delete default constructor
        TransposeView() = delete;

        /// @brief constructor
        /// @param matrix the matrix to transpose
        TransposeView(Matrix<T, S> &matrix) : matrix(matrix) {};

        /// @brief virtual destructor
        virtual ~TransposeView() = default;

        /// @brief set an element in the matrix (dynamic construction of the matrix)
        /// @param row row index
        /// @param col column index
        /// @param value value to set
        void set(size_t row, size_t col, const T &value)
        {
            matrix.set(col, row, value);
        };

        /// @brief call operator non-const version
        /// @param row row index
        /// @param col column index
        /// @return a proxy to the matrix element at (row, col)
        Proxy<T, S> operator()(size_t row, size_t col)
        {
            return matrix(col, row);
        };

        /// @brief call operator const version
        /// @param row row index
        /// @param col column index
        /// @return the matrix element at (row, col)
        T operator()(size_t row, size_t col) const
        {
            return matrix(col, row);
        };

        /// @brief get the number of rows
        /// @return number of rows
        size_t get_rows() const { return matrix.get_cols(); };

        /// @brief get the number of columns
        /// @return number of columns
        size_t get_cols() const { return matrix.get_rows(); };

        /// @brief get the number of non-zero elements
        /// @return number of non-zero elements
        size_t get_nnz() const { return matrix.get_nnz(); };

        /// @brief calculate the norm of the matrix
        /// @tparam N type of the norm (One, Infinity, Frobenius)
        /// @return value of the norm
        template <NormType N>
        double norm() const
        {
            if constexpr (N == NormType::One)
            {
                return matrix.template norm<NormType::Infinity>();
            }
            else if constexpr (N == NormType::Infinity)
            {
                return matrix.template norm<NormType::One>();
            }
            else
            {
                return matrix.template norm<N>();
            }
        };

        // friend functions
        // multiply with a std::vector
        template <AddMulType U, StorageOrder V>
        friend std::vector<U> operator*(const TransposeView<U, V> &m, const std::vector<U> &v);

        // multiply with another matrix
        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const TransposeView<U, V> &m1, const TransposeView<U, V> &m2);

        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const TransposeView<U, V> &m1, const Matrix<U, V> &m2);

        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const Matrix<U, V> &m1, const TransposeView<U, V> &m2);
    };

    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class DiagonalView final : public AbstractMatrix<T, S>
    {
    public:
        // reference to the matrix
        SquareMatrix<T, S> &matrix;

        /// @brief delete default constructor
        DiagonalView() = delete;

        /// @brief constructor
        /// @param matrix the matrix to see as diagonal: all off-diagonal elements are ignored
        DiagonalView(SquareMatrix<T, S> &matrix) : matrix(matrix) {}

        /// @brief virtual destructor
        virtual ~DiagonalView() = default;

        /// @brief set an element in the matrix (dynamic construction of the matrix)
        /// @param row row index
        /// @param col column index
        /// @param value value to set
        void set(size_t row, size_t col, const T &value)
        {
            if (row == col)
            {
                matrix.set(row, col, value);
            }
            else
            {
                throw std::invalid_argument("Cannot set off-diagonal elements in a diagonal matrix");
            }
        }

        /// @brief delete call operator non-const version
        /// @param row row index
        /// @param col column index
        /// @return a proxy to the matrix element at (row, col)
        /// @note this function is used to set the diagonal element at (row, col)
        /// @note if the element is not on the diagonal, an exception is thrown
        Proxy<T, S> operator()(size_t row, size_t col){
            if (row == col)
            {
                return matrix(row, col);
            }
            else
            {
                throw std::invalid_argument("Cannot get off-diagonal elements in a diagonal matrix");
            }
        }

        /// @brief call operator const version
        /// @param row row index
        /// @param col column index
        /// @return the matrix element at (row, col)
        /// @note this function is used to get the diagonal element at (row, col)
        /// @note if the element is not on the diagonal, an exception is thrown
        T operator()(size_t row, size_t col) const {
            if (row == col)
            {
                return matrix(row, col);
            }
            else
            {
                throw std::invalid_argument("Cannot get off-diagonal elements in a diagonal matrix");
            }
        };

        /// @brief get the number of rows
        /// @return number of rows
        size_t get_rows() const { return matrix.get_rows(); };

        /// @brief get the number of columns
        /// @return number of columns
        size_t get_cols() const { return matrix.get_cols(); };

        /// @brief get the number of non-zero elements in the diagonal
        /// @return number of non-zero elements in the diagonal
        size_t get_nnz() const
        {
            T sum = 0;
            for (size_t i = 0; i < matrix.get_rows(); i++)
            {
                sum += matrix(i, i);
            }
            return sum;
        };

        /// @brief calculate the norm of the matrix
        /// @tparam N type of the norm (One, Infinity, Frobenius)
        /// @return value of the norm
        template <NormType N>
        double norm() const
        {
            if constexpr (N = NormType::Frobenius)
            {
                double sum{0};
                for (size_t i = 0; i < matrix.get_rows(); i++)
                {
                    sum += std::abs(matrix(i, i)) * std::abs(matrix(i, i));
                }
                return std::sqrt(sum);
            }
            else
            { // One or Infinity are equivalent for diagonal matrices
                std::vector<double> diag(matrix.get_rows(), 0);
                for (size_t i = 0; i < matrix.get_rows(); i++)
                {
                    diag[i] = std::abs(matrix(i, i));
                }
                return *std::max_element(std::execution::par_unseq, diag.begin(), diag.end());
            }
        };

        // friend functions
        // multiply with a std::vector
        template <AddMulType U, StorageOrder V>
        friend std::vector<U> operator*(const DiagonalView<U, V> &m, const std::vector<U> &v);

        // multiply with another matrix
        template <AddMulType U, StorageOrder V>
        friend SquareMatrix<U, V> operator*(const DiagonalView<U, V> &m1, const DiagonalView<U, V> &m2);

        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const Matrix<U, V> &m1, const DiagonalView<U, V> &m2);

        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const DiagonalView<U, V> &m1, const Matrix<U, V> &m2);
    };

    template <AddMulType T, StorageOrder S>
    std::vector<T> operator*(const TransposeView<T, S> &m, const std::vector<T> &v)
    {
        auto matrix = m.matrix;
        if (matrix.rows != v.size())
        {
            throw std::invalid_argument("Matrix and vector dimensions do not match for multiplication");
        }
        std::vector<T> result(matrix.cols, 0);
        if (typeid(matrix) == typeid(SquareMatrix<T, S>))
        {
            auto this_square = static_cast<const SquareMatrix<T, S> *>(&matrix);
            if (this_square->is_modified())
            {
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    size_t col, start, end;

                    // off-diagonal elements
                    for (col = 0; col < m.cols - 1; ++col)
                    {
                        // iterate over non-zero elements in the column "col" of m
                        start = m.compressed_format_mod.bind[col];
                        end = m.compressed_format_mod.bind[col + 1];
                        for (size_t j = start; j < end; ++j)
                        {
                            // row index of the non-zero element
                            size_t row = m.compressed_format_mod.bind[j];
                            // add off-diagonal elements
                            result[col] += m.compressed_format_mod.values[j] * v[row];
                        }
                    }

                    // handle last column
                    start = m.compressed_format_mod.bind[col];
                    end = m.compressed_format_mod.values.size();
                    for (size_t j = start; j < end; ++j)
                    {
                        // row index of the non-zero element
                        size_t row = m.compressed_format_mod.bind[j];
                        // add off-diagonal elements
                        result[col] += m.compressed_format_mod.values[j] * v[row];
                    }

                    // add diagonal elements
                    for (size_t i = 0; i < m.cols; ++i)
                    {
                        result[i] += m.compressed_format_mod.values[i] * v[i];
                    }
                }
                else
                {
                    size_t row, start, end;

                    // off-diagonal elements
                    for (row = 0; row < m.rows - 1; ++row)
                    {
                        // iterate over non-zero elements in the row "row" of m
                        start = m.compressed_format_mod.bind[row];
                        end = m.compressed_format_mod.bind[row + 1];
                        for (size_t j = start; j < end; ++j)
                        {
                            // col index of the non-zero element
                            size_t col = m.compressed_format_mod.bind[j];
                            // add off-diagonal elements
                            result[col] += m.compressed_format_mod.values[j] * v[row];
                        }
                    }

                    // handle last row
                    start = m.compressed_format_mod.bind[row];
                    end = m.compressed_format_mod.values.size();
                    for (size_t j = start; j < end; ++j)
                    {
                        // col index of the non-zero element
                        size_t col = m.compressed_format_mod.bind[j];
                        // add off-diagonal elements
                        result[col] += m.compressed_format_mod.values[j] * v[row];
                    }

                    // add diagonal elements
                    for (size_t i = 0; i < m.cols; ++i)
                    {
                        result[i] += m.compressed_format_mod.values[i] * v[i];
                    }
                }
            }
        }
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
                for (size_t col = 0; col < m.cols; col++)
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
                for (size_t row = 0; row < matrix.rows; row++)
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

    /*
        //FRIENDS: MULTIPLICATIONS WITH TRANSPOSED VIEWS

        template <AddMulType T, StorageOrder S>
        Matrix<T, S> operator*(const TransposeView<T, S> &m1, const TransposeView<T, S> &m2)
        {
            TransposeView<T,S> result(m2.matrix * m1.matrix);
            return Matrix<T,S>(result);
        }

        template <AddMulType T, StorageOrder S>
        Matrix<T, S> operator*(const TransposeView<T, S> &m1, const Matrix<T, S> &m2)
        {
            // TO DO

            Matrix<T,S> result(m1.get_rows, m2.get_cols);

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
            return Matrix<T,S>();
        }


        //FRIENDS: MULTIPLICATIONS WITH DIAGONAL VIEWS
        template <AddMulType T, StorageOrder S>
        std::vector<T> operator*(const DiagonalView<T, S> &m, const std::vector<T> &v)
        {
            if (m.get_size() != v.size())
            {
                throw std::invalid_argument("Matrix and vector dimensions do not match for multiplication");
            }
            std::vector<T> result(m.get_size(), 0);

            if (m.is_modified()){
                for (size_t i = 0; i < m.get_size(); i++)
                {
                    result[i] = m.compressed_format_mod.values[i] * v[i];
                }
                return result;
            }
            else if (m.is_compressed()){
                for (size_t i = 0; i < m.get_size(); i++)
                {
                    // for this format the call operator is the best option, we are forced to go by each stored row/column
                    result[i] = m(i, i) * v[i];
                }
                return result;
            }
            else{
                for (size_t i = 0; i < m.get_size(); i++)
                {
                    result[i] = m.uncompressed_format[{i,i}] * v[i];
                }
                return result;
            }
        };

        template <AddMulType T, StorageOrder S>
        SquareMatrix<T, S> operator*(const DiagonalView<T, S> &m1, const DiagonalView<T, S> &m2)
        {
            if (m1.get_size() != m2.get_size())
            {
                throw std::invalid_argument("Matrix dimensions do not match for multiplication");
            }
            if ((m1.is_modified() != m2.is_modified()) || (m1.is_compressed() != m2.is_compressed()))
            {
                throw std::invalid_argument("Matrix compression formats do not match");
            }

            SquareMatrix<T, S> result(m1.get_size());

            if (m1.is_modified())


            return result;

        };

        template <AddMulType T, StorageOrder S>
        Matrix<T, S> operator*(const Matrix<T, S> &m1, const DiagonalView<T, S> &m2)
        {

        };

        template <AddMulType T, StorageOrder S>
        Matrix<T, S> operator*(const DiagonalView<T, S> &m1, const Matrix<T, S> &m2)
        {

        };
    */
}
#endif // MATRIX_VIEWS_HPP