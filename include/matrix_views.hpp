#ifndef MATRIX_VIEWS_HPP
#define MATRIX_VIEWS_HPP

#include "matrix.hpp"
#include "proxy.hpp"

#include <execution>

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

        /// @brief call operator const version
        /// @param row row index
        /// @param col column index
        /// @return the matrix element at (row, col)
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

        /// @brief check if the matrix is in a compressed format
        /// @return true if the matrix is compressed, false otherwise
        bool is_compressed() const { return matrix.is_compressed(); };

        /// @brief get the number of rows
        /// @return number of rows
        size_t get_rows() const { return matrix.get_cols(); };
        
        /// @brief get the number of columns
        /// @return number of columns
        size_t get_cols() const { return matrix.get_rows(); };

        /// @brief get the number of non-zero elements
        /// @return number of non-zero elements
        virtual size_t get_nnz() const { return matrix.get_nnz(); };

        /// @brief calculate the norm of the matrix
        /// @tparam N type of the norm (One, Infinity, Frobenius)
        /// @return value of the norm
        template <NormType N>
        double norm() const{
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
        friend std::vector<U> operator*(const MatrixTransposeView<U, V> &m, const std::vector<U> &v);

        // multiply with another matrix
        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const MatrixTransposeView<U, V> &m1, const MatrixTransposeView<U, V> &m2);
        
        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const MatrixTransposeView<U, V> &m1, const Matrix<U, V> &m2);

        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const Matrix<U, V> &m1, const MatrixTransposeView<U, V> &m2);

        private:
        Matrix<T, S> &matrix;
    };

    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class MatrixDiagonalView 
    {
        public:
            /// @brief delete default constructor
            MatrixDiagonalView() = delete;

            /// @brief constructor
            /// @param matrix the matrix to see as diagonal: all off-diagonal elements are ignored
            MatrixDiagonalView(SquareMatrix<T, S> &matrix) : matrix(matrix) {}
        
            /// @brief call operator non-const version
            /// @param idx row and column index
            /// @return a proxy to the matrix element at (idx, idx)
            Proxy<T,S> operator()(size_t idx)
            {
                return matrix(idx, idx);
            }

            /// @brief call operator const version
            /// @param idx row and column index
            /// @return the matrix element at (idx, idx)
            T operator()(size_t idx) const
            {
                return matrix(idx, idx);
            }

            /// @brief  get the matrix
            /// @return the matrix
            SquareMatrix<T, S> &get_matrix() const
            {
                return matrix;
            }
            
            /// @brief check if the matrix is in a compressed format
            /// @return true if the matrix is compressed, false otherwise
            bool is_compressed() const { return matrix.is_compressed(); };

            /// @brief check if the matrix is in a modified compressed format
            /// @return true if the matrix is modified compressed, false otherwise 
            bool is_modified() const { return matrix.is_modified(); };

            /// @brief get the matrix dimension
            /// @return number of rows or columns (the matrix is square)
            size_t get_size() const { return matrix.get_rows(); };

            /// @brief get the number of non-zero elements in the diagonal
            /// @return number of non-zero elements in the diagonal
            virtual size_t get_nnz() const {
                T sum = 0;
                for(size_t i = 0; i < matrix.get_rows(); i++)
                {
                    sum += matrix(i,i);
                }
                return sum;
            };

            /// @brief calculate the norm of the matrix
            /// @tparam N type of the norm (One, Infinity, Frobenius)
            /// @return value of the norm
            template <NormType N>
            double norm() const {
                if constexpr (N = NormType::Frobenius){
                    double sum{0};
                    for(size_t i = 0; i < matrix.get_rows(); i++)
                    {
                        sum += std::abs(matrix(i,i)) * std::abs(matrix(i,i));
                    }
                    return std::sqrt(sum);
                }
                else{ // One or Infinity are equivalent for diagonal matrices
                    std::vector<double> diag(matrix.get_rows(), 0);
                    for(size_t i = 0; i < matrix.get_rows(); i++)
                    {
                        diag[i] = std::abs(matrix(i,i));
                    }
                    return *std::max_element(std::execution::par_unseq, diag.begin(), diag.end());
                }
            };

            // friend functions
            // multiply with a std::vector
            template <AddMulType U, StorageOrder V>
            friend std::vector<U> operator*(const MatrixDiagonalView<U, V> &m, const std::vector<U> &v);
            
            // multiply with another matrix
            template <AddMulType U, StorageOrder V>
            friend Matrix<U, V> operator*(const MatrixDiagonalView<U, V> &m1, const MatrixDiagonalView<U, V> &m2);

            template <AddMulType U, StorageOrder V>
            friend Matrix<U, V> operator*(const Matrix<U, V> &m1, const MatrixDiagonalView<U, V> &m2);

            template <AddMulType U, StorageOrder V>
            friend Matrix<U, V> operator*(const MatrixDiagonalView<U, V> &m1, const Matrix<U, V> &m2);

        private:
            SquareMatrix<T, S> &matrix;
    };

    //FRIENDS: MULTIPLICATIONS WITH TRANSPOSED VIEWS
    template <AddMulType T, StorageOrder S>
    std::vector<T> operator*(const MatrixTransposeView<T, S> &m, const std::vector<T> &v)
    {

        if (m.rows != v.size())
        {
            throw std::invalid_argument("Matrix and vector dimensions do not match for multiplication");
        }
        std::vector<T> result(m.cols, 0);
        if (not m.is_compressed())
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

    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const MatrixTransposeView<T, S> &m1, const Matrix<T, S> &m2)
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
    Matrix<T, S> operator*(const Matrix<T, S> &m1, const MatrixTransposeView<T, S> &m2)
    {
        return Matrix<T,S>();
    }


    //FRIENDS: MULTIPLICATIONS WITH DIAGONAL VIEWS
    template <AddMulType T, StorageOrder S>
    std::vector<T> operator*(const MatrixDiagonalView<T, S> &m, const std::vector<T> &v)
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
    Matrix<T, S> operator*(const MatrixDiagonalView<T, S> &m1, const MatrixDiagonalView<T, S> &m2)
    {
        if (m1.get_size() != m2.get_size())
        {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }
        if (m1.is_compressed() != m2.is_compressed())
        {
            throw std::invalid_argument("Matrix compression formats do not match");
        }

        Matrix<T, S> result(m1.get_size(), m1.get_size());

        return result;

    };

    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const Matrix<T, S> &m1, const MatrixDiagonalView<T, S> &m2)
    {

    };
    
    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const MatrixDiagonalView<T, S> &m1, const Matrix<T, S> &m2)  
    {

    };

}
#endif // MATRIX_VIEWS_HPP