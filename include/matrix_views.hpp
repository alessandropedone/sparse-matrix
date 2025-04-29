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

        /// @brief initializing constructor
        /// @param rows number of rows
        /// @param cols number of columns
        TransposeView(size_t rows, size_t cols) : matrix(*new Matrix<T, S>(rows, cols)) {};

        /// @brief constructor
        /// @param matrix the matrix to transpose
        TransposeView(Matrix<T, S> &matrix) : matrix(matrix) {};

        /// @brief virtual destructor
        virtual ~TransposeView() = default;

        /// @brief set an element in the matrix (dynamic construction of the matrix)
        /// @param row row index
        /// @param col column index
        /// @param value value to set
        void set(size_t row, size_t col, const T &value) override { matrix.set(col, row, value); };

        /// @brief call operator non-const version
        /// @param row row index
        /// @param col column index
        /// @return a proxy to the matrix element at (row, col)
        Proxy<T, S> operator()(size_t row, size_t col) override { return matrix(col, row); };

        /// @brief call operator const version
        /// @param row row index
        /// @param col column index
        /// @return the matrix element at (row, col)
        T operator()(size_t row, size_t col) const override { return matrix(col, row); };

        /// @brief check if the matrix is in a compressed format
        virtual bool is_compressed() const override { return matrix.is_compressed(); };

        /// @brief compress the matrix if it is in an uncompressed format
        virtual void compress() override { matrix.compress(); };

        /// @brief uncompress the matrix if it is in a compressed format
        virtual void uncompress() override { matrix.uncompress(); };

        /// @brief Function to read a matrix in Matrix Market format
        /// @param filename input file name
        /// @note this function is not implemented for the TransposeView class
        void reader(const std::string &filename) override { matrix.reader(filename); };

        /// @brief get the number of rows
        /// @return number of rows
        size_t get_rows() const override { return matrix.get_cols(); };

        /// @brief get the number of columns
        /// @return number of columns
        size_t get_cols() const override { return matrix.get_rows(); };

        /// @brief get the number of non-zero elements
        /// @return number of non-zero elements
        size_t get_nnz() const override { return matrix.get_nnz(); };

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
    };

    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class DiagonalView final : public AbstractMatrix<T, S>
    {
    public:
        // reference to the matrix
        SquareMatrix<T, S> &matrix;

        /// @brief delete default constructor
        DiagonalView() = delete;

        /// @brief initializing constructor
        /// @param rows number of rows
        /// @param cols number of columns
        DiagonalView(size_t rows, size_t cols) : matrix(*new SquareMatrix<T, S>(rows, cols))
        {
            if (rows != cols)
            {
                throw std::invalid_argument("Matrix is not square");
            }
        };

        /// @brief constructor
        /// @param matrix the matrix to see as diagonal: all off-diagonal elements are ignored
        DiagonalView(SquareMatrix<T, S> &matrix) : matrix(matrix) {}

        /// @brief virtual destructor
        virtual ~DiagonalView() = default;

        /// @brief set an element in the matrix (dynamic construction of the matrix)
        /// @param row row index
        /// @param col column index
        /// @param value value to set
        void set(size_t row, size_t col, const T &value) override
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
        Proxy<T, S> operator()(size_t row, size_t col) override
        {
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
        T operator()(size_t row, size_t col) const override
        {
            if (row == col)
            {
                return matrix(row, col);
            }
            else
            {
                throw std::invalid_argument("Cannot get off-diagonal elements in a diagonal matrix");
            }
        };

        /// @brief check if the matrix is in a compressed format
        /// @return true if the matrix is (modified) compressed, false otherwise
        virtual bool is_compressed() const override { return matrix.is_compressed(); };

        /// @brief check if the matrix is in a modified compressed format
        /// @return true if the matrix is (modified) compressed, false otherwise
        virtual bool is_modified() const { return matrix.is_modified(); };

        /// @brief compress the matrix if it is in an uncompressed format
        virtual void compress() override { matrix.compress(); };

        /// @brief uncompress the matrix if it is in a compressed format
        virtual void uncompress() override { matrix.uncompress(); };

        /// @brief Function to read a matrix in Matrix Market format
        /// @param filename input file name
        virtual void reader(const std::string &filename) override { matrix.reader(filename); };

        /// @brief get the number of rows
        /// @return number of rows
        size_t get_rows() const override { return matrix.get_rows(); };

        /// @brief get the number of columns
        /// @return number of columns
        size_t get_cols() const override { return matrix.get_cols(); };

        /// @brief get the number of non-zero elements in the diagonal
        /// @return number of non-zero elements in the diagonal
        size_t get_nnz() const override
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
            if constexpr (N == NormType::Frobenius)
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
    };
}
#endif // MATRIX_VIEWS_HPP