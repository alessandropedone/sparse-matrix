/**
 * @file matrix_views.hpp
 * @brief Provides view classes for matrix operations such as transpose and diagonal extraction.
 *
 * This header defines the `TransposeView` and `DiagonalView` classes, which are lightweight wrappers
 * around existing matrix objects to provide alternative views (such as the transpose or diagonal) without
 * copying or modifying the underlying data. These views inherit from `AbstractMatrix` and can be used
 * interchangeably with other matrix types in the algebra namespace.
 *
 * - `TransposeView` allows access to the transpose of a matrix, redirecting element access and modification
 *   to the corresponding transposed indices of the original matrix.
 * - `DiagonalView` provides a view of only the diagonal elements of a square matrix, treating all off-diagonal
 *   elements as zero and preventing their modification.
 *
 * Both classes support cloning, norm calculation, and compressed storage operations, and are compatible
 * with the matrix interface defined in `AbstractMatrix`.
 *
 * @author
 * @date
 */
#ifndef MATRIX_VIEWS_HPP
#define MATRIX_VIEWS_HPP

#include "matrix.hpp"
#include "square_matrix.hpp"
#include "proxy.hpp"
#include "abstract_matrix.hpp"

#include <execution>

namespace algebra
{
    /**
     * @brief A view that represents the transpose of a given matrix.
     *
     * The TransposeView class provides a transposed view of a matrix without copying its data.
     * All operations are redirected to the underlying matrix with row and column indices swapped.
     * This class supports both general and square matrices, and transparently handles compressed formats.
     *
     * @tparam T The type of the matrix elements.
     * @tparam S The storage type or additional matrix traits.
     *
     * @note The TransposeView does not own the underlying matrix unless constructed with dimensions,
     *       in which case it creates a new matrix.
     *
     * @see AbstractMatrix
     * @see Matrix
     * @see SquareMatrix
     */
    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class TransposeView final : public AbstractMatrix<T, S>
    {
    public:
        Matrix<T, S> &matrix; /// reference to the matrix

        /// @brief delete default constructor
        TransposeView() = delete;

        /// @brief initializing constructor
        /// @param rows number of rows
        /// @param cols number of columns
        TransposeView(size_t rows, size_t cols) : matrix(*new Matrix<T, S>(rows, cols)) {};

        /// @brief constructor
        /// @param matrix the matrix to transpose
        TransposeView(Matrix<T, S> &matrix) : matrix(matrix) {};

        /// @brief clone method
        /// @return a pointer to the cloned object
        virtual std::unique_ptr<AbstractMatrix<T, S>> clone() const override
        {
            if (typeid(matrix) == typeid(SquareMatrix<T, S>))
            {
                auto &square_matrix = static_cast<SquareMatrix<T, S> &>(matrix);
                auto cloned_matrix = square_matrix.clone();
                return std::make_unique<TransposeView<T, S>>(*static_cast<SquareMatrix<T, S> *>(cloned_matrix.release()));
            }
            else
            {
                auto &general_matrix = static_cast<Matrix<T, S> &>(matrix);
                auto cloned_matrix = general_matrix.clone();
                return std::make_unique<TransposeView<T, S>>(*static_cast<Matrix<T, S> *>(cloned_matrix.release()));
            }
        };

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
        T operator()(size_t row, size_t col) const override
        {
            if (typeid(matrix) == typeid(SquareMatrix<T, S>))
            {
                const auto &square_matrix = static_cast<const SquareMatrix<T, S> &>(matrix);
                return square_matrix(col, row);
            }
            else if (typeid(matrix) == typeid(Matrix<T, S>))
            {
                const auto &general_matrix = static_cast<const Matrix<T, S> &>(matrix);
                return general_matrix(col, row);
            }
            else
            {
                throw std::invalid_argument("Matrix type not supported");
            }
        };

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

    /**
     * @brief A view that represents the diagonal of a square matrix.
     *
     * The DiagonalView class provides a view of the diagonal elements of a square matrix,
     * treating all off-diagonal elements as zero. This class supports both general and square
     * matrices, and transparently handles compressed formats.
     *
     * @tparam T The type of the matrix elements.
     * @tparam S The storage type or additional matrix traits.
     *
     * @note The DiagonalView does not own the underlying matrix unless constructed with dimensions,
     *       in which case it creates a new matrix.
     *
     * @see AbstractMatrix
     * @see Matrix
     * @see SquareMatrix
     */
    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class DiagonalView final : public AbstractMatrix<T, S>
    {
    public:
        SquareMatrix<T, S> &matrix; /// reference to the matrix

        /// @brief delete default constructor
        DiagonalView() = delete;

        /// @brief initializing constructor
        /// @param rows number of rows
        /// @param cols number of columns
        DiagonalView(size_t rows, size_t cols) : matrix(*new SquareMatrix<T, S>(rows))
        {
            if (rows != cols)
            {
                throw std::invalid_argument("Matrix must be square");
            }
        };

        /// @brief constructor
        /// @param matrix the matrix to see as diagonal: all off-diagonal elements are ignored
        DiagonalView(SquareMatrix<T, S> &matrix) : matrix(matrix) {}

        /// @brief clone method
        /// @return a pointer to the cloned object
        virtual std::unique_ptr<AbstractMatrix<T, S>> clone() const override
        {
            auto cloned_matrix = matrix.clone();
            return std::make_unique<DiagonalView<T, S>>(*static_cast<SquareMatrix<T, S> *>(cloned_matrix.release()));
        };

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
                const auto temp = matrix;
                return temp(row, col);
            }
            else
            {
                return T(0);
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
            size_t sum{0};
            for (size_t i = 0; i < matrix.get_rows(); i++)
            {
                sum += std::abs(static_cast<T>(matrix(i, i))) > std::numeric_limits<AbsReturnType_t<T>>::epsilon();
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
                    sum += std::abs(static_cast<T>(matrix(i, i))) * std::abs(static_cast<T>(matrix(i, i)));
                }
                return std::sqrt(sum);
            }
            else
            { // One or Infinity are equivalent for diagonal matrices
                std::vector<double> diag(matrix.get_rows(), 0);
                for (size_t i = 0; i < matrix.get_rows(); i++)
                {
                    diag[i] = std::abs(static_cast<T>(matrix(i, i)));
                }
                return *std::max_element(std::execution::par_unseq, diag.begin(), diag.end());
            }
        };
    };
}
#endif // MATRIX_VIEWS_HPP