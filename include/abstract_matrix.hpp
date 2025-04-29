#ifndef ABSTRACT_MATRIX_HPP
#define ABSTRACT_MATRIX_HPP

#include "storage.hpp"
#include "proxy.hpp"

#include <memory>

namespace algebra
{

    /// @brief type of norm
    enum class NormType
    {
        One,
        Infinity,
        Frobenius
    };

    // forward declarations
    // Matrix
    template <AddMulType T, StorageOrder S>
    class Matrix;
    // SquareMatrix
    template <AddMulType T, StorageOrder S>
    class SquareMatrix;
    // TransposeView
    template <AddMulType T, StorageOrder S>
    class TransposeView;
    // DiagonalView
    template <AddMulType T, StorageOrder S>
    class DiagonalView;


    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class AbstractMatrix
    {
    public:

        /// @brief clone method
        /// @return a pointer to the cloned object
        virtual std::unique_ptr<AbstractMatrix<T, S>> clone() const = 0;

        /// @brief default destructor
        virtual ~AbstractMatrix() = default;

        /// @brief set an element in the matrix (dynamic construction of the matrix)
        /// @param row row index
        /// @param col column index
        /// @param value value to set
        virtual void set(size_t row, size_t col, const T &value) = 0;

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

        /// @brief calculate the norm of the matrix
        /// @tparam N type of the norm (One, Infinity, Frobenius)
        /// @return value of the norm
        template <NormType N>
        double norm() const { 

            if (typeid(*this) == typeid(AbstractMatrix<T, S>))
            {
                throw std::invalid_argument("Cannot calculate norm of abstract matrix");
            }
            else if (typeid(*this) == typeid(SquareMatrix<T, S>))
            {
                return static_cast<const SquareMatrix<T, S> *>(this)->template norm<N>();
            }
            else if (typeid(*this) == typeid(DiagonalView<T, S>))
            {
                return static_cast<const DiagonalView<T, S> *>(this)->template norm<N>();
            }
            else if (typeid(*this) == typeid(TransposeView<T, S>))
            {
                return static_cast<const TransposeView<T, S> *>(this)->template norm<N>();
            }
            else if (typeid(*this) == typeid(Matrix<T, S>))
            {
                return static_cast<const Matrix<T, S> *>(this)->template norm<N>();
            }
            else
            {
                throw std::invalid_argument("Cannot calculate norm of unknown matrix type");
            }
         };

        /// @brief check if the matrix is in a compressed format
        virtual bool is_compressed() const = 0;

        /// @brief compress the matrix if it is in an uncompressed format
        virtual void compress() = 0;

        /// @brief uncompress the matrix if it is in a compressed format
        virtual void uncompress() = 0;

        /// @brief Function to read a matrix in Matrix Market format
        /// @param filename input file name
        virtual void reader(const std::string &filename) = 0;

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