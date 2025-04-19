#ifndef VIEWS_HPP
#define VIEWS_HPP

#include "matrix.hpp"

namespace algebra
{
    template <typename T, StorageOrder S = StorageOrder::RowMajor>
    class DiagonalView : public Matrix<T, S>
    {
    public:
        DiagonalView(Matrix<T, S> &matrix) : matrix(matrix) {}
        T &operator()(size_t row, size_t col)
        {
            return (row == col) ? matrix(row, col) : T(0);
        }
        T &operator()(size_t row, size_t col) const
        {
            return (row == col) ? matrix(row, col) : T(0);
        }

        // friend functions
        // multiply with a std::vector
        template <AddMulType U, StorageOrder V>
        friend std::vector<U> operator*(DiagonalView<U, V> &m, const std::vector<U> &v);

    private:
        Matrix<T,S> &matrix;
    };

    template <AddMulType T, StorageOrder V>
    std::vector<T> operator*(DiagonalView<T, V> &m, const std::vector<T> &v)
    {
        std::vector<T> result(m.rows, 0);
        for (size_t i = 0; i < m.rows; i++)
        {
            result[i] = m(i, i) * v[i];
        }
        return result;
    };

    template <typename T, StorageOrder S = StorageOrder::RowMajor>
    class TransposeView : public Matrix<T, S>
    {
    public:
        TransposeView(Matrix<T, S> &matrix) : matrix(matrix) {}
        T operator()(size_t row, size_t col)
        {
            return matrix(col, row);
        }
        T operator()(size_t row, size_t col) const
        {
            return matrix(col, row);
        }

    private:
        Matrix<T,S> &matrix;
    };
}

#endif // VIEWS_HPP