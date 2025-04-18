#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "storage.hpp"
#include <vector>

namespace algebra
{

    /// @brief type of norm
    enum class NormType
    {
        One,
        Infinity,
        Frobenius
    };

    /// @brief Matrix class
    /// @tparam T type of the matrix elements
    /// @tparam S storage order of the matrix (RowMajor or ColumnMajor)
    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class Matrix
    {
    public:
        // delete default constructor
        Matrix() = delete;

        /// @brief constructor with size
        /// @param rows number of rows
        /// @param cols number of columns
        Matrix(size_t rows, size_t cols)
        {
            this->rows = rows;
            this->cols = cols;
            this->compressed = false;
        };

        /// @brief set an element in the matrix (dynamic construction of the matrix)
        /// @param row row index
        /// @param col column index
        /// @param value value to set
        void set(size_t row, size_t col, const T &value)
        {
            // check if the index is out of range
            if (row >= rows || col >= cols)
            {
                throw std::out_of_range("Index out of range");
            }
            if (compressed)
            {
                uncompress();
            }
            uncompressed[std::make_pair(row, col)] = value;
        }

        /// @brief check if the matrix is in a compressed format
        /// @return true if the matrix is compressed, false otherwise
        bool is_compressed() const { return compressed; };

        /// @brief compress the matrix if it is in an uncompressed format
        void compress()
        {
            if (!compressed)
            {
                compressed.inner.clear();
                compressed.outer.clear();
                compressed.values.clear();

                // resize the vectors of the compressed format
                compressed_format.inner.resize(rows+1);
                compressed_format.outer.resize(cols+1);
                compressed_format.values.resize(rows * cols);

                for (const auto &it : uncompressed)
                {
                    compressed.inner.push_back(it.first.first);
                    compressed.outer.push_back(it.first.second);
                    compressed.values.push_back(it.second);
                }
                compressed.inner.shrink_to_fit();
                compressed.outer.shrink_to_fit();
                compressed.values.shrink_to_fit();
                compressed = CompressedStorage<T>(compressed.inner, compressed.outer, compressed.values);
                uncompressed.clear();
                compressed = true;
            }
        };

        /// @brief uncompress the matrix if it is in a compressed format
        void uncompress()
        {
            if (compressed)
            {
                compressed.inner.clear();
                compressed.outer.clear();
                compressed.values.clear();

                for (const auto &it : uncompressed)
                {
                    compressed.inner.push_back(it.first.first);
                    compressed.outer.push_back(it.first.second);
                    compressed.values.push_back(it.second);
                }
                compressed.inner.shrink_to_fit();
                compressed.outer.shrink_to_fit();
                compressed.values.shrink_to_fit();
                compressed = CompressedStorage<T>(compressed.inner, compressed.outer, compressed.values);
                compressed = false;
            }
        };

        /// @brief call operator() const version
        /// @param row row index
        /// @param col column index
        /// @return element at (row, col)
        T operator()(size_t row, size_t col) const
        {
            // check if the index is in range
            if (row >= rows || col >= cols)
            {
                throw std::out_of_range("Index out of range");
            }
            // check if the matrix is compressed
            if (!compressed)
            {
                auto it = uncompressed_format.find(std::make_pair(row, col));
                if (it != uncompressed_format.end())
                {
                    return it->second;
                }
                else
                {
                    return T(0);
                }
            }
            else
            {
                auto it = std::find(compressed.outer.begin(), compressed.outer.end(), col);
                if (it != compressed.outer.end())
                {
                    size_t index = std::distance(compressed.outer.begin(), it);
                    return compressed.values[index];
                }
                else
                {
                    return T(0);
                }
            }
        };

        /// @brief call operator() non-const version
        /// @param row row index
        /// @param col column index
        /// @return reference to the element at (row, col)
        T &operator()(size_t row, size_t col)
        {
            // check if the index is in range
            if (row >= rows || col >= cols)
            {
                throw std::out_of_range("Index out of range");
            }
            if (compressed)
            {
                uncompress();
            }
            return uncompressed_format[std::make_pair(row, col)];
        }

        /// @brief resize the matrix
        /// @param rows number of rows
        /// @param cols number of columns
        void resize_and_clear(size_t rows, size_t cols)
        {
            this->rows = rows;
            this->cols = cols;
            if (compressed)
            {
                uncompress();
            }
            uncompressed_format.clear();
            compressed_format.inner.clear();
            compressed_format.outer.clear();
            compressed_format.values.clear();
        }

        /// @brief calculate the norm of the matrix
        /// @tparam N type of the norm
        /// @return value of the norm
        template <NormType N>
        T norm() const {
            // be careful with complex numbers
            // smart/efficient way to calculate the norm?
        };

        // friend functions
        // multiply with a std::vector
        friend std::vector<T> operator*(const Matrix &m, const std::vector<T> &v);
        // multiply with another matrix
        friend Matrix operator*(const Matrix &m1, const Matrix &m2);

    private:
        size_t rows;             // number of rows
        size_t cols;             // number of columns
        bool compressed = false; // flag to check if the matrix is compressed

        // storage for the matrix
        // uncompressed matrix
        UncompressedStorage<T, S> uncompressed_format; // COO format
        // compressed matrix
        CompressedStorage<T> compressed_format; // CSR or CSC format
    };

    template <AddMulType T, StorageOrder S>
    std::vector<T> operator*(const Matrix<T, S> &m, const std::vector<T> &v)
    {
        // check if the matrix is compressed
        if (m.is_compressed())
        {
            // uncompress the matrix
            m.uncompress();
        }
        // multiply the matrix with the vector
        std::vector<T> result(m.rows());
        for (size_t i = 0; i < m.rows(); i++)
        {
            result[i] = 0;
            for (size_t j = 0; j < m.cols(); j++)
            {
                result[i] += m(i, j) * v[j];
            }
        }
        return result;
    }

    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const Matrix<T, S> &m1, const Matrix<T, S> &m2)
    {
        // TBD
    }

}

#endif // MATRIX_HPP
