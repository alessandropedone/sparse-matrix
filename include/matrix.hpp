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
                std::cout << "Matrix is compressed, uncompressing..." << std::endl;
                uncompress();
            }
            uncompressed[{row, col}] = value;
        }

        /// @brief check if the matrix is in a compressed format
        /// @return true if the matrix is compressed, false otherwise
        bool is_compressed() const { return compressed; };

        /// @brief compress the matrix if it is in an uncompressed format
        void compress()
        {
            if (!compressed)
            {
                // clear the compressed matrix
                compressed.inner.clear();
                compressed.outer.clear();
                compressed.values.clear();

                // reserve space for the compressed matrix
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    compressed_format.inner.resize(cols + 1);
                }
                else
                {
                    compressed_format.inner.resize(rows + 1);
                }
                std::fill(compressed_format.inner.begin(), compressed_format.inner.end(), 0);

                // fill the compressed matrix
                size_t index = 0;
                for (const auto &it : uncompressed)
                {
                    if constexpr (S == StorageOrder::ColumnMajor)
                    {
                        if (it.first.col > index)
                        {
                            index = it.first.col;
                            compressed_format.inner[index] = compressed.outer.size();
                        }
                        compressed.outer.push_back(it.first.row);
                    }
                    else
                    {
                        if (it.first.row > index)
                        {
                            index = it.first.row;
                            compressed_format.inner[index] = compressed.outer.size();
                        }
                        compressed.outer.push_back(it.first.col);
                    }
                    compressed.values.push_back(it.second);
                }
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    compressed_format.inner[cols] = compressed.outer.size();
                }
                else
                {
                    compressed_format.inner[rows] = compressed.outer.size();
                }

                // clear the uncompressed matrix
                uncompressed.clear();

                // update the compressed flag
                compressed = true;
            }
        };

        /// @brief uncompress the matrix if it is in a compressed format
        void uncompress()
        {
            if (compressed)
            {
                // clear the uncompressed matrix
                uncompressed.clear();

                // fill the uncompressed matrix
                if constexpr (S == StorageOrder::ColumnMajor)
                {

                    for (size_t col_idx = 0; col_idx < cols; col_idx++)
                    {
                        size_t start = compressed.inner[col_idx];
                        size_t end = compressed.inner[col_idx + 1];
                        for (size_t j = start; j < end; j++)
                        {
                            size_t row_idx = compressed.outer[j];
                            uncompressed[{row_idx, col_idx}] = compressed.values[j];
                        }
                    }
                }
                else
                {
                    for (size_t row_idx = 0; row_idx < rows; row_idx++)
                    {
                        size_t start = compressed.inner[row_idx];
                        size_t end = compressed.inner[row_idx + 1];
                        for (size_t j = start; j < end; j++)
                        {
                            size_t col_idx = compressed.outer[j];
                            uncompressed[{row_idx, col_idx}] = compressed.values[j];
                        }
                    }
                }

                // clear the compressed matrix
                compressed_format.inner.clear();
                compressed_format.outer.clear();
                compressed_format.values.clear();

                // update the compressed flag
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
                auto it = uncompressed_format.find({row, col});
                if (it != uncompressed_format.end())
                {
                    return it->second;
                }
            }
            else
            {
                if constexpr (S == StorageOrder::ColumnMajor)
                {

                    size_t start = compressed.inner[col];
                    size_t end = compressed.inner[col + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        if (compressed.outer[j] == row)
                        {
                            return compressed.values[j];
                        }
                    }
                }
                else
                {
                    size_t start = compressed.inner[row];
                    size_t end = compressed.inner[row + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        if (compressed.outer[j] == col)
                        {
                            return compressed.values[j];
                        }
                    }
                }
            }
            return T(0);
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
            return uncompressed_format[{row, col}];
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
