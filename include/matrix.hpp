#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "storage.hpp"
#include <vector>
#include <iostream>
#include <fstream>
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
            uncompressed_format[{row, col}] = value;
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
                compressed_format.inner.clear();
                compressed_format.outer.clear();
                compressed_format.values.clear();

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
                for (const auto &it : uncompressed_format)
                {
                    if constexpr (S == StorageOrder::ColumnMajor)
                    {
                        if (it.first.col > index)
                        {
                            index = it.first.col;
                            compressed_format.inner[index] = compressed_format.outer.size();
                        }
                        compressed_format.outer.push_back(it.first.row);
                    }
                    else
                    {
                        if (it.first.row > index)
                        {
                            index = it.first.row;
                            compressed_format.inner[index] = compressed_format.outer.size();
                        }
                        compressed_format.outer.push_back(it.first.col);
                    }
                    compressed_format.values.push_back(it.second);
                }
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    compressed_format.inner[cols] = compressed_format.outer.size();
                }
                else
                {
                    compressed_format.inner[rows] = compressed_format.outer.size();
                }

                // clear the uncompressed matrix
                uncompressed_format.clear();

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
                uncompressed_format.clear();

                // fill the uncompressed matrix
                if constexpr (S == StorageOrder::ColumnMajor)
                {

                    for (size_t col_idx = 0; col_idx < cols; col_idx++)
                    {
                        size_t start = compressed_format.inner[col_idx];
                        size_t end = compressed_format.inner[col_idx + 1];
                        for (size_t j = start; j < end; j++)
                        {
                            size_t row_idx = compressed_format.outer[j];
                            uncompressed_format[{row_idx, col_idx}] = compressed_format.values[j];
                        }
                    }
                }
                else
                {
                    for (size_t row_idx = 0; row_idx < rows; row_idx++)
                    {
                        size_t start = compressed_format.inner[row_idx];
                        size_t end = compressed_format.inner[row_idx + 1];
                        for (size_t j = start; j < end; j++)
                        {
                            size_t col_idx = compressed_format.outer[j];
                            uncompressed_format[{row_idx, col_idx}] = compressed_format.values[j];
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

                    size_t start = compressed_format.inner[col];
                    size_t end = compressed_format.inner[col + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        if (compressed_format.outer[j] == row)
                        {
                            return compressed_format.values[j];
                        }
                    }
                }
                else
                {
                    size_t start = compressed_format.inner[row];
                    size_t end = compressed_format.inner[row + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        if (compressed_format.outer[j] == col)
                        {
                            return compressed_format.values[j];
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

        /// @brief Function to read a matrix in Matrix Market format
        /// @param filename input file name
        void reader(const std::string &filename)
        {
            std::ifstream file(filename);
            if (!file.is_open()) {
                throw std::runtime_error("Unable to open file: " + filename);
            }

            std::string line;
            // Skip Matrix Market header and comments (first lines starting with %%)
            while (std::getline(file, line)) {
                if (line.substr(0, 2) == "%%" or line.substr(0, 1) == "%") {
                    continue; // Skip comment lines
                } else {
                    break;
                }
            }

            // Read matrix dimensions (rows, columns) and number of non-zero elements
            std::istringstream sizes(line);
            size_t row_read, col_read, nnz;
            sizes >> row_read >> col_read >> nnz;
            // Resize the matrix
            rows = row_read;
            cols = col_read;

            // Read matrix values            
            while (std::getline(file, line)) {
                std::istringstream iss(line);
                T value;
                size_t row, col;
                iss >> row >> col >> value;
                //I traslate the row and column indices to 0-based format and set the element
                // in the matrix
                set(row-1, col-1, value);
            }
        
            file.close();

        };

        // friend functions
        // multiply with a std::vector
        template <AddMulType U, StorageOrder V>
        friend std::vector<U> operator*(const Matrix<U, V> &m, const std::vector<U> &v);
        // multiply with another matrix
        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(const Matrix<U, V> &m1, const Matrix<U, V> &m2);

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
