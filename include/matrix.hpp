#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "storage.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <execution>
#include <cmath>

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
                            while (it.first.col > index)
                            {
                                index++;
                                compressed_format.inner[index] = compressed_format.outer.size();
                            }
                        }
                        compressed_format.outer.push_back(it.first.row);
                    }
                    else
                    {
                        if (it.first.row > index)
                        {
                            while (it.first.row > index)
                            {
                                index++;
                                compressed_format.inner[index] = compressed_format.outer.size();
                            }
                        }
                        compressed_format.outer.push_back(it.first.col);
                    }
                    compressed_format.values.push_back(it.second);
                }
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    while (cols > index)
                    {
                        index++;
                        compressed_format.inner[index] = compressed_format.outer.size();
                    }
                }
                else
                {
                    while (rows > index)
                    {
                        index++;
                        compressed_format.inner[index] = compressed_format.outer.size();
                    }
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
                std::cout << "Matrix is compressed, uncompressing..." << std::endl;
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
                std::cout << "Matrix is compressed, uncompressing..." << std::endl;
                uncompress();
            }
            uncompressed_format.clear();
            compressed_format.inner.clear();
            compressed_format.outer.clear();
            compressed_format.values.clear();
        }

        /// @brief calculate the norm of the matrix
        /// @tparam N type of the norm (One, Infinity, Frobenius)
        /// @return value of the norm
        template <NormType N>
        double norm() {
            // be careful with complex numbers
            // smart/efficient way to calculate the norm?
            if (!compressed)
            {
                std::cout << "Matrix is uncompressed, compressing..." << std::endl;
                compress();
            }
            if constexpr (N == NormType::One)
            {               
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    double norm = 0;
                    for (size_t col = 0; col < cols; col++)
                    {
                        double sum = 0;
                        size_t start = compressed_format.inner[col];
                        size_t end = compressed_format.inner[col + 1];
                        for (size_t j = start; j < end; j++)
                        {
                            sum += std::abs(compressed_format.values[j]);
                        }
                        norm = std::max(norm, sum);
                    }
                    return norm;
                }
                else
                {
                    std::vector<double> col_sums(cols, 0);
                    for (size_t row = 0; row < rows; row++) // exploit locality (cache)
                    {
                        size_t start = compressed_format.inner[row];
                        size_t end = compressed_format.inner[row + 1];
                        for (size_t j = start; j < end; j++)
                        {
                            size_t col = compressed_format.outer[j];
                            col_sums[col] += std::abs(compressed_format.values[j]);
                        }
                    }
                    return *std::max_element(std::execution::par_unseq, col_sums.begin(), col_sums.end());                   ;
                }
            }
            else if constexpr (N == NormType::Infinity)
            {
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    std::vector<double> row_sums(rows, 0);
                    for (size_t col = 0; col < cols; col++)
                    {
                        size_t start = compressed_format.inner[col];
                        size_t end = compressed_format.inner[col + 1];
                        for (size_t j = start; j < end; j++)
                        {
                            size_t row = compressed_format.outer[j];
                            row_sums[row] += std::abs(compressed_format.values[j]);
                        }
                    }
                    return *std::max_element(std::execution::par_unseq, row_sums.begin(), row_sums.end());
                }
                else
                {
                    double norm = 0;
                    for (size_t row = 0; row < rows; row++)
                    {
                        double sum = 0;
                        size_t start = compressed_format.inner[row];
                        size_t end = compressed_format.inner[row + 1];
                        for (size_t j = start; j < end; j++)
                        {
                            sum += std::abs(compressed_format.values[j]);
                        }
                        norm = std::max(norm, sum);
                    }
                    return norm;
                }
            }
            else
            {
                double norm = 0;
                for (const auto &it : compressed_format.values)
                {
                    norm += std::abs(it) * std::abs(it);
                }
                return std::sqrt(std::abs(norm));
            }
        };

        /// @brief Function to read a matrix in Matrix Market format
        /// @param filename input file name
        void reader(const std::string &filename)
        {
            std::ifstream file(filename);
            if (!file.is_open())
            {
                throw std::runtime_error("Unable to open file: " + filename);
            }

            std::string line;
            // Skip Matrix Market header and comments (first lines starting with %%)
            while (std::getline(file, line))
            {
                if (line.substr(0, 2) == "%%" or line.substr(0, 1) == "%")
                {
                    continue; // Skip comment lines
                }
                else
                {
                    break;
                }
            }

            // Read matrix dimensions (rows, columns) and number of non-zero elements
            std::istringstream sizes(line);
            size_t row_read, col_read, nnz;
            sizes >> row_read >> col_read >> nnz;
            // Resize the matrix
            resize_and_clear(row_read, col_read);

            // Read matrix values
            while (std::getline(file, line))
            {
                std::istringstream iss(line);
                T value;
                size_t row, col;
                iss >> row >> col >> value;
                // I traslate the row and column indices to 0-based format and set the element
                //  in the matrix
                set(row - 1, col - 1, value);
            }

            file.close();
        };

        // friend functions
        // multiply with a std::vector
        template <AddMulType U, StorageOrder V>
        friend std::vector<U> operator*(Matrix<U, V> &m, const std::vector<U> &v);
        // multiply with another matrix
        template <AddMulType U, StorageOrder V>
        friend Matrix<U, V> operator*(Matrix<U, V> &m1, Matrix<U, V> &m2);

    protected:
        size_t rows;             // number of rows
        size_t cols;             // number of columns
        bool compressed = false; // flag to check if the matrix is compressed

        // storage for the matrix
        // uncompressed matrix
        UncompressedStorage<T, S> uncompressed_format; // COO format
        // compressed matrix
        CompressedStorage<T> compressed_format; // CSR or CSC format
    };

    /// @brief multiply a matrix with a vector
    /// @tparam T type of the matrix elements
    /// @tparam S storage order of the matrix (RowMajor or ColumnMajor)
    /// @param m matrix
    /// @param v vector
    /// @return vector result
    /// @note the vector must have the same number of rows as the matrix
    /// @note the matrix must be compressed (if not, it will be compressed)
    template <AddMulType T, StorageOrder S>
    std::vector<T> operator*(Matrix<T, S> &m, const std::vector<T> &v)
    {
        if (!m.is_compressed())
        {
            std::cout << "Matrix is uncompressed, compressing..." << std::endl;
            m.compress();
        }
        std::vector<T> result(m.rows, 0);
        if constexpr (S == StorageOrder::ColumnMajor)
        {
            // iterate over columns of m
            for (size_t col = 0; col < m.col; col++)
            {
                // iterate over rows of m that are non-zero in the column "col" of m
                for (size_t j = m.compressed_format.inner[col]; j < m.compressed_format.inner[col + 1]; j++)
                {
                    // row = row of m that we are currently processing
                    size_t row = m.compressed_format.outer[j];

                    // add the product of the non-zero elements to the "result" vector
                    result[row] += m.compressed_format.values[j] * v[col];
                }
            }
        }
        else
        {
            // iterate over rows of m
            for (size_t row = 0; row < m.rows; row++)
            {
                // iterate over columns of m that are non-zero in the row "row" of m
                for (size_t j = m.compressed_format.inner[row]; j < m.compressed_format.inner[row + 1]; j++)
                {
                    // col = column of m that we are currently processing
                    size_t col = m.compressed_format.outer[j];

                    // add the product of the non-zero elements to the "result" vector
                    result[row] += m.compressed_format.values[j] * v[col];
                }
            }
        }
        return result;
    }

    /// @brief mutliply two matrices with the same storage order
    /// @tparam T type of the matrix elements
    /// @tparam S storage order of the matrix (RowMajor or ColumnMajor)
    /// @param m1 first matrix
    /// @param m2 second matrix
    /// @return matrix result
    /// @note the result matrix will be compressed
    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(Matrix<T, S> &m1, Matrix<T, S> &m2)
    {
        if (!m1.is_compressed())
        {
            std::cout << "First matrix is uncompressed, compressing..." << std::endl;
            m1.compress();
        }
        if (!m2.is_compressed())
        {
            std::cout << "Second matrix is uncompressed, compressing..." << std::endl;
            m2.compress();
        }
        if (m1.cols != m2.rows)
        {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }

        Matrix<T, S> result(m1.rows, m2.cols);
        if constexpr (S == StorageOrder::ColumnMajor)
        {
            // iterate over columns of m2
            for (size_t col = 0; col < m2.cols; col++)
            {
                // iterate over rows of m2 (and columns of m1) that are non-zero in the column "col" of m2
                for (size_t k = m2.compressed_format.inner[col]; k < m2.compressed_format.inner[col + 1]; k++)
                {
                    // j = row of m2 (or column of m1) that we are currently processing
                    size_t j = m2.compressed_format.outer[k];

                    // iterate over rows of m1 that are non-zero in the column j of m1
                    for (size_t i = m1.compressed_format.inner[j]; i < m1.compressed_format.inner[j + 1]; i++)
                    {
                        // row = row of m1 corresponding to the index i
                        size_t row = m1.compressed_format.outer[i];

                        // add the product of the non-zero elements to the "result" matrix
                        result(row, col) += m1.compressed_format.values[i] * m2.compressed_format.values[k];
                    }
                }
            }
        }
        else
        {
            // iterate over rows of m1
            for (size_t row = 0; row < m1.rows; row++)
            {
                // iterate over columns of m1 (and rows of m2) that are non-zero in the row "row" of m1
                for (size_t j = m1.compressed_format.inner[row]; j < m1.compressed_format.inner[row + 1]; j++)
                {
                    // k = column of m1 (or row of m2) that we are currently processing
                    size_t k = m1.compressed_format.outer[j];

                    // iterate over columns of m2 that are non-zero in the row k of m2
                    for (size_t i = m2.compressed_format.inner[k]; i < m2.compressed_format.inner[k + 1]; i++)
                    {
                        // col = column of m2 corresponding to the index i
                        size_t col = m2.compressed_format.outer[i];

                        // add the product of the non-zero elements to the "result" matrix
                        result(row, col) += m1.compressed_format.values[j] * m2.compressed_format.values[i];
                    }
                }
            }
        }
        // compress the result matrix
        result.compress();
        return result;
    }

}

#endif // MATRIX_HPP
