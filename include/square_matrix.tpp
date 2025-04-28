#ifndef SQUARE_MATRIX_TPP
#define SQUARE_MATRIX_TPP

// For more verbose error messages
#include <cstring> // for strerror
#include <cerrno>  // for errno
#include <cassert>
#include <execution>

#include "square_matrix.hpp"

namespace algebra
{
    template <AddMulType T, StorageOrder S>
    const size_t SquareMatrix<T, S>::get_mod_size() const
    {
        size_t size = 0;
        for (size_t i = 0; i < this->rows; ++i)
        {
            if ((*this)(i, i) == 0)
            {
                ++size;
            }
        }
        size += this->get_nnz();
        return size;
    };

    /// @brief compress the matrix in modified format
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::compress_mod()
    {
        if (modified)
            return;

        // clear the modified compressed matrix
        compressed_format_mod.values.clear();
        compressed_format_mod.bind.clear();

        // reserve space for modified compressed structure
        const size_t size = get_mod_size(); // nnz + extra space for possible zero diagonal elements
        compressed_format_mod.values.resize(size);
        compressed_format_mod.bind.resize(size);
        std::fill(std::execution::par_unseq, compressed_format_mod.values.begin(), compressed_format_mod.values.end(), 0);
        std::fill(std::execution::par_unseq, compressed_format_mod.bind.begin(), compressed_format_mod.bind.end(), 0);

        // to store correctly the pointers in the bind vector (that don't account for the diagonal elements)
        size_t off_diag_idx = 0;

        if (this->compressed)
        {
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                // iterate over columns of m
                for (size_t col_idx = 0; col_idx < this->cols; col_idx++)
                {
                    // set col pointer
                    compressed_format_mod.bind[col_idx] = off_diag_idx + this->rows;

                    // iterate over rows of m that are non-zero in the column "col" of m
                    size_t start = this->compressed_format.inner[col_idx];
                    size_t end = this->compressed_format.inner[col_idx + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        // get the row index of the non-zero element
                        size_t row_idx = this->compressed_format.outer[j];

                        // distinguish diagonal and off-diagonal elements
                        if (row_idx == col_idx)
                        { // diagonal element
                            compressed_format_mod.values[row_idx] = this->compressed_format.values[j];
                        }
                        else
                        { // off-diagonal element
                            compressed_format_mod.values[this->rows + off_diag_idx] = this->compressed_format.values[j];
                            compressed_format_mod.bind[this->rows + off_diag_idx] = row_idx;
                            ++off_diag_idx;
                        }
                    }
                }
            }
            else
            {
                // iterate over rows of m
                for (size_t row_idx = 0; row_idx < this->rows; row_idx++)
                {
                    // set row pointer
                    compressed_format_mod.bind[row_idx] = off_diag_idx + this->rows;

                    // iterate over columns of m that are non-zero in the row "row" of m
                    size_t start = this->compressed_format.inner[row_idx];
                    size_t end = this->compressed_format.inner[row_idx + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        // get the col index of the non-zero element
                        size_t col_idx = this->compressed_format.outer[j];

                        // distinguish diagonal and off-diagonal elements
                        if (row_idx == col_idx)
                        { // diagonal element
                            compressed_format_mod.values[row_idx] = this->compressed_format.values[j];
                        }
                        else
                        { // off-diagonal element
                            compressed_format_mod.values[this->rows + off_diag_idx] = this->compressed_format.values[j];
                            compressed_format_mod.bind[this->rows + off_diag_idx] = col_idx;
                            ++off_diag_idx;
                        }
                    }
                }
            }

            // clear the compressed matrix
            this->compressed_format.inner.clear();
            this->compressed_format.outer.clear();
            this->compressed_format.values.clear();
        }
        else
        {
            if constexpr (S == StorageOrder::ColumnMajor)
            {

                for (const auto &it : this->uncompressed_format)
                {
                    if (it.first.col == it.first.row)
                    { // diagonal element
                        compressed_format_mod.values[it.first.row] = it.second;
                    }
                    else
                    { // off-diagonal element
                        compressed_format_mod.values[this->rows + off_diag_idx] = it.second;
                        compressed_format_mod.bind[this->rows + off_diag_idx] = it.first.row;
                        ++off_diag_idx;
                        // set column pointers
                        if (it.first.col < this->cols - 1)
                            ++compressed_format_mod.bind[it.first.col + 1];
                    }
                }
            }
            else
            {
                for (const auto &it : this->uncompressed_format)
                {
                    // std::cout << "\nInserting element-> Row: " << it.first.row << ", Col: " << it.first.col << ", Value: " << it.second << std::endl;
                    if (it.first.col == it.first.row)
                    { // diagonal element
                        compressed_format_mod.values[it.first.row] = it.second;
                    }
                    else
                    { // off-diagonal element
                        // std::cout << "Setting values vector[" << this->rows + off_diag_idx << "] to " << it.second << std::endl;
                        compressed_format_mod.values[this->rows + off_diag_idx] = it.second;
                        // std::cout << "Setting bind vector[" << this->rows + off_diag_idx << "] to " << it.first.col << std::endl;
                        compressed_format_mod.bind[this->rows + off_diag_idx] = it.first.col;
                        ++off_diag_idx;
                        if (it.first.row < this->rows - 1)
                            // set row pointers
                            ++compressed_format_mod.bind[it.first.row + 1];
                        // std::cout << "Updating row pointer: " << it.first.row + 1 << " to " << compressed_format_mod.bind[it.first.row + 1] << std::endl;
                    }
                }
            }
            // std::cout << "Row pointers vector: " << std::endl;
            for (size_t i = 1; i < this->rows; ++i)
            {
                compressed_format_mod.bind[i] += compressed_format_mod.bind[i - 1];
                compressed_format_mod.bind[i - 1] += this->rows; // add shift of diagonal elements
                // std::cout << compressed_format_mod.bind[i - 1] << "\t";
            }
            compressed_format_mod.bind[this->rows - 1] += this->rows;
            // std::cout << compressed_format_mod.bind[this->rows - 1] << std::endl;
            //  clear the uncompressed matrix
            this->uncompressed_format.clear();
        }

        // update flags
        this->compressed = false;
        this->modified = true;
        return;
    };

    /// @brief set an element in the matrix (dynamic construction of the matrix)
    /// @param row row index
    /// @param col column index
    /// @param value value to set
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::set(size_t row, size_t col, const T &value)
    {
        if (modified)
        {
            std::cout << "Matrix is in modified compressed format, uncompressing..." << std::endl;
            // uncompress the matrix
            uncompress();
        }
        Matrix<T, S>::set(row, col, value);
    };

    /// @brief compress the matrix if it is in an uncompressed format
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::compress()
    {
        if (this->compressed)
            return;
        if (modified)
        {
            // clear the compressed matrix
            this->compressed_format.inner.clear();
            this->compressed_format.outer.clear();
            this->compressed_format.values.clear();

            // reserve space for the compressed matrix
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                this->compressed_format.inner.resize(this->cols + 1);
            }
            else
            {
                this->compressed_format.inner.resize(this->rows + 1);
            }
            std::fill(std::execution::par_unseq, this->compressed_format.inner.begin(), this->compressed_format.inner.end(), 0);

            // fill the compressed matrix
            size_t index = 0; // keeps track of nnz elements
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                for (size_t i = 0; i < this->cols; ++i)
                {
                    // std::cout << "Col: " << i << std::endl;
                    bool flag = 0; // flag to check if the diagonal element has been inserted
                    size_t start = compressed_format_mod.bind[i];
                    // std::cout << "Col pointer start: " << start << std::endl;
                    size_t end = compressed_format_mod.bind[i + 1];
                    // handle last column
                    if (i + 1 == this->cols)
                    {
                        // std::cout << "Last column\n" << std::endl;
                        end = compressed_format_mod.values.size();
                    }
                    // std::cout << "Col pointer end: " << end << std::endl;
                    for (size_t j = start; j < end; ++j)
                    {
                        size_t rowidx = compressed_format_mod.bind[j];
                        // std::cout << "Row index: " << rowidx << std::endl;

                        // Insert diagonal element if it is not already inserted, !0 and I am after the diagonal
                        if (not flag and compressed_format_mod.values[i] != 0 and rowidx > i)
                        {
                            // std::cout << "Adding diagonal element: " << compressed_format_mod.values[i] << std::endl;
                            this->compressed_format.values.push_back(compressed_format_mod.values[i]);
                            this->compressed_format.outer.push_back(i);
                            ++index;
                            flag = true;
                        }

                        // Add off-diagonal elements
                        this->compressed_format.values.push_back(compressed_format_mod.values[j]);
                        this->compressed_format.outer.push_back(rowidx);
                        // std::cout << "Adding off-diagonal element after diagonal: " << compressed_format_mod.values[j] << std::endl;
                    }
                    // if there are no other elements in the col or no elements after the diagonal one, add the diagonal element
                    if (not flag and compressed_format_mod.values[i] != 0)
                    {
                        // std::cout << "Adding diagonal element: " << compressed_format_mod.values[i] << std::endl;
                        this->compressed_format.values.push_back(compressed_format_mod.values[i]);
                        this->compressed_format.outer.push_back(i);
                        ++index;
                        // updating the flag is useless here
                    }
                    index += end - start;
                    // std::cout << "Index of nnz elements for now: " << index << std::endl;
                    this->compressed_format.inner[i + 1] = index;
                }
            }
            else
            { // ROWORDER
                for (size_t i = 0; i < this->rows; ++i)
                {
                    bool flag = 0; // flag to check if the diagonal element has been inserted
                    size_t start = compressed_format_mod.bind[i];
                    size_t end = compressed_format_mod.bind[i + 1];
                    // handle last row
                    if (i + 1 == this->rows)
                    {
                        end = compressed_format_mod.values.size();
                    }
                    for (size_t j = start; j < end; ++j)
                    {
                        size_t colidx = compressed_format_mod.bind[j];

                        // Insert diagonal element if it is not already inserted, !0 and I am after the diagonal
                        if (not flag and compressed_format_mod.values[i] != 0 and colidx > i)
                        {
                            this->compressed_format.values.push_back(compressed_format_mod.values[i]);
                            this->compressed_format.outer.push_back(i);
                            ++index;
                            flag = true;
                        }

                        // Add off-diagonal elements
                        this->compressed_format.values.push_back(compressed_format_mod.values[j]);
                        this->compressed_format.outer.push_back(colidx);
                    }
                    // if there are no other elements in the row, add the diagonal element
                    if (not flag and compressed_format_mod.values[i] != 0)
                    {
                        this->compressed_format.values.push_back(compressed_format_mod.values[i]);
                        this->compressed_format.outer.push_back(i);
                        ++index;
                    }
                    index += end - start;
                    this->compressed_format.inner[i + 1] = index;
                }
            }

            // clear the modified compressed matrix
            compressed_format_mod.values.clear();
            compressed_format_mod.bind.clear();

            // update the flags
            this->modified = false;
            this->compressed = true;
            return;
        }
        Matrix<T, S>::compress();
        return;
    };

    /// @brief uncompress the matrix if it is in a compressed format
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::uncompress()
    {
        if (modified)
        {
            // clear the uncompressed format
            this->uncompressed_format.clear();

            // fill the uncompressed matrix
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                for (size_t col_idx = 0; col_idx < this->cols; ++col_idx)
                {
                    // add diagonal element
                    if (compressed_format_mod.values[col_idx] != 0)
                    {
                        this->uncompressed_format[{col_idx, col_idx}] = compressed_format_mod.values[col_idx];
                    }
                    size_t start = compressed_format_mod.bind[col_idx];
                    size_t end = compressed_format_mod.bind[col_idx + 1];
                    if (col_idx + 1 == this->cols)
                        end = compressed_format_mod.values.size();
                    for (size_t j = start; j < end; ++j)
                    {
                        size_t row_idx = compressed_format_mod.bind[j];
                        this->uncompressed_format[{row_idx, col_idx}] = compressed_format_mod.values[j];
                    }
                }
            }
            else
            {
                for (size_t row_idx = 0; row_idx < this->rows; ++row_idx)
                {
                    // add diagonal element
                    if (compressed_format_mod.values[row_idx] != 0)
                    {
                        this->uncompressed_format[{row_idx, row_idx}] = compressed_format_mod.values[row_idx];
                    }
                    size_t start = compressed_format_mod.bind[row_idx];
                    size_t end = compressed_format_mod.bind[row_idx + 1];
                    if (row_idx + 1 == this->rows)
                        end = compressed_format_mod.values.size();
                    for (size_t j = start; j < end; ++j)
                    {
                        size_t col_idx = compressed_format_mod.bind[j];
                        this->uncompressed_format[{row_idx, col_idx}] = compressed_format_mod.values[j];
                    }
                }
            }
            modified = false;
            return;
        }
        Matrix<T, S>::uncompress();
        return;
    };

    /// @brief call operator() const version
    /// @param row row index
    /// @param col column index
    /// @return element at (row, col)
    template <AddMulType T, StorageOrder S>
    T SquareMatrix<T, S>::operator()(size_t row, size_t col) const
    {
        // check if the index is in range
        if (row >= this->rows or col >= this->cols)
        {
            throw std::out_of_range("Index out of range");
        }
        if (modified)
        {
            if (row == col)
            {
                return compressed_format_mod.values[row];
            }
            else
            {
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    size_t start = compressed_format_mod.bind[col];
                    size_t end = (col != this->cols - 1) ? compressed_format_mod.bind[col + 1] : compressed_format_mod.values.size();
                    for (size_t j = start; j < end; ++j)
                    {
                        if (compressed_format_mod.bind[j] == row)
                        {
                            return compressed_format_mod.values[j];
                        }
                    }
                }
                else
                {
                    size_t start = compressed_format_mod.bind[row];
                    size_t end = (row != this->rows - 1) ? compressed_format_mod.bind[row + 1] : compressed_format_mod.values.size();
                    for (size_t j = start; j < end; ++j)
                    {
                        if (compressed_format_mod.bind[j] == col)
                        {
                            return compressed_format_mod.values[j];
                        }
                    }
                }
                return T(0);
            }
        }
        else
            return Matrix<T, S>::operator()(row, col);
    };

    /// @brief call operator() non-const version
    /// @param row row index
    /// @param col column index
    /// @return reference to the element at (row, col) with proxy (to avoid setting zero values)
    template <AddMulType T, StorageOrder S>
    Proxy<T, S> SquareMatrix<T, S>::operator()(size_t row, size_t col)
    {
        if (row >= this->rows or col >= this->cols)
            throw std::out_of_range("Index out of range");

        if (modified or this->compressed)
        {
            std::cout << "Matrix is compressed, uncompressing..." << std::endl;
            uncompress();
        }
        return Proxy<T, S>{this->uncompressed_format, row, col};
    };

    /// @brief resize the matrix
    /// @param rows number of rows
    /// @param cols number of columns
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::resize_and_clear(size_t dim)
    {
        this->rows = dim;
        this->cols = dim;
        this->compressed = false;
        this->modified = false;
        this->uncompressed_format.clear();
        this->compressed_format.inner.clear();
        this->compressed_format.outer.clear();
        this->compressed_format.values.clear();
        this->compressed_format_mod.values.clear();
        this->compressed_format_mod.bind.clear();
    };

    /// @brief calculate the norm of the matrix
    /// @tparam N type of the norm (One, Infinity, Frobenius)
    /// @return value of the norm
    template <AddMulType T, StorageOrder S>
    template <NormType N>
    double SquareMatrix<T, S>::norm() const
    {
        if (modified)
        {
            if constexpr (N == NormType::One)
            {
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    std::vector<double> col_sums(this->cols, 0);
                    for (size_t i = 0; i < this->cols; ++i)
                    {
                        size_t start = compressed_format_mod.bind[i];
                        size_t end = (i != this->cols - 1) ? compressed_format_mod.bind[i + 1] : compressed_format_mod.values.size();
                        for (size_t j = start; j < end; ++j)
                        {
                            col_sums[i] += std::abs(compressed_format_mod.values[j]);
                        }
                        col_sums[i] += std::abs(compressed_format_mod.values[i]);
                    }
                    return *std::max_element(std::execution::par_unseq, col_sums.begin(), col_sums.end());
                }
                else
                {
                    std::vector<double> col_sums(this->cols, 0);
                    for (size_t i = 0; i < this->rows; ++i)
                    {
                        size_t start = compressed_format_mod.bind[i];
                        size_t end = (i != this->rows - 1) ? compressed_format_mod.bind[i + 1] : compressed_format_mod.values.size();
                        for (size_t j = start; j < end; ++j)
                        {
                            size_t col = compressed_format_mod.bind[j];
                            col_sums[col] += std::abs(compressed_format_mod.values[j]);
                        }
                        col_sums[i] += std::abs(compressed_format_mod.values[i]);
                    }
                    return *std::max_element(std::execution::par_unseq, col_sums.begin(), col_sums.end());
                }
            }
            else if constexpr (N == NormType::Infinity)
            {
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    std::vector<double> row_sums(this->rows, 0);
                    for (size_t i = 0; i < this->cols; ++i)
                    {
                        size_t start = compressed_format_mod.bind[i];
                        size_t end = (i != this->cols - 1) ? compressed_format_mod.bind[i + 1] : compressed_format_mod.values.size();
                        for (size_t j = start; j < end; ++j)
                        {
                            size_t row = compressed_format_mod.bind[j];
                            row_sums[row] += std::abs(compressed_format_mod.values[j]);
                        }
                        row_sums[i] += std::abs(compressed_format_mod.values[i]);
                    }
                    return *std::max_element(std::execution::par_unseq, row_sums.begin(), row_sums.end());
                }
                else
                {
                    std::vector<double> row_sums(this->rows, 0);
                    for (size_t i = 0; i < this->rows; ++i)
                    {
                        size_t start = compressed_format_mod.bind[i];
                        size_t end = (i != this->rows - 1) ? compressed_format_mod.bind[i + 1] : compressed_format_mod.values.size();
                        for (size_t j = start; j < end; ++j)
                        {
                            row_sums[i] += std::abs(compressed_format_mod.values[j]);
                        }
                        row_sums[i] += std::abs(compressed_format_mod.values[i]);
                    }
                    return *std::max_element(std::execution::par_unseq, row_sums.begin(), row_sums.end());
                }
            }
            else // Frobenius
            {
                double norm{0};
                for (const auto it : compressed_format_mod.values)
                {
                    norm += std::abs(it) * std::abs(it);
                }
                return std::sqrt(norm);
            }
        }
        else
        {
            return Matrix<T, S>::template norm<N>();
        }
    };

    /// @brief reader method for the modified compressed matrix
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::reader(const std::string &filename)
    {
        std::ifstream file(filename);
        if (not file.is_open())
        {
            // more verbose error message
            throw std::runtime_error("Unable to open file '" + filename + "': " + strerror(errno));
        }

        std::string line;
        // Skip Matrix Market header and comments (first lines starting with %% or %)
        while (std::getline(file, line))
        {
            if (line.substr(0, 2) == "%%" or line.substr(0, 1) == "%")
            {
                continue;
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

        if (row_read != col_read)
        {
            throw std::invalid_argument("Matrix is not square");
        }

        // Resize the matrix
        auto dim = row_read;
        resize_and_clear(dim);

        // Read matrix values
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            T value;
            size_t row, col;
            iss >> row >> col >> value;
            assert(row <= this->rows and col <= this->cols);
            // I traslate the row and column indices to 0-based format and set the element
            //  in the matrix
            set(row - 1, col - 1, value);
        }

        file.close();
    };

    template <AddMulType T, StorageOrder S>
    size_t SquareMatrix<T, S>::get_nnz() const
    {
        if (modified)
        {
            size_t size = compressed_format_mod.values.size(); // this doesn't account for zeros in the diagonal
            for (size_t i = 0; i < this->rows; ++i)
            {
                if (compressed_format_mod.values[i] == 0)
                {
                    --size;
                }
            }
            return size;
        }
        else
        {
            return Matrix<T, S>::get_nnz();
        }
    };

    /// @brief multiply with a std::vector
    /// @param m matrix
    /// @param v vector
    /// @return the result of the multiplication
    /// @note this function is a friend of the Matrix class, so it can access the private members
    template <AddMulType T, StorageOrder S>
    std::vector<T> operator*(const SquareMatrix<T, S> &m, const std::vector<T> &v)
    {
        if (m.modified)
        {
            if (v.size() != m.cols)
            {
                throw std::invalid_argument("Matrix and vector dimensions do not match");
            }
            std::vector<T> result(m.rows, 0);
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
                        result[row] += m.compressed_format_mod.values[j] * v[col];
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
                    result[row] += m.compressed_format_mod.values[j] * v[col];
                }

                // add diagonal elements
                for (size_t i = 0; i < m.rows; ++i)
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
                        result[row] += m.compressed_format_mod.values[j] * v[col];
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
                    result[row] += m.compressed_format_mod.values[j] * v[col];
                }

                // add diagonal elements
                for (size_t i = 0; i < m.rows; ++i)
                {
                    result[i] += m.compressed_format_mod.values[i] * v[i];
                }
            }
            return result;
        }
        return static_cast<const Matrix<T, S> &>(m) * v;
    };

    /// @brief multiply with another matrix
    /// @param m1 first matrix
    /// @param m2 second matrix
    /// @return the result of the multiplication
    /// @note this function is a friend of the Matrix class, so it can access the private members
    template <AddMulType T, StorageOrder S>
    SquareMatrix<T, S> operator*(const SquareMatrix<T, S> &m1, const SquareMatrix<T, S> &m2)
    {
        if (m1.modified or m2.modified)
        {
            if (not m1.modified or not m2.modified)
            {
                std::runtime_error("Matrix multiplication between compressed and uncompressed matrix is not supported");
            }
            if (m1.cols != m2.rows)
            {
                throw std::invalid_argument("Matrix dimensions do not match");
            }
            SquareMatrix<T, S> result(m1.rows);
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                size_t col;

                ////// ITERATE OVER COLUMNS OF m2 //////

                for (col = 0; col < m2.cols - 1; ++col)
                {
                    // ADD OFF-DIAGONAL ELEMENTS OF m2 AND m1
                    size_t start = m2.compressed_format_mod.bind[col];
                    size_t end = m2.compressed_format_mod.bind[col + 1];
                    size_t k;
                    // iterate over rows of m2 (and columns of m1) that are non-zero in the column "col" of m2
                    for (k = start; k < end - 1; ++k)
                    {
                        // j = row of m2 (or column of m1) that we are currently processing
                        size_t j = m2.compressed_format_mod.bind[k];

                        // iterate over rows of m1 that are non-zero in the column j of m1
                        size_t s = m1.compressed_format_mod.bind[j];
                        size_t e = m1.compressed_format_mod.bind[j + 1];
                        for (size_t i = s; i < e; ++i)
                        {
                            // row = row of m1 corresponding to the index i
                            size_t row = m1.compressed_format_mod.bind[i];

                            // add the product of the non-zero off-diagonal elements to the "result" matrix
                            result(row, col) += m1.compressed_format_mod.values[i] * m2.compressed_format_mod.values[k];
                        }
                    }
                    // handle last column of m1 (last row of m2)
                    size_t j = m2.compressed_format_mod.bind[k];
                    size_t s = m1.compressed_format_mod.bind[j];
                    size_t e = (j + 1 == m1.cols) ? m1.compressed_format_mod.values.size() : m1.compressed_format_mod.bind[j + 1];
                    for (size_t i = s; i < e; ++i)
                    {
                        // row = row of m1 corresponding to the index i
                        size_t row = m1.compressed_format_mod.bind[i];

                        // add the product of the non-zero off-diagonal elements to the "result" matrix
                        result(row, col) += m1.compressed_format_mod.values[i] * m2.compressed_format_mod.values[k];
                    }

                    // ADD DIAGONAL ELEMENTS OF m1
                    // iterate over rows of m2 (and columns of m1) that are non-zero in the column "col" of m2
                    for (k = start; k < end; ++k)
                    {
                        // j = row of m2 (or column of m1) that we are currently processing
                        size_t j = m2.compressed_format_mod.bind[k];

                        // add the product between the diagonal element of m1 and current non-zero element of m2 to the "result" matrix
                        result(j, col) += m1.compressed_format_mod.values[j] * m2.compressed_format_mod.values[k];
                    }
                }

                ////// HANDLE LAST COLUMN OF m2 //////

                // ADD OFF-DIAGONAL ELEMENTS OF m2 AND m1
                size_t start = m2.compressed_format_mod.bind[col];
                size_t end = m2.compressed_format_mod.values.size();
                size_t k;
                // iterate over rows of m2 (and columns of m1) that are non-zero in the column "col" of m2
                for (k = start; k < end - 1; ++k)
                {
                    // j = row of m2 (or column of m1) that we are currently processing
                    size_t j = m2.compressed_format_mod.bind[k];

                    // iterate over rows of m1 that are non-zero in the column j of m1
                    size_t s = m1.compressed_format_mod.bind[j];
                    size_t e = m1.compressed_format_mod.bind[j + 1];
                    for (size_t i = s; i < e; ++i)
                    {
                        // row = row of m1 corresponding to the index i
                        size_t row = m1.compressed_format_mod.bind[i];

                        // add the product of the non-zero off-diagonal elements to the "result" matrix
                        result(row, col) += m1.compressed_format_mod.values[i] * m2.compressed_format_mod.values[k];
                    }
                }
                // handle last column of m1 (last row of m2)
                size_t j = m2.compressed_format_mod.bind[k];
                size_t s = m1.compressed_format_mod.bind[j];
                size_t e = (j + 1 == m1.cols) ? m1.compressed_format_mod.values.size() : m1.compressed_format_mod.bind[j + 1];
                for (size_t i = s; i < e; ++i)
                {
                    // row = row of m1 corresponding to the index i
                    size_t row = m1.compressed_format_mod.bind[i];

                    // add the product of the non-zero off-diagonal elements to the "result" matrix
                    result(row, col) += m1.compressed_format_mod.values[i] * m2.compressed_format_mod.values[k];
                }

                // ADD DIAGONAL ELEMENTS OF m1
                // iterate over rows of m2 (and columns of m1) that are non-zero in the column "col" of m2
                for (size_t k = start; k < end; ++k)
                {
                    // j = row of m2 (or column of m1) that we are currently processing
                    size_t j = m2.compressed_format_mod.bind[k];

                    // add the product between the diagonal element of m1 and current non-zero element of m2 to the "result" matrix
                    result(j, col) += m1.compressed_format_mod.values[j] * m2.compressed_format_mod.values[k];
                }

                ////// ADD DIAGONAL ELEMENTS OF m2 //////

                for (col = 0; col < m2.cols - 1; ++col)
                {
                    // multiply each column of m1 with the diagonal element of m2
                    size_t s = m1.compressed_format_mod.bind[col];
                    size_t e = m1.compressed_format_mod.bind[col + 1];
                    // iterate over the non-zero elements of m1 that are in the column "col"
                    for (size_t i = s; i < e; ++i)
                    {
                        // row = row of m1 corresponding to the index i
                        size_t row = m1.compressed_format_mod.bind[i];

                        // add the product between the diagonal element of m2 and current non-zero element of m1 to the "result" matrix
                        result(row, col) += m1.compressed_format_mod.values[i] * m2.compressed_format_mod.values[col];
                    }
                }
                // HANDLE LAST COLUMN OF m1
                // multiply each column of m1 with the diagonal element of m2
                s = m1.compressed_format_mod.bind[col];
                e = m1.compressed_format_mod.values.size();
                // iterate over the non-zero elements of m1 that are in the column "col"
                for (size_t i = s; i < e; ++i)
                {
                    // row = row of m1 corresponding to the index i
                    size_t row = m1.compressed_format_mod.bind[i];

                    // add the product between the diagonal element of m2 and current non-zero element of m1 to the "result" matrix
                    result(row, col) += m1.compressed_format_mod.values[i] * m2.compressed_format_mod.values[col];
                }

                ////// ADD PRODUCTS OF DIAGONAL ELEMENTS //////

                // iterate over the diagonal elements of m1 and m2
                for (size_t i = 0; i < m1.rows; ++i)
                {
                    // add the product between the diagonal elements of m1 and m2 to the "result" matrix
                    result(i, i) += m1.compressed_format_mod.values[i] * m2.compressed_format_mod.values[i];
                }
            }
            else
            {
                size_t row;

                ////// ITERATE OVER COLUMNS OF m1 //////

                for (row = 0; row < m1.rows - 1; ++row)
                {
                    // ADD OFF-DIAGONAL ELEMENTS OF m1 AND m2
                    size_t start = m1.compressed_format_mod.bind[row];
                    size_t end = m1.compressed_format_mod.bind[row + 1];
                    size_t k;
                    // iterate over columns of m1 (and rows of m2) that are non-zero in the row "row" of m1
                    for (k = start; k < end - 1; ++k)
                    {
                        // j = column of m1 (or row of m2) that we are currently processing
                        size_t j = m1.compressed_format_mod.bind[k];

                        // iterate over columns of m2 that are non-zero in the row j of m2
                        size_t s = m2.compressed_format_mod.bind[j];
                        size_t e = m2.compressed_format_mod.bind[j + 1];
                        for (size_t i = s; i < e; ++i)
                        {
                            // col = column of m2 corresponding to the index i
                            size_t col = m2.compressed_format_mod.bind[i];

                            // add the product of the non-zero off-diagonal elements to the "result" matrix
                            result(row, col) += m1.compressed_format_mod.values[k] * m2.compressed_format_mod.values[i];
                        }
                    }
                    // handle last row of m2 (last column of m1)
                    size_t j = m1.compressed_format_mod.bind[k];
                    size_t s = m2.compressed_format_mod.bind[j];
                    size_t e = (j + 1 == m2.rows) ? m2.compressed_format_mod.values.size() : m2.compressed_format_mod.bind[j + 1];
                    for (size_t i = s; i < e; ++i)
                    {
                        // col = column of m2 corresponding to the index i
                        size_t col = m2.compressed_format_mod.bind[i];

                        // add the product of the non-zero off-diagonal elements to the "result" matrix
                        result(row, col) += m1.compressed_format_mod.values[k] * m2.compressed_format_mod.values[i];
                    }

                    // ADD DIAGONAL ELEMENTS OF m2
                    // iterate over columns of m1 (and rows of m2) that are non-zero in the row "row" of m1
                    for (k = start; k < end; ++k)
                    {
                        // j = column of m1 (or row of m2) that we are currently processing
                        size_t j = m1.compressed_format_mod.bind[k];

                        // add the product between the diagonal element of m2 and current non-zero element of m1 to the "result" matrix
                        result(row, j) += m1.compressed_format_mod.values[k] * m2.compressed_format_mod.values[j];
                    }
                }

                ///// HANDLE LAST ROW OF m1 /////

                // ADD OFF-DIAGONAL ELEMENTS OF m1 AND m2
                size_t start = m1.compressed_format_mod.bind[row];
                size_t end = m1.compressed_format_mod.values.size();
                size_t k;
                // iterate over columns of m1 (and rows of m2) that are non-zero in the row "row" of m1
                for (k = start; k < end - 1; ++k)
                {
                    // j = column of m1 (or row of m2) that we are currently processing
                    size_t j = m1.compressed_format_mod.bind[k];

                    // iterate over columns of m2 that are non-zero in the row j of m2
                    size_t s = m2.compressed_format_mod.bind[j];
                    size_t e = m2.compressed_format_mod.bind[j + 1];
                    for (size_t i = s; i < e; ++i)
                    {
                        // col = column of m2 corresponding to the index i
                        size_t col = m2.compressed_format_mod.bind[i];

                        // add the product of the non-zero off-diagonal elements to the "result" matrix
                        result(row, col) += m1.compressed_format_mod.values[k] * m2.compressed_format_mod.values[i];
                    }
                }
                // handle last row of m2 (last column of m1)
                size_t j = m1.compressed_format_mod.bind[k];
                size_t s = m2.compressed_format_mod.bind[j];
                size_t e = (j + 1 == m2.rows) ? m2.compressed_format_mod.values.size() : m2.compressed_format_mod.bind[j + 1];
                for (size_t i = s; i < e; ++i)
                {
                    // col = column of m2 corresponding to the index i
                    size_t col = m2.compressed_format_mod.bind[i];

                    // add the product of the non-zero off-diagonal elements to the "result" matrix
                    result(row, col) += m1.compressed_format_mod.values[k] * m2.compressed_format_mod.values[i];
                }

                // ADD DIAGONAL ELEMENTS OF m2
                // iterate over columns of m1 (and rows of m2) that are non-zero in the row "row" of m1
                for (size_t k = start; k < end; ++k)
                {
                    // j = column of m1 (or row of m2) that we are currently processing
                    size_t j = m1.compressed_format_mod.bind[k];

                    // add the product between the diagonal element of m2 and current non-zero element of m1 to the "result" matrix
                    result(row, j) += m1.compressed_format_mod.values[k] * m2.compressed_format_mod.values[j];
                }

                ////// ADD DIAGONAL ELEMENTS OF m1 //////

                for (row = 0; row < m1.rows - 1; ++row)
                {
                    // multiply each row of m2 with the diagonal element of m1
                    size_t s = m2.compressed_format_mod.bind[row];
                    size_t e = m2.compressed_format_mod.bind[row + 1];
                    // iterate over the non-zero elements of m2 that are in the row "row"
                    for (size_t i = s; i < e; ++i)
                    {
                        // col = column of m2 corresponding to the index i
                        size_t col = m2.compressed_format_mod.bind[i];

                        // add the product between the diagonal element of m1 and current non-zero element of m2 to the "result" matrix
                        result(row, col) += m1.compressed_format_mod.values[row] * m2.compressed_format_mod.values[i];
                    }
                }
                // HANDLE LAST ROW OF m2
                // multiply each row of m2 with the diagonal element of m1
                s = m2.compressed_format_mod.bind[row];
                e = m2.compressed_format_mod.values.size();
                // iterate over the non-zero elements of m2 that are in the row "row"
                for (size_t i = s; i < e; ++i)
                {
                    // col = column of m2 corresponding to the index i
                    size_t col = m2.compressed_format_mod.bind[i];

                    // add the product between the diagonal element of m1 and current non-zero element of m2 to the "result" matrix
                    result(row, col) += m1.compressed_format_mod.values[row] * m2.compressed_format_mod.values[i];
                }

                ////// ADD PRODUCTS OF DIAGONAL ELEMENTS //////

                // iterate over the diagonal elements of m1 and m2
                for (size_t i = 0; i < m1.rows; ++i)
                {
                    // add the product between the diagonal elements of m1 and m2 to the "result" matrix
                    result(i, i) += m1.compressed_format_mod.values[i] * m2.compressed_format_mod.values[i];
                }
            }
            return result;
        }
        return static_cast<const Matrix<T, S> &>(m1) * static_cast<const Matrix<T, S> &>(m2);
    };
};

#endif // SQUARE_MATRIX_TPP