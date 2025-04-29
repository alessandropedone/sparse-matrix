#ifndef MATRIX_TPP
#define MATRIX_TPP

#include "matrix.hpp"
#include "square_matrix.hpp"

// For more verbose error messages
#include <cstring> // for strerror
#include <cerrno>  // for errno
#include <cassert>

namespace algebra
{
    /// @brief constructor from a TransposeView
    /// @note the constructed matrix is in uncompressed format
    /// @param view transposed view of matrix to copy
    template <AddMulType T, StorageOrder S>
    Matrix<T, S>::Matrix(const TransposeView<T, S> &view)
    {
        auto matrix = view.matrix;
        if (matrix.is_compressed())
        {
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                for (size_t col = 0; col < matrix.get_cols(); col++)
                {
                    size_t start = matrix.compressed_format.inner[col];
                    size_t end = matrix.compressed_format.inner[col + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        size_t row = matrix.compressed_format.outer[j];
                        T value = matrix.compressed_format.values[j];
                        this->set(col, row, value);
                    }
                }
            }
            else
            {
                for (size_t row = 0; row < matrix.get_rows(); row++)
                {
                    size_t start = matrix.compressed_format.inner[row];
                    size_t end = matrix.compressed_format.inner[row + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        size_t col = matrix.compressed_format.outer[j];
                        T value = matrix.compressed_format.values[j];
                        this->set(col, row, value);
                    }
                }
            }
            // set the number of rows and columns
            this->rows = matrix.get_cols();
            this->cols = matrix.get_rows();
            // set the compressed flag
            this->compressed = true;
        }
        else
        {
            for (const auto &it : matrix.uncompressed_format)
            {
                // set the element in the matrix
                this->set(it.first.col, it.first.row, it.second);
            }
            // set the number of rows and columns
            this->rows = matrix.get_cols();
            this->cols = matrix.get_rows();
            // set the compressed flag
            this->compressed = false;
        }
    }

    /// @brief constructor from a DiagonalView
    /// @note the constructed matrix is in uncompressed format
    /// @param view diagonal view of matrix to copy
    template <AddMulType T, StorageOrder S>
    Matrix<T, S>::Matrix(const DiagonalView<T, S> &view)
    {
        auto matrix = view.matrix;
        if (matrix.is_compressed())
        {
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                for (size_t col = 0; col < matrix.get_cols(); col++)
                {
                    size_t start = matrix.compressed_format.inner[col];
                    size_t end = matrix.compressed_format.inner[col + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        size_t row = matrix.compressed_format.outer[j];
                        if (row == col)
                        {
                            T value = matrix.compressed_format.values[j];
                            this->set(row, col, value);
                        }
                    }
                }
            }
            else
            {
                for (size_t row = 0; row < matrix.get_rows(); row++)
                {
                    size_t start = matrix.compressed_format.inner[row];
                    size_t end = matrix.compressed_format.inner[row + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        size_t col = matrix.compressed_format.outer[j];
                        if (row == col)
                        {
                            T value = matrix.compressed_format.values[j];
                            this->set(row, col, value);
                        }
                    }
                }
            }
            // set the number of rows and columns
            this->rows = matrix.get_rows();
            this->cols = matrix.get_cols();
            // set the compressed flag
            this->compressed = true;
        }
        else
        {
            for (const auto &it : matrix.uncompressed_format)
            {
                if (it.first.row == it.first.col)
                {
                    // set the element in the matrix
                    this->set(it.first.row, it.first.col, it.second);
                }
            }
            // set the number of rows and columns
            this->rows = matrix.get_rows();
            this->cols = matrix.get_cols();
            // set the compressed flag
            this->compressed = false;
        }
    }

    template <AddMulType T, StorageOrder S>
    void Matrix<T, S>::set(size_t row, size_t col, const T &value)
    {
        // check if the index is out of range
        if (row >= rows or col >= cols)
        {
            throw std::out_of_range("Index out of range");
        }
        if (compressed)
        {
            std::cout << "Matrix is compressed, uncompressing..." << std::endl;
            uncompress();
        }
        if (value != 0)
        {
            uncompressed_format[{row, col}] = value;
        }
        else
        {
            auto it = uncompressed_format.find({row, col});
            if (it != uncompressed_format.end())
            {
                uncompressed_format.erase(it);
            }
        }
    }

    template <AddMulType T, StorageOrder S>
    void Matrix<T, S>::compress()
    {
        if (compressed)
            return;

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
                // until the index reaches the value of the column index
                while (it.first.col > index)
                {
                    // increment the index
                    index++;

                    // set the inner index to the size of the outer vector
                    compressed_format.inner[index] = compressed_format.outer.size();
                }
                // add the row index to the outer vector
                compressed_format.outer.push_back(it.first.row);
            }
            else
            {
                // until the index reaches the value of the row index
                while (it.first.row > index)
                {
                    // increment the index
                    index++;

                    // set the inner index to the size of the outer vector
                    compressed_format.inner[index] = compressed_format.outer.size();
                }
                // add the column index to the outer vector
                compressed_format.outer.push_back(it.first.col);
            }
            // add the value to the values vector
            compressed_format.values.push_back(it.second);
        }
        if constexpr (S == StorageOrder::ColumnMajor)
        {
            // check if you have passed the column index
            while (cols > index)
            {
                // until the index reaches the value of the column index
                index++;

                // set the inner index to the size of the outer vector
                compressed_format.inner[index] = compressed_format.outer.size();
            }
        }
        else
        {
            // check if you have passed the row index
            while (rows > index)
            {
                // until the index reaches the value of the row index
                index++;

                // set the inner index to the size of the outer vector
                compressed_format.inner[index] = compressed_format.outer.size();
            }
        }

        // clear the uncompressed matrix
        uncompressed_format.clear();

        // update the compressed flag
        compressed = true;
    };

    template <AddMulType T, StorageOrder S>
    void Matrix<T, S>::compress_parallel()
    {
        if (compressed)
            return;

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
        std::fill(std::execution::par_unseq, compressed_format.inner.begin(), compressed_format.inner.end(), 0);
        const size_t nnz = get_nnz();
        compressed_format.outer.resize(nnz);
        compressed_format.values.resize(nnz);

        // count non‑zeros per major index
        std::vector<size_t> counts;
        if constexpr (S == StorageOrder::ColumnMajor)
        {
            // resize the counts vector to the number of columns
            counts.resize(cols, 0);

            // count non-zeros per column using atomic operations
            std::for_each(
                std::execution::par,
                uncompressed_format.begin(), uncompressed_format.end(),
                [&](auto const &entry)
                {
                    size_t idx = entry.first.col;
                    __atomic_fetch_add(&counts[idx], 1, __ATOMIC_RELAXED);
                });
        }
        else
        {
            // resize the counts vector to the number of rows
            counts.resize(rows, 0);

            // count non-zeros per row using atomic operations
            std::for_each(
                std::execution::par,
                uncompressed_format.begin(), uncompressed_format.end(),
                [&](auto const &entry)
                {
                    size_t idx = entry.first.row;
                    __atomic_fetch_add(&counts[idx], 1, __ATOMIC_RELAXED);
                });
        }

        // exclusive prefix‐sum to build the "inner" index array
        compressed_format.inner[0] = 0;
        // Create a temporary vector to hold non-atomic values for the scan
        std::vector<size_t> temp_counts(counts.size());
        std::transform(std::execution::par, counts.begin(), counts.end(), temp_counts.begin(),
                       [](const std::atomic<size_t> &atomic_val)
                       { return atomic_val.load(); });

        // Perform the inclusive scan on the temporary vector
        std::inclusive_scan(std::execution::par, temp_counts.begin(), temp_counts.end(), compressed_format.inner.begin() + 1);

        // allocate storage for the "outer" indices and values
        compressed_format.outer.resize(nnz);
        compressed_format.values.resize(nnz);

        // fill the "outer" indices and values
        std::vector<size_t> indices(uncompressed_format.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::for_each(
            std::execution::par,
            indices.begin(), indices.end(),
            [&](size_t index)
            {
                auto it = std::next(uncompressed_format.begin(), index);
                if constexpr (S == StorageOrder::ColumnMajor)
                {
                    compressed_format.outer[index] = it->first.row;
                }
                else
                {
                    compressed_format.outer[index] = it->first.col;
                }
                compressed_format.values[index] = it->second;
            });

        // clear the compressed matrix
        uncompressed_format.clear();

        // update the compressed flag
        compressed = true;
    }

    template <AddMulType T, StorageOrder S>
    void Matrix<T, S>::uncompress()
    {
        if (not compressed)
            return;

        // clear the uncompressed matrix
        uncompressed_format.clear();

        // fill the uncompressed matrix
        if constexpr (S == StorageOrder::ColumnMajor)
        {
            // iterate over columns of m
            for (size_t col_idx = 0; col_idx < cols; col_idx++)
            {
                // iterate over rows of m that are non-zero in the column "col" of m
                size_t start = compressed_format.inner[col_idx];
                size_t end = compressed_format.inner[col_idx + 1];
                for (size_t j = start; j < end; j++)
                {
                    // get the row index of the non-zero element
                    size_t row_idx = compressed_format.outer[j];

                    // add the non-zero element to the uncompressed matrix
                    uncompressed_format[{row_idx, col_idx}] = compressed_format.values[j];
                }
            }
        }
        else
        {
            // iterate over rows of m
            for (size_t row_idx = 0; row_idx < rows; row_idx++)
            {
                // iterate over columns of m that are non-zero in the row "row" of m
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
    };

    template <AddMulType T, StorageOrder S>
    T Matrix<T, S>::operator()(size_t row, size_t col) const
    {
        // check if the index is in range
        if (row >= rows or col >= cols)
        {
            throw std::out_of_range("Index out of range");
        }
        // check if the matrix is compressed
        if (not compressed)
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

    template <AddMulType T, StorageOrder S>
    Proxy<T, S> Matrix<T, S>::operator()(size_t row, size_t col)
    {
        if (row >= rows or col >= cols)
        {
            throw std::out_of_range("Index out of range");
        }
        if (compressed)
        {
            std::cout << "Matrix is compressed, uncompressing...\n";
            uncompress();
        }
        return Proxy<T, S>{uncompressed_format, row, col};
    }

    template <AddMulType T, StorageOrder S>
    void Matrix<T, S>::resize_and_clear(size_t rows, size_t cols)
    {
        this->rows = rows;
        this->cols = cols;

        compressed = false; // default value
        uncompressed_format.clear();
        compressed_format.inner.clear();
        compressed_format.outer.clear();
        compressed_format.values.clear();
    }

    template <AddMulType T, StorageOrder S>
    template <NormType N>
    double Matrix<T, S>::norm() const
    {
        if (typeid(*this) == typeid(SquareMatrix<T, S>))
        {
            auto *this_square = static_cast<const SquareMatrix<T, S> *>(this);
            if (this_square->is_modified())
            {
                return this_square->template norm<N>();
            }
        }
        if (not compressed)
        {
            if constexpr (N == NormType::One)
            {
                std::vector<double> col_sums(cols, 0);
                for (const auto &it : uncompressed_format)
                {
                    col_sums[it.first.col] += std::abs(it.second);
                }
                return *std::max_element(std::execution::par_unseq, col_sums.begin(), col_sums.end());
            }
            else if constexpr (N == NormType::Infinity)
            {
                std::vector<double> row_sums(rows, 0);
                for (const auto &it : uncompressed_format)
                {
                    row_sums[it.first.row] += std::abs(it.second);
                }
                return *std::max_element(std::execution::par_unseq, row_sums.begin(), row_sums.end());
            }
            else
            {
                double norm = 0;
                for (const auto &it : uncompressed_format)
                {
                    norm += std::abs(it.second) * std::abs(it.second);
                }
                return std::sqrt(std::abs(norm));
            }
        }
        if constexpr (N == NormType::One)
        {
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                std::vector<double> col_sums(cols, 0);
                for (size_t col = 0; col < cols; col++)
                {
                    size_t start = compressed_format.inner[col];
                    size_t end = compressed_format.inner[col + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        col_sums[col] += std::abs(compressed_format.values[j]);
                    }
                }
                return *std::max_element(std::execution::par_unseq, col_sums.begin(), col_sums.end());
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
                return *std::max_element(std::execution::par_unseq, col_sums.begin(), col_sums.end());
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
                std::vector<double> row_sums(rows, 0);
                for (size_t row = 0; row < rows; row++)
                {
                    size_t start = compressed_format.inner[row];
                    size_t end = compressed_format.inner[row + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        row_sums[row] += std::abs(compressed_format.values[j]);
                    }
                }
                return *std::max_element(std::execution::par_unseq, row_sums.begin(), row_sums.end());
            }
        }
        else
        {
            double norm = 0;
            for (const auto &it : compressed_format.values)
            {
                norm += std::abs(it) * std::abs(it);
            }
            return std::sqrt(norm);
        }
    };

    template <AddMulType T, StorageOrder S>
    void Matrix<T, S>::reader(const std::string &filename)
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

        // Resize the matrix
        resize_and_clear(row_read, col_read);

        // Read matrix values
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            T value;
            size_t row, col;
            iss >> row >> col >> value;
            assert(row <= rows and col <= cols);
            // I traslate the row and column indices to 0-based format and set the element
            //  in the matrix
            set(row - 1, col - 1, value);
        }

        file.close();
    };

    template <AddMulType T, StorageOrder S>
    std::vector<T> operator*(const Matrix<T, S> &m, const std::vector<T> &v)
    {
        if (m.cols != v.size())
        {
            throw std::invalid_argument("Matrix and vector dimensions do not match for multiplication");
        }
        std::vector<T> result(m.rows, 0);
        if (not m.is_compressed())
        {
            for (const auto &it : m.uncompressed_format)
            {
                result[it.first.row] += it.second * v[it.first.col];
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
                        result[row] += m.compressed_format.values[j] * v[col];
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
                        result[row] += m.compressed_format.values[j] * v[col];
                    }
                }
            }
        }
        return result;
    }

    template <AddMulType T, StorageOrder S>
    Matrix<T, S> operator*(const Matrix<T, S> &m1, const Matrix<T, S> &m2)
    {
        if (m1.cols != m2.rows)
        {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }
        if (m1.is_compressed() != m2.is_compressed())
        {
            throw std::invalid_argument("Matrix compression formats do not match");
        }

        Matrix<T, S> result(m1.rows, m2.cols);

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
        else
        {
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
        }
        return result;
    }

    template <AddMulType T, StorageOrder S>
    size_t Matrix<T, S>::get_nnz() const
    {
        if (compressed)
        {
            return compressed_format.values.size();
        }
        else
        {
            return uncompressed_format.size();
        }
    };
}

#endif // MATRIX_TPP
