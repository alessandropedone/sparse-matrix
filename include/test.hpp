#ifndef TEST_HPP
#define TEST_HPP

#include <string>
#include <chrono>
#include <random>

#include "json_utility.hpp"
#include "storage.hpp"
#include "matrix.hpp"
#include "square_matrix.hpp"

using namespace json_utility;

using MyClock = std::chrono::high_resolution_clock;
using MyTimePoint = std::chrono::time_point<MyClock>;

namespace algebra
{
    /// @brief print a vector
    /// @tparam T type of the vector elements
    /// @param v vector to print
    template <AddMulType T>
    void print(const std::vector<T> &v)
    {
        for (size_t i = 0; i < v.size(); i++)
        {
            std::cout << std::setw(15) << v[i] << std::endl;
        }
        std::cout << std::endl;
    }

    /// @brief print a matrix
    /// @tparam T type of the matrix elements
    /// @tparam S type of the storage order (RowMajor or ColumnMajor)
    /// @param m matrix to print
    template <AddMulType T, StorageOrder S>
    void print(const Matrix<T, S> &m)
    {
        for (size_t i = 0; i < m.get_rows(); i++)
        {
            for (size_t j = 0; j < m.get_cols(); j++)
            {
                std::cout << std::setw(15) << m(i, j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    /// @brief test if two matrices are equal
    /// @tparam T type of the matrix elements
    /// @tparam S type of the storage order (RowMajor or ColumnMajor)
    /// @param m1 first matrix
    /// @param m2 second matrix
    /// @return true if the matrices are equal, false otherwise
    template <AddMulType T, StorageOrder S>
    bool are_equal(const Matrix<T, S> &m1, const Matrix<T, S> &m2)
    {
        if (m1.get_rows() != m2.get_rows() or m1.get_cols() != m2.get_cols())
        {
            return false;
        }
        for (size_t i = 0; i < m1.get_rows(); i++)
        {
            for (size_t j = 0; j < m1.get_cols(); j++)
            {
                if (std::abs(m1(i, j) - m2(i, j)) > std::numeric_limits<T>::epsilon())
                {
                    return false;
                }
            }
        }
        return true;
    }

    /// @brief test if the matrix is compressed and uncompressed correctly
    /// @tparam T type of the matrix elements
    /// @tparam S type of the storage order (RowMajor or ColumnMajor)
    /// @param m matrix
    /// @return true if the test passed, false otherwise
    template <AddMulType T, StorageOrder S>
    bool test_compression_square_matrix(const Matrix<T, S> &m)
    {
        // Read matrix
        SquareMatrix<T, S> compare_matrix(m);

        compare_matrix.compress_mod();
        if (not are_equal(m, compare_matrix))
        {
            throw std::runtime_error("Error passing from uncompressed to modified compressed format");
            return false;
        }

        compare_matrix.compress();
        if (not are_equal(m, compare_matrix))
        {
            throw std::runtime_error("Error passing from modified compressed to compressed format");
            return false;
        }

        compare_matrix.uncompress();
        if (not are_equal(m, compare_matrix))
        {
            throw std::runtime_error("Error passing from compressed to uncompressed format");
            return false;
        }

        compare_matrix.compress_mod();
        if (not are_equal(m, compare_matrix))
        {
            throw std::runtime_error("Error passing from uncompressed to modified compressed format");
            return false;
        }

        compare_matrix.compress();
        if (not are_equal(m, compare_matrix))
        {
            throw std::runtime_error("Error passing from modified compressed to compressed format");
            return false;
        }

        std::cout << "Compression test passed" << std::endl;
        std::cout << std::endl;
        return true;
    }

    /// @brief test if the matrix is compressed and uncompressed correctly
    /// @tparam T type of the matrix elements
    /// @tparam S type of the storage order (RowMajor or ColumnMajor)
    /// @param m matrix
    /// @return true if the test passed, false otherwise
    template <AddMulType T, StorageOrder S>
    bool test_compression_matrix(const Matrix<T, S> &m)
    {
        // Read matrix
        Matrix<T, S> compare_matrix(m);

        compare_matrix.compress();
        if (not are_equal(m, compare_matrix))
        {
            throw std::runtime_error("Error passing from uncompressed to compressed format");
            return false;
        }

        compare_matrix.uncompress();
        if (not are_equal(m, compare_matrix))
        {
            throw std::runtime_error("Error passing from compressed to uncompressed format");
            return false;
        }

        std::cout << "Compression test passed" << std::endl;
        std::cout << std::endl;
        return true;
    }

    /// @brief compute the norm of a matrix
    /// @tparam T type of the matrix elements
    /// @tparam S type of the storage order (RowMajor or ColumnMajor)
    /// @param m matrix
    template <AddMulType T, StorageOrder S>
    void norm_test(const Matrix<T, S> &m)
    {
        // Read matrix
        std::cout << "Norm test" << std::endl;
        std::cout << "One norm:       " << std::setw(14) << m.template norm<NormType::One>() << std::endl;
        std::cout << "Infinity norm:  " << std::setw(14) << m.template norm<NormType::Infinity>() << std::endl;
        std::cout << "Frobenius norm: " << std::setw(14) << m.template norm<NormType::Frobenius>() << std::endl;
        std::cout << std::endl;
        return;
    }

    /// @brief test the execution time of matrix-matrix and matrix-vector products
    /// @tparam S type of the storage order (RowMajor or ColumnMajor)
    /// @param matrix_names vector of matrix names
    /// @note the matrix names must be in the data.json file inside the data folder
    template <StorageOrder S>
    void test(const std::vector<std::string> &matrix_names)
    {

        for (const auto &matrix_name : matrix_names)
        {

            std::cout << std::endl;
            Matrix<double, S> testMatrix(0, 0);
            SquareMatrix<double, S> testSquareMatrix(0);

            std::cout << "------------------------------------" << std::endl;
            std::cout << "Test with Matrix class" << std::endl;
            std::cout << "------------------------------------" << std::endl;
            execute_test(testMatrix, matrix_name);

            std::cout << "------------------------------------" << std::endl;
            std::cout << "Test with SquareMatrix class" << std::endl;
            std::cout << "------------------------------------" << std::endl;
            execute_test(testSquareMatrix, matrix_name);
        };
    }

    template <typename T, StorageOrder S>
    void execute_test(Matrix<T, S> &testMatrix, const std::string &matrix_name)
    {
        // Import matrix
        testMatrix.reader(static_cast<std::string>("data/" + matrix_name));

        std::cout << "Test matrix " << matrix_name << std::endl;
        std::cout << std::endl;

        // Test compression
        test_compression_matrix(testMatrix);

        // Test norm
        norm_test(testMatrix);

        // Test excution time of products
        std::cout << "Test for execution time of products" << std::endl;

        // Initialize random distribution to create a random filled vector
        std::random_device seed;
        std::default_random_engine gen(seed());
        std::uniform_real_distribution<double> distr(-1., 1.);

        // Generate vector to perform the matrix - vector product
        std::vector<double> vec(testMatrix.get_cols(), 0);
        for (auto &val : vec)
        {
            val = distr(gen);
        }

        Matrix<double, S> res1(testMatrix.get_rows(), testMatrix.get_cols());
        Matrix<double, S> res2(testMatrix.get_rows(), testMatrix.get_cols());
        std::vector<double> res3(testMatrix.get_rows(), 0);
        std::vector<double> res4(testMatrix.get_rows(), 0);

        MyTimePoint start, stop;
        std::string filename = "data/execution_time.json";
        json time_info = read_json(filename);

        // matrix - vector product in compressed format
        testMatrix.compress();

        start = MyClock::now();
        res1 = testMatrix * testMatrix;
        stop = MyClock::now();
        auto time_span_mu = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        time_info[matrix_name + " (compressed_format_matrix_matrix_product_mus)"] = time_span_mu.count();

        start = MyClock::now();
        res3 = testMatrix * vec;
        stop = MyClock::now();
        auto time_span_n = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        time_info[matrix_name + " (compressed_format_matrix_vector_product_ns)"] = time_span_n.count();

        // matrix - vector product in uncompressed format
        testMatrix.uncompress();

        start = MyClock::now();
        res2 = testMatrix * testMatrix;
        stop = MyClock::now();
        time_span_mu = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        time_info[matrix_name + " (uncompressed_format_matrix_matrix_product_mus)"] = time_span_mu.count();

        start = MyClock::now();
        res4 = testMatrix * vec;
        stop = MyClock::now();
        time_span_n = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        time_info[matrix_name + " (uncompressed_format_matrix_vector_product_ns)"] = time_span_n.count();

        // save json
        save_json(filename, time_info);

        // Print execution times and speedups
        std::cout << std::endl;

        int compressed_matrix_vector_time = time_info[matrix_name + " (compressed_format_matrix_vector_product_ns)"];
        int uncompressed_matrix_vector_time = time_info[matrix_name + " (uncompressed_format_matrix_vector_product_ns)"];

        std::cout << "Compressed format matrix-vector product time: " << compressed_matrix_vector_time << " ns" << std::endl;
        std::cout << "Uncompressed format matrix-vector product time: " << uncompressed_matrix_vector_time << " ns" << std::endl;
        std::cout << "Matrix-vector product speedup: "
                  << static_cast<double>(uncompressed_matrix_vector_time) / compressed_matrix_vector_time << std::endl;

        std::cout << std::endl;

        int compressed_matrix_matrix_time = time_info[matrix_name + " (compressed_format_matrix_matrix_product_mus)"];
        int uncompressed_matrix_matrix_time = time_info[matrix_name + " (uncompressed_format_matrix_matrix_product_mus)"];

        std::cout << "Compressed format matrix-matrix product time: " << compressed_matrix_matrix_time << " µs" << std::endl;
        std::cout << "Uncompressed format matrix-matrix product time: " << uncompressed_matrix_matrix_time << " µs" << std::endl;
        std::cout << "Matrix-matrix product speedup: "
                  << static_cast<double>(uncompressed_matrix_matrix_time) / compressed_matrix_matrix_time << std::endl;

        std::cout << std::endl;
    }

}

#endif // TEST_HPP