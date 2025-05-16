/**
 * @file test.hpp
 * @brief Utility functions for testing and benchmarking matrix classes and operations.
 *
 * This header provides a collection of templated utility functions for:
 * - Printing vectors and matrices.
 * - Comparing matrices for equality.
 * - Testing compression and decompression of matrices (including SquareMatrix, Matrix, TransposeView, DiagonalView).
 * - Computing matrix norms (One, Infinity, Frobenius).
 * - Generating random vectors for testing.
 * - Running comprehensive tests on 5x5 matrices, including matrix-vector and matrix-matrix products.
 * - Measuring and reporting execution times for matrix operations in both compressed and uncompressed formats.
 * - Saving timing results to JSON files for further analysis.
 *
 * The utilities are designed to work with custom matrix types and storage orders, supporting both real and complex types.
 * 
 * @author 
 * @date 
 */
#ifndef TEST_HPP
#define TEST_HPP

#include <string>
#include <chrono>
#include <random>

#include "json_utility.hpp"
#include "storage.hpp"
#include "abstract_matrix.hpp"
#include "matrix_views.hpp"
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
            std::cout << std::setw(20) << v[i] << std::endl;
        }
        std::cout << std::endl;
    }

    /// @brief print a matrix
    /// @tparam T type of the matrix elements
    /// @tparam S type of the storage order (RowMajor or ColumnMajor)
    /// @param m matrix to print
    template <AddMulType T, StorageOrder S>
    void print(const AbstractMatrix<T, S> &m)
    {
        for (size_t i = 0; i < m.get_rows(); i++)
        {
            for (size_t j = 0; j < m.get_cols(); j++)
            {
                std::cout << std::setw(20) << m(i, j) << " ";
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
    bool are_equal(const AbstractMatrix<T, S> &m1, const AbstractMatrix<T, S> &m2)
    {
        if (m1.get_rows() != m2.get_rows() or m1.get_cols() != m2.get_cols())
        {
            return false;
        }
        for (size_t i = 0; i < m1.get_rows(); i++)
        {
            for (size_t j = 0; j < m1.get_cols(); j++)
            {
                if (std::abs(static_cast<T>(m1(i, j) - m2(i, j))) > std::numeric_limits<AbsReturnType_t<T>>::epsilon())
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
    bool test_compression_square_matrix(const SquareMatrix<T, S> &m)
    {
        // Read matrix
        SquareMatrix<T, S> compare_matrix(m);

        compare_matrix.compress();
        if (not are_equal(m, compare_matrix))
        {
            throw std::runtime_error("Error passing from uncompressed to compressed format");
            return false;
        }

        compare_matrix.compress_mod();
        if (not are_equal(m, compare_matrix))
        {
            throw std::runtime_error("Error passing from compressed to modified compressed format");
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

        compare_matrix.uncompress();
        if (not are_equal(m, compare_matrix))
        {
            throw std::runtime_error("Error passing from modified compressed to uncompressed format");
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
    bool test_compression_matrix(const AbstractMatrix<T, S> &m)
    {
        auto compare_matrix = m.clone();

        compare_matrix->compress();
        if (not are_equal(m, *compare_matrix))
        {
            throw std::runtime_error("Error passing from uncompressed to compressed format");
            return false;
        }

        compare_matrix->uncompress();
        if (not are_equal(m, *compare_matrix))
        {
            throw std::runtime_error("Error passing from compressed to uncompressed format");
            return false;
        }

        std::cout << "Compression test passed" << std::endl;
        std::cout << std::endl;
        return true;
    }

    /// @brief test compression of a matrix
    /// @tparam T type of the matrix elements
    /// @tparam S type of the storage order (RowMajor or ColumnMajor)
    /// @param m matrix
    /// @note this function is a friend of the Matrix class, so it can access the private members
    template <AddMulType T, StorageOrder S>
    void test_compression(const AbstractMatrix<T, S> &m)
    {
        if (typeid(m) == typeid(SquareMatrix<T, S>))
        {
            test_compression_square_matrix(static_cast<const SquareMatrix<T, S> &>(m));
        }
        else
        {
            test_compression_matrix(m);
        }
    }

    /// @brief compute the norm of a matrix
    /// @tparam T type of the matrix elements
    /// @tparam S type of the storage order (RowMajor or ColumnMajor)
    /// @param m matrix
    template <AddMulType T, StorageOrder S>
    void norm_test(const AbstractMatrix<T, S> &m)
    {
        // Read matrix
        std::cout << "Matrix norms" << std::endl;
        std::cout << "One norm:       " << std::setw(14) << m.template norm<NormType::One>() << std::endl;
        std::cout << "Infinity norm:  " << std::setw(14) << m.template norm<NormType::Infinity>() << std::endl;
        std::cout << "Frobenius norm: " << std::setw(14) << m.template norm<NormType::Frobenius>() << std::endl;
        std::cout << std::endl;
        return;
    }

    template <AddMulType T, StorageOrder S>
    void print_result(const AbstractMatrix<T, S> &m, const std::vector<T> &v)
    {
        std::cout << "M*v" << std::endl;
        print(v);
        std::cout << "M^2 " << std::endl;
        print(m);
    }

    template <AddMulType T>
    void generateRandomVector(std::vector<T> &vec)
    {
        std::random_device seed;
        unsigned int constant_seed = 3; // Set a constant seed for reproducibility
        std::default_random_engine gen(constant_seed);

        if constexpr (std::is_floating_point<T>::value)
        {
            std::uniform_real_distribution<T> distr(-1.0, 1.0);
            for (auto &val : vec)
            {
                val = distr(gen);
            }
        }
        else if constexpr (std::is_integral<T>::value)
        {
            std::uniform_int_distribution<T> distr(-1, 1);
            for (auto &val : vec)
            {
                val = distr(gen);
            }
        }
        else if constexpr (is_complex<T>::value)
        {
            using RealType = typename T::value_type; // Get the underlying type of the complex number
            std::uniform_real_distribution<RealType> distr(-1.0, 1.0);
            for (auto &val : vec)
            {
                val.real(distr(gen)); // Generate random value for the real part
                val.imag(distr(gen)); // Generate random value for the imaginary part
            }
        }
        else
        {
            static_assert(std::is_arithmetic<T>::value, "Unsupported type T for testing");
        }
    }

    /// @brief 5x5 matrix for testing
    /// @tparam T type of the matrix elements
    /// @tparam S type of the storage order (RowMajor or ColumnMajor)
    template <AddMulType T, StorageOrder S>
    void test5x5(AbstractMatrix<T, S> &m, const std::string &name)
    {
        m.reader(static_cast<std::string>("data/" + name));

        if (typeid(m) == typeid(SquareMatrix<T, S>))
        {
            std::cout << "----------------------------" << std::endl;
            std::cout << "Test with SquareMatrix class" << std::endl;
            std::cout << "----------------------------" << std::endl;
        }
        else if (typeid(m) == typeid(TransposeView<T, S>))
        {
            std::cout << "------------------------" << std::endl;
            std::cout << "Test with TransposeView" << std::endl;
            std::cout << "------------------------" << std::endl;
        }
        else if (typeid(m) == typeid(DiagonalView<T, S>))
        {
            std::cout << "-------------------------" << std::endl;
            std::cout << "Test with DiagonalView" << std::endl;
            std::cout << "-------------------------" << std::endl;
        }
        else
        {
            std::cout << "----------------------" << std::endl;
            std::cout << "Test with Matrix class" << std::endl;
            std::cout << "----------------------" << std::endl;
        }
        // Print the matrix
        std::cout << "Test 5x5 matrix" << std::endl;
        print(m);

        // Test compression
        test_compression(m);

        // Test norm functions
        norm_test(m);

        // Initialize a vector
        std::vector<T> v(m.get_cols(), T(0));
        generateRandomVector(v);

        // Print the vector
        std::cout << "Test vector" << std::endl;
        print(v);

        if (typeid(m) == typeid(SquareMatrix<T, S>))
        {
            auto sm = static_cast<const SquareMatrix<T, S> &>(m);
            sm.compress_mod();
            // Do the matrix - vector product
            auto result = sm * v;
            // Do the matrix - matrix product
            auto m2 = sm * sm;
            // Print the result
            print_result(m2, result);
        }
        else if (typeid(m) == typeid(TransposeView<T, S>))
        {
            auto tm = static_cast<const TransposeView<T, S> &>(m);
            if (typeid(tm.matrix) == typeid(SquareMatrix<T, S>))
            {
                auto sm = static_cast<const SquareMatrix<T, S> &>(tm.matrix);
                sm.compress_mod();
            }
            else
            {
                tm.compress();
            }
            // Do the matrix - vector product
            auto result = tm * v;
            // Do the matrix - matrix product
            auto m2 = tm * tm;
            // Print the result
            print_result(m2, result);
        }
        else if (typeid(m) == typeid(DiagonalView<T, S>))
        {
            auto dm = static_cast<const DiagonalView<T, S> &>(m);
            dm.compress();
            // Do the matrix - vector product
            auto result = dm * v;
            // Do the matrix - matrix product
            auto m2 = dm * dm;
            // Print the result
            print_result(m2, result);
        }
        else
        {
            m.compress();
            auto mm = static_cast<const Matrix<T, S> &>(m);
            // Do the matrix - vector product
            auto result = mm * v;
            // Do the matrix - matrix product
            auto m2 = mm * mm;
            // Print the result
            print_result(m2, result);
        }
        m.uncompress();
    }

    /// @brief test the execution time of matrix-matrix and matrix-vector products
    /// @tparam S type of the storage order (RowMajor or ColumnMajor)
    /// @param matrix_names vector of matrix names
    /// @note the matrix names must be in the data.json file inside the data folder
    template <AddMulType T, StorageOrder S>
    void test(const std::vector<std::string> &matrix_names)
    {

        for (const auto &matrix_name : matrix_names)
        {

            std::cout << std::endl;
            Matrix<T, S> testMatrix(0, 0);
            SquareMatrix<T, S> testSquareMatrix(0);
            Matrix<T, S> temp(0, 0);
            SquareMatrix<T, S> temp2(0);
            TransposeView<T, S> testTransposeView(temp);
            DiagonalView<T, S> testDiagonalView(temp2);

            std::cout << "------------------------------------" << std::endl;
            std::cout << "Test with Matrix class" << std::endl;
            std::cout << "------------------------------------" << std::endl;
            execute_test(testMatrix, matrix_name);

            std::cout << "------------------------------------" << std::endl;
            std::cout << "Test with SquareMatrix class" << std::endl;
            std::cout << "------------------------------------" << std::endl;
            execute_test(testSquareMatrix, matrix_name);

            std::cout << "------------------------------------" << std::endl;
            std::cout << "Test with TransposeView class" << std::endl;
            std::cout << "------------------------------------" << std::endl;
            execute_test(testTransposeView, matrix_name);

            std::cout << "------------------------------------" << std::endl;
            std::cout << "Test with DiagonalView class" << std::endl;
            std::cout << "------------------------------------" << std::endl;
            execute_test(testDiagonalView, matrix_name);
        };
    }

    template <AddMulType T, StorageOrder S>
    void execute_test(AbstractMatrix<T, S> &testMatrix, const std::string &matrix_name)
    {
        // Import matrix
        testMatrix.reader(static_cast<std::string>("data/" + matrix_name));

        std::cout << "Test matrix " << matrix_name << std::endl;
        std::cout << std::endl;

        // Test compression
        test_compression(testMatrix);

        // Test norm
        norm_test(testMatrix);

        // Test excution time of products
        std::cout << "Test for execution time of products" << std::endl;

        // Generate vector to perform the matrix - vector product
        std::vector<T> vec(testMatrix.get_cols(), T(0));
        generateRandomVector(vec);

        Matrix<T, S> res1(testMatrix.get_rows(), testMatrix.get_cols());
        Matrix<T, S> res2(testMatrix.get_rows(), testMatrix.get_cols());
        std::vector<T> res3(testMatrix.get_rows(), 0);
        std::vector<T> res4(testMatrix.get_rows(), 0);

        MyTimePoint start, stop;
        std::string filename = "data/execution_time.json";
        json time_info = read_json(filename);

        if (typeid(testMatrix) == typeid(SquareMatrix<T, S>))
        {
            auto testSquareMatrix = static_cast<SquareMatrix<T, S> &>(testMatrix);

            testSquareMatrix.compress_mod();

            start = MyClock::now();
            res1 = testSquareMatrix * testSquareMatrix;
            stop = MyClock::now();
            auto time_span_mu = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (compressed_format_matrix_matrix_product_mus)"] = time_span_mu.count();

            start = MyClock::now();
            res3 = testSquareMatrix * vec;
            stop = MyClock::now();
            auto time_span_n = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (compressed_format_matrix_vector_product_ns)"] = time_span_n.count();

            testSquareMatrix.uncompress();

            start = MyClock::now();
            res2 = testSquareMatrix * testSquareMatrix;
            stop = MyClock::now();
            auto time_span_mu_uncompressed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (uncompressed_format_matrix_matrix_product_mus)"] = time_span_mu_uncompressed.count();

            start = MyClock::now();
            res4 = testSquareMatrix * vec;
            stop = MyClock::now();
            auto time_span_n_uncompressed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (uncompressed_format_matrix_vector_product_ns)"] = time_span_n_uncompressed.count();
        }
        else if (typeid(testMatrix) == typeid(TransposeView<T, S>))
        {
            auto testTransposeView = static_cast<TransposeView<T, S> &>(testMatrix);

            testTransposeView.compress();

            start = MyClock::now();
            res1 = testTransposeView * testTransposeView;
            stop = MyClock::now();
            auto time_span_mu = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (compressed_format_matrix_matrix_product_mus)"] = time_span_mu.count();

            start = MyClock::now();
            res3 = testTransposeView * vec;
            stop = MyClock::now();
            auto time_span_n = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (compressed_format_matrix_vector_product_ns)"] = time_span_n.count();

            testTransposeView.uncompress();

            start = MyClock::now();
            res2 = testTransposeView * testTransposeView;
            stop = MyClock::now();
            auto time_span_mu_uncompressed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (uncompressed_format_matrix_matrix_product_mus)"] = time_span_mu_uncompressed.count();

            start = MyClock::now();
            res4 = testTransposeView * vec;
            stop = MyClock::now();
            auto time_span_n_uncompressed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (uncompressed_format_matrix_vector_product_ns)"] = time_span_n_uncompressed.count();
        }
        else if (typeid(testMatrix) == typeid(DiagonalView<T, S>))
        {
            auto testDiagonalView = static_cast<DiagonalView<T, S> &>(testMatrix);

            testDiagonalView.compress();

            start = MyClock::now();
            res1 = testDiagonalView * testDiagonalView;
            stop = MyClock::now();
            auto time_span_mu = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (compressed_format_matrix_matrix_product_mus)"] = time_span_mu.count();

            start = MyClock::now();
            res3 = testDiagonalView * vec;
            stop = MyClock::now();
            auto time_span_n = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (compressed_format_matrix_vector_product_ns)"] = time_span_n.count();

            testDiagonalView.uncompress();

            start = MyClock::now();
            res2 = testDiagonalView * testDiagonalView;
            stop = MyClock::now();
            auto time_span_mu_uncompressed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (uncompressed_format_matrix_matrix_product_mus)"] = time_span_mu_uncompressed.count();

            start = MyClock::now();
            res4 = testDiagonalView * vec;
            stop = MyClock::now();
            auto time_span_n_uncompressed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (uncompressed_format_matrix_vector_product_ns)"] = time_span_n_uncompressed.count();
        }
        else
        {
            auto dynamicMatrix = dynamic_cast<Matrix<T, S> *>(&testMatrix);
            if (!dynamicMatrix)
            {
                throw std::runtime_error("No matrix type found to perform the test");
            }

            dynamicMatrix->compress();

            start = MyClock::now();
            res1 = (*dynamicMatrix) * (*dynamicMatrix);
            stop = MyClock::now();
            auto time_span_mu = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (compressed_format_matrix_matrix_product_mus)"] = time_span_mu.count();

            start = MyClock::now();
            res3 = (*dynamicMatrix) * vec;
            stop = MyClock::now();
            auto time_span_n = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (compressed_format_matrix_vector_product_ns)"] = time_span_n.count();

            dynamicMatrix->uncompress();

            start = MyClock::now();
            res2 = (*dynamicMatrix) * (*dynamicMatrix);
            stop = MyClock::now();
            auto time_span_mu_uncompressed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (uncompressed_format_matrix_matrix_product_mus)"] = time_span_mu_uncompressed.count();

            start = MyClock::now();
            res4 = (*dynamicMatrix) * vec;
            stop = MyClock::now();
            auto time_span_n_uncompressed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            time_info[matrix_name + " " + typeid(testMatrix).name() + " (uncompressed_format_matrix_vector_product_ns)"] = time_span_n_uncompressed.count();
        }

        // save json
        save_json(filename, time_info);

        // Print execution times and speedups
        std::cout << std::endl;

        int compressed_matrix_vector_time = time_info[matrix_name + " " + typeid(testMatrix).name() + " (compressed_format_matrix_vector_product_ns)"];
        int uncompressed_matrix_vector_time = time_info[matrix_name + " " + typeid(testMatrix).name() + " (uncompressed_format_matrix_vector_product_ns)"];

        std::cout << "Compressed format matrix-vector product time: " << compressed_matrix_vector_time << " ns" << std::endl;
        std::cout << "Uncompressed format matrix-vector product time: " << uncompressed_matrix_vector_time << " ns" << std::endl;
        std::cout << "Matrix-vector product speedup: "
                  << static_cast<double>(uncompressed_matrix_vector_time) / compressed_matrix_vector_time << std::endl;

        std::cout << std::endl;

        int compressed_matrix_matrix_time = time_info[matrix_name + " " + typeid(testMatrix).name() + " (compressed_format_matrix_matrix_product_mus)"];
        int uncompressed_matrix_matrix_time = time_info[matrix_name + " " + typeid(testMatrix).name() + " (uncompressed_format_matrix_matrix_product_mus)"];

        std::cout << "Compressed format matrix-matrix product time: " << compressed_matrix_matrix_time << " µs" << std::endl;
        std::cout << "Uncompressed format matrix-matrix product time: " << uncompressed_matrix_matrix_time << " µs" << std::endl;
        std::cout << "Matrix-matrix product speedup: "
                  << static_cast<double>(uncompressed_matrix_matrix_time) / compressed_matrix_matrix_time << std::endl;

        std::cout << std::endl;
    }

}

#endif // TEST_HPP