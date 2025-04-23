#include "matrix.hpp"
#include "square_matrix.hpp"
#include "matrix_views.hpp"
#include "json_utility.hpp"

#include <random>
#include <iostream>
#include <chrono>

using namespace algebra;
using namespace json_utility;
using MyClock = std::chrono::high_resolution_clock;
using MyTimePoint = std::chrono::time_point<MyClock>;

template <StorageOrder storage_order>
void test_storage_order(const std::vector<std::string> &matrix_names);

template <typename T, StorageOrder storage_order>
void test_square_matrix(const std::vector<std::string> &matrix_name);

template <typename T, StorageOrder S>
void print_matrix(const Matrix<T, S> &m);

int main()
{

    std::cout << "Test with 5x5 matrix" << std::endl;
    /*
    // Import matrix
    Matrix<double, StorageOrder::ColumnMajor> m(0, 0);
    m.reader(static_cast<std::string>("data/read_test_5x5.mtx"));

    m.compress_parallel();
    // Print the matrix
    std::cout << "Matrix M" << std::endl;
    const auto &ref = m;
    for (size_t i = 0; i < m.get_rows(); i++)
    {
        for (size_t j = 0; j < m.get_cols(); j++)
        {
            std::cout << std::setw(14) << ref(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Compute the norm of the matrix
    std::cout << "One norm:       " << std::setw(14) << m.norm<NormType::One>() << std::endl;
    std::cout << "Infinity norm:  " << std::setw(14) << m.norm<NormType::Infinity>() << std::endl;
    std::cout << "Frobenius norm: " << std::setw(14) << m.norm<NormType::Frobenius>() << std::endl;
    std::cout << std::endl;

    // Initialize a vector
    std::vector<double> v(m.get_cols(), 0);
    std::random_device seed;
    std::default_random_engine gen(seed());
    std::uniform_real_distribution<double> distr(-1., 1.);
    for (auto &val : v)
    {
        val = distr(gen);
    }

    // Print the vector
    std::cout << "Vector v" << std::endl;
    for (size_t i = 0; i < v.size(); i++)
    {
        std::cout << std::setw(14) << v[i] << std::endl;
    }
    std::cout << std::endl;

    // Do the matrix - vector product
    auto result = m * v;
    std::cout << "M*v" << std::endl;

    // Print the result
    for (size_t i = 0; i < result.size(); i++)
    {
        std::cout << std::setw(14) << result[i] << std::endl;
    }
    std::cout << std::endl;

    // Do the matrix - matrix product
    std::cout << "M^2 " << std::endl;

    // Print the result
    auto m2 = m * m;
    for (size_t i = 0; i < m2.get_rows(); i++)
    {
        for (size_t j = 0; j < m2.get_rows(); j++)
        {
            std::cout << std::setw(14) << m2(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    json data = read_json(static_cast<std::string>("data/data.json"));
    const std::vector<std::string> matrix_names = data["matrix_name"];

    const std::array<StorageOrder, 2> storage_orders = {StorageOrder::RowMajor, StorageOrder::ColumnMajor};
    for (const auto &storage_order : storage_orders)
    {
        if (storage_order == StorageOrder::RowMajor)
        {
            std::cout << "---------------------------------" << std::endl;
            std::cout << "Test with storage order: RowMajor" << std::endl;
            std::cout << "---------------------------------" << std::endl;
            test_storage_order<StorageOrder::RowMajor>(matrix_names);
        }
        else
        {
            std::cout << "------------------------------------" << std::endl;
            std::cout << "Test with storage order: ColumnMajor" << std::endl;
            std::cout << "------------------------------------" << std::endl;
            test_storage_order<StorageOrder::ColumnMajor>(matrix_names);
        }
        
    }*/

    test_square_matrix<double, StorageOrder::RowMajor>({"data/read_test_5x5.mtx"});

    test_square_matrix<double, StorageOrder::ColumnMajor>({"data/read_test_5x5.mtx"});

    return 0;
}

// test function with storage order
template <StorageOrder storage_order>
void test_storage_order(const std::vector<std::string> &matrix_names)
{

    for (const auto &matrix_name : matrix_names)
    {

        std::cout << std::endl;
        std::cout << "Test for execution time with matrix " << matrix_name << std::endl;
        Matrix<double, storage_order> testMatrix(0, 0);
        std::vector<double> v(testMatrix.get_cols(), 0);
        std::random_device seed;
        std::default_random_engine gen(seed());
        std::uniform_real_distribution<double> distr(-1., 1.);
        for (auto &val : v)
        {
            val = distr(gen);
        }

        // Import matrix
        testMatrix.reader(static_cast<std::string>("data/" + matrix_name));

        // Generate vector
        std::vector<double> vec(testMatrix.get_cols(), 0);
        Matrix<double, storage_order> res1(testMatrix.get_rows(), testMatrix.get_cols());
        Matrix<double, storage_order> res2(testMatrix.get_rows(), testMatrix.get_cols());
        std::vector<double> res3(testMatrix.get_rows(), 0);
        std::vector<double> res4(testMatrix.get_rows(), 0);
        for (auto &val : vec)
        {
            val = distr(gen);
        }

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

        // compress parallel vs compress
        testMatrix.uncompress();
        start = MyClock::now();
        testMatrix.compress();
        stop = MyClock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "Time for compress (μs): " << time_span.count() << std::endl;

        testMatrix.uncompress();
        start = MyClock::now();
        testMatrix.compress_parallel();
        stop = MyClock::now();
        time_span = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "Time for compress parallel (μs): " << time_span.count() << std::endl;

        std::cout << std::endl;
    };

}

template <typename T, StorageOrder storage_order>
void test_square_matrix(const std::vector<std::string> &matrix_name){
    // Read matrix
    SquareMatrix<T, storage_order> m(0);
    m.reader(static_cast<std::string>(matrix_name));
    print_matrix(m);

    // See all the data structures
    m.compress_mod();
    std::cout << "From uncompressed to modified compressed format\n" << std::endl;
    print_matrix(m);

    m.compress();
    std::cout << "From modified compressed to compressed format" << std::endl;
    print_matrix(m);

    m.uncompress();
    std::cout << "From compressed to uncompressed format" << std::endl;
    print_matrix(m);

    m.compress_mod();
    std::cout << "From uncompressed to modified compressed format" << std::endl;
    print_matrix(m);
}


template <typename T, StorageOrder S>
void print_matrix(const Matrix<T, S> &m){
    for (size_t i = 0; i < m.get_rows(); i++)
    {
        for (size_t j = 0; j < m.get_rows(); j++)
        {
            std::cout << std::setw(5) << m(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}