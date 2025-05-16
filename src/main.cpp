/**
 * @file main.cpp
 * @brief Entry point for testing various matrix types and storage orders.
 *
 * This program performs a series of tests on different matrix classes, including real and complex matrices,
 * square matrices, and matrix views (transpose and diagonal). It reads matrix data from files and a JSON
 * configuration, and runs tests for both row-major and column-major storage orders.
 *
 * The following matrix types are tested:
 * - Matrix<double, StorageOrder::ColumnMajor>
 * - SquareMatrix<double, StorageOrder::ColumnMajor>
 * - TransposeView<double, StorageOrder::ColumnMajor>
 * - DiagonalView<double, StorageOrder::ColumnMajor>
 * - Matrix<std::complex<double>, StorageOrder::ColumnMajor>
 * - SquareMatrix<std::complex<double>, StorageOrder::ColumnMajor>
 * - TransposeView<std::complex<double>, StorageOrder::ColumnMajor>
 * - DiagonalView<std::complex<double>, StorageOrder::ColumnMajor>
 *
 * The program also reads a list of matrix names from a JSON file and runs tests on all matrices for both
 * row-major and column-major storage orders.
 *
 * @return int Returns 0 upon successful completion.
 */
#include "abstract_matrix.hpp"
#include "matrix.hpp"
#include "square_matrix.hpp"
#include "matrix_views.hpp"
#include "json_utility.hpp"
#include "square_matrix.hpp"
#include "test.hpp"

#include <random>
#include <iostream>
#include <chrono>
#include <complex>

using namespace algebra;
using namespace json_utility;

int main()
{
    // Test with a 5x5 matrix
    std::cout << "------------------------------------" << std::endl;
    std::cout << "Test with a 5x5 real matrix" << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::string filename = "real_test_5x5.mtx";
    Matrix<double, StorageOrder::ColumnMajor> m(0, 0);
    test5x5(m, filename);
    SquareMatrix<double, StorageOrder::ColumnMajor> sm(0);
    test5x5(sm, filename);
    TransposeView<double, StorageOrder::ColumnMajor> tv(0, 0);
    test5x5(tv, filename);
    DiagonalView<double, StorageOrder::ColumnMajor> dv(0, 0);
    test5x5(dv, filename);

    std::cout << "------------------------------------" << std::endl;
    std::cout << "Test with a 5x5 complex matrix" << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::string filename2 = "complex_test_5x5.mtx";
    Matrix<std::complex<double>, StorageOrder::ColumnMajor> cm(0, 0);
    test5x5(cm, filename2);
    SquareMatrix<std::complex<double>, StorageOrder::ColumnMajor> csm(0);
    test5x5(csm, filename2);
    TransposeView<std::complex<double>, StorageOrder::ColumnMajor> ctv(0, 0);
    test5x5(ctv, filename2);
    DiagonalView<std::complex<double>, StorageOrder::ColumnMajor> cdv(0, 0);
    test5x5(cdv, filename2);

    // Run the tests on all matrices provided in the json file for both storage orders
    json data = read_json(static_cast<std::string>("data/data.json"));
    const std::vector<std::string> matrix_names = data["matrix_name"];
    const std::array<StorageOrder, 2> storage_orders = {StorageOrder::RowMajor, StorageOrder::ColumnMajor};
    for (const auto &storage_order : storage_orders)
    {
        if (storage_order == StorageOrder::RowMajor)
        {
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << "---------------------------------" << std::endl;
            std::cout << "Test with storage order: RowMajor" << std::endl;
            std::cout << "---------------------------------" << std::endl;
            test<double, StorageOrder::RowMajor>(matrix_names);
        }
        else
        {
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << "------------------------------------" << std::endl;
            std::cout << "Test with storage order: ColumnMajor" << std::endl;
            std::cout << "------------------------------------" << std::endl;
            test<double, StorageOrder::ColumnMajor>(matrix_names);
        }
    }

    return 0;
}