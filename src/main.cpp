#include "matrix.hpp"
#include "square_matrix.hpp"
#include "matrix_views.hpp"
#include "json_utility.hpp"
#include "square_matrix.hpp"
#include "test.hpp"

#include <random>
#include <iostream>
#include <chrono>

using namespace algebra;
using namespace json_utility;

int main()
{
    // Import matrix
    Matrix<double, StorageOrder::ColumnMajor> m(0, 0);
    m.reader(static_cast<std::string>("data/read_test_5x5.mtx"));

    // Print the matrix
    std::cout << "Test 5x5 matrix" << std::endl;
    print(m);

    // Test compression
    test_compression_matrix(m);

    // Test norm functions
    norm_test(m);

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
    std::cout << "Test vector" << std::endl;
    print(v);

    // Do the matrix - vector product
    auto result = m * v;
    std::cout << "M*v" << std::endl;
    print(result);

    // Do the matrix - matrix product
    std::cout << "M^2 " << std::endl;

    // Print the result
    auto m2 = m * m;
    print(m2);

    // Run the tests on all matrices provided in the json file for both storage orders
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
            test<StorageOrder::RowMajor>(matrix_names);
        }
        else
        {
            std::cout << "------------------------------------" << std::endl;
            std::cout << "Test with storage order: ColumnMajor" << std::endl;
            std::cout << "------------------------------------" << std::endl;
            test<StorageOrder::ColumnMajor>(matrix_names);
        }
        
    }
    
    return 0;
}