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
    SquareMatrix<double, StorageOrder::ColumnMajor> sm(0);
    test5x5(m);
    test5x5(sm);

    // Run the tests on all matrices provided in the json file for both storage orders
    /*     json data = read_json(static_cast<std::string>("data/data.json"));
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

        } */

    return 0;
}