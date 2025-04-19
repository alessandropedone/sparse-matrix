#include "matrix.hpp"
#include "matrix_views.hpp"

using namespace algebra;

int main() {
    Matrix<int, StorageOrder::RowMajor> m1(3, 3);
    m1.set(0, 0, 1);
    m1(2,2) = 4;
    std::vector<int> v = {1, 2, 3};

    std::cout << "Matrix M" << std::endl;
    m1.compress();
    const auto & ref = m1;
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            std::cout << ref(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl; 

    std::cout << "Vector v" << std::endl;
    for (size_t i = 0; i < 3; i++)
    {
        std::cout << v[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "M*v" << std::endl;
    auto result = m1 * v;
    for (size_t i = 0; i < 3; i++)
    {
        std::cout << result[i] << std::endl;
    }
    std::cout << std::endl;
    return 0;
}