#include "matrix.hpp"
#include "matrix_views.hpp"

using namespace algebra;

int main()
{
    Matrix<int, StorageOrder::RowMajor> m1(3, 3);
    m1.set(0, 0, 1);
    m1.set(0, 1, 2);
    m1.set(0, 2, 3);
    m1.set(1, 0, 4);
    m1.set(1, 1, 5);
    m1.set(1, 2, 6);
    m1.set(2, 0, 3);
    m1.set(2, 1, 3);
    m1.set(2, 2, 3);
    std::vector<int> v = {1, 2, 3};

    std::cout << "Matrix M" << std::endl;
    m1.compress();
    //m1.uncompress();
    const auto &ref = m1;
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            std::cout << ref(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "One norm: " << m1.norm<NormType::One>() << std::endl;
    std::cout << "Infinity norm: " << m1.norm<NormType::Infinity>() << std::endl;
    std::cout << "Frobenius norm: " << m1.norm<NormType::Frobenius>() << std::endl;
    std::cout << std::endl;

    std::cout << "Vector v" << std::endl;
    for (size_t i = 0; i < 3; i++)
    {
        std::cout << v[i] << std::endl;
    }
    std::cout << std::endl;

    auto result = m1 * v;
    std::cout << "M*v" << std::endl;
    for (size_t i = 0; i < 3; i++)
    {
        std::cout << result[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "M^2 " << std::endl;
    auto m2 = m1 * m1;
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            std::cout << m2(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    return 0;
}