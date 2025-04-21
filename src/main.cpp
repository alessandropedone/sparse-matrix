#include "matrix.hpp"
#include "matrix_views.hpp"

#include <random>
#include <GetPot>
#include <iostream>

#include "json_utility.hpp"

using namespace algebra;
using namespace json_utility;

int main()
{

    std::cout << "Test for read method" << std::endl;

    Matrix<double, StorageOrder::RowMajor> m_r(0, 0);

    m_r.reader(static_cast<std::string>("read_test_5x5.mtx"));

    std::cout << "Matrix M_r" << std::endl;
    const auto &ref = m_r;
    for (size_t i = 0; i < m_r.get_rows(); i++)
    {
        for (size_t j = 0; j < m_r.get_cols(); j++)
        {
            std::cout << ref(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    Matrix<int, StorageOrder::RowMajor> m1(3, 3);
    m1.set(0, 0, 1);
    m1.set(0, 1, 2);
    m1.set(0, 2, 3);
    m1.set(1, 0, 0);
    m1.set(1, 1, 0);
    m1.set(1, 2, 0);
    m1.set(2, 0, 3);
    m1.set(2, 1, 3);
    m1.set(2, 2, 3);
    m1.set(2, 2, 0);
    m1(2, 2) = m1(0, 0) + 1;
    std::vector<int> v = {1, 2, 3};
    m1.compress_parallel();
    m1.uncompress();
    m1.compress();

    std::cout << "Matrix M" << std::endl;
    const auto &ref2 = m1;
    for (size_t i = 0; i < m1.get_cols(); i++)
    {
        for (size_t j = 0; j < m1.get_cols(); j++)
        {
            std::cout << ref2(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "One norm: " << m1.norm<NormType::One>() << std::endl;
    std::cout << "Infinity norm: " << m1.norm<NormType::Infinity>() << std::endl;
    std::cout << "Frobenius norm: " << m1.norm<NormType::Frobenius>() << std::endl;
    std::cout << std::endl;

    std::cout << "Vector v" << std::endl;
    for (size_t i = 0; i < v.size(); i++)
    {
        std::cout << v[i] << std::endl;
    }
    std::cout << std::endl;

    auto result = m1 * v;
    std::cout << "M*v" << std::endl;
    for (size_t i = 0; i < result.size(); i++)
    {
        std::cout << result[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "M^2 " << std::endl;
    auto m2 = m1 * m1;
    for (size_t i = 0; i < m2.get_rows(); i++)
    {
        for (size_t j = 0; j < m2.get_rows(); j++)
        {
            std::cout << m2(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Test for execution time with matrix lnsp_131" << std::endl;

    // Import matrix
    Matrix<double, StorageOrder::RowMajor> testMatrix(0, 0);
    testMatrix.reader(static_cast<std::string>("lnsp_131.mtx"));

    // Generate vector
    std::vector<double> vec(testMatrix.get_rows(), 0), res(testMatrix.get_rows(), 0);
    std::random_device seed;
    std::default_random_engine gen(seed());
    std::uniform_real_distribution<double> distr(-1., 1.);
    for (auto &val : vec)
    {
        val = distr(gen);
    }

    testMatrix.compress();
    res = testMatrix * vec;

    testMatrix.uncompress();
    res = testMatrix * vec;

    std::string filename = "execution_time.json";
    json data = read_json(filename);
    save_json(filename, data);
    std::cout << "Uncompressed format execution time: " << data["uncompressed_format_product"] << " nanoseconds."
              << std::endl;
    std::cout << "Compressed format execution time: " << data["compressed_format_product"] << " nanoseconds."
              << std::endl;
    // speedup
    if (data["compressed_format_product"] == 0)
    {
        std::cout << "Compressed format execution time is 0, cannot calculate speedup." << std::endl;
        return 1;
    }
    else
    {
        double speedup = static_cast<double>(data["uncompressed_format_product"]) / static_cast<double>(data["compressed_format_product"]);
        std::cout << "Speedup: " << speedup << std::endl;
    }
    return 0;
}
