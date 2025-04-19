#ifndef STORAGE_HPP
#define STORAGE_HPP

#include <utility>
#include <cstddef>
#include <vector>
#include <map>
#include <iostream>
#include <concepts>
namespace algebra
{
    /// @brief type of storage order
    enum class StorageOrder
    {
        RowMajor,
        ColumnMajor
    };

    template <typename T>
    concept AddMulType = requires(T a, T b) {
        { a + b } -> std::convertible_to<T>;
        { a * b } -> std::convertible_to<T>;
        { abs(a) } -> std::convertible_to<T>;
    };

    /// @brief matrix storage in compressed format
    /// @tparam T type of the matrix elements
    template <AddMulType T>
    struct CompressedStorage
    {
        std::vector<size_t> inner; // Starting index for each row (for CSR) or column (for CSC)
        std::vector<size_t> outer; // Column (for CSR) or row (for CSC) indices of non-zero elements
        std::vector<T> values;     // Non-zero values
    };

    /// @brief alias for the index of the matrix
    struct Index {
        size_t row;
        size_t col;
    };

    /// @brief comparison operators for the index if the order is RowMajor
    struct RowMajor
    {
        bool operator()(const Index &a, const Index &b) const
        {
            return (a.row < b.row) || (a.row == b.row && a.col < b.col);
        }
    };

    /// @brief comparison operators for the index if the order is ColumnMajor
    struct ColMajor
    {
        bool operator()(const Index &a, const Index &b) const
        {
            return (a.col < b.col) || (a.col == b.col && a.row < b.row);
        }
    };

    /// @brief selector for the comparator based on the storage order
    /// @tparam S storage order
    template <StorageOrder S>
    struct ComparatorSelector;

    /// @brief specialization for RowMajor
    template <>
    struct ComparatorSelector<StorageOrder::RowMajor>
    {
        using type = RowMajor;
    };

    /// @brief specialization for ColumnMajor
    template <>
    struct ComparatorSelector<StorageOrder::ColumnMajor>
    {
        using type = ColMajor;
    };

    /// @brief matrix storage in uncompressed format
    /// @tparam T type of the matrix elements
    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    using UncompressedStorage = std::map<Index, T, typename ComparatorSelector<S>::type>;

}
#endif // STORAGE_HPP