#ifndef STORAGE_HPP
#define STORAGE_HPP

#include <utility>
#include <cstddef>
#include <vector>
#include <map>
#include <iostream>
#include <concepts>
#include <complex>

namespace algebra
{
    /// @brief type of storage order
    enum class StorageOrder
    {
        RowMajor,
        ColumnMajor
    };

    // Primary template: false for all types
    template <typename T>
    struct is_complex : std::false_type
    {
    };

    // Specialization: true for std::complex<U>
    template <typename U>
    struct is_complex<std::complex<U>> : std::true_type
    {
    };

    template <typename T>
    struct AbsReturnType
    {
        using type = T;
    };

    template <typename T>
    struct AbsReturnType<std::complex<T>>
    {
        using type = T;
    };

    template <typename T>
    using AbsReturnType_t = typename AbsReturnType<T>::type;

    template <typename T>
    concept AddMulType = requires(T a, T b) {
        { a + b } -> std::convertible_to<T>;
        { a * b } -> std::convertible_to<T>;
        { std::abs(a) } -> std::convertible_to<AbsReturnType_t<T>>;
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

    /// @brief matrix storage in modified compressed format
    /// @tparam T type of the matrix elements
    template <AddMulType T>
    struct ModifiedCompressedStorage
    {
        // let nnz = number of non-zero elements, considering the whole principal diagonal NON-zero
        std::vector<T> values;
        // from 0 to n-1-> diagonal elements
        // from n to nnz-1 -> off-diagonal elements in row or column major order

        std::vector<size_t> bind;
        // from 0 to n-1 -> row or column pointer
        //(cumulative sum of nnz that are OFF the diagonal up to that row/col + size of matrix(first n elements are the diagonal ones))
        // from n to nnz - 1 -> column or row index of the off-diagonal elements
    };

    /// @brief alias for the index of the matrix
    struct Index
    {
        size_t row;
        size_t col;
    };

    /// @brief comparison operators for the index if the order is RowMajor
    struct RowMajor
    {
        bool operator()(const Index &a, const Index &b) const
        {
            return (a.row < b.row) or (a.row == b.row and a.col < b.col);
        }
    };

    /// @brief comparison operators for the index if the order is ColumnMajor
    struct ColMajor
    {
        bool operator()(const Index &a, const Index &b) const
        {
            return (a.col < b.col) or (a.col == b.col and a.row < b.row);
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