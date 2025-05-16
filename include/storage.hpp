/**
 * @file storage.hpp
 * @author 
 * @brief Defines storage formats and utilities for sparse matrices.
 * 
 * This header provides various storage structures and utilities for representing
 * sparse matrices in different formats, including compressed, modified compressed,
 * and uncompressed (map-based) storage. It also defines concepts and type traits
 * for handling numeric and complex types, as well as storage order (row-major or column-major).
 * 
 * @details
 * The main components of this file are:
 * - @ref algebra::StorageOrder : Enum for specifying matrix storage order.
 * - @ref algebra::is_complex : Type trait to detect std::complex types.
 * - @ref algebra::AbsReturnType : Type trait to determine the return type of std::abs.
 * - @ref algebra::AddMulType : Concept for types supporting addition, multiplication, and absolute value.
 * - @ref algebra::CompressedStorage : Structure for compressed sparse matrix storage (CSR/CSC).
 * - @ref algebra::ModifiedCompressedStorage : Structure for modified compressed storage with explicit diagonal.
 * - @ref algebra::Index : Struct representing a matrix index (row, col).
 * - @ref algebra::RowMajor, @ref algebra::ColMajor : Comparators for index ordering.
 * - @ref algebra::ComparatorSelector : Selector for comparator based on storage order.
 * - @ref algebra::UncompressedStorage : Alias for map-based sparse matrix storage.
 * 
 * @copyright
 * Copyright (c) 
 */
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
    /**
     * @enum StorageOrder
     * @brief Enum class to specify the storage order of a matrix.
     *
     * This enum is used to indicate whether the matrix is stored in row-major or column-major order.
     * - RowMajor: Matrix elements are stored row by row.
     * - ColumnMajor: Matrix elements are stored column by column.
     */
    enum class StorageOrder
    {
        RowMajor,
        ColumnMajor
    };

    /// @brief check if the type is a complex number
    /// @tparam T type to check
    /// @note primary template is false for all types
    template <typename T>
    struct is_complex : std::false_type
    {
    };

    /// @brief specialization for std::complex<T>
    /// @tparam U type of the complex number
    /// @note specialization for std::complex<U> is true
    template <typename U>
    struct is_complex<std::complex<U>> : std::true_type
    {
    };

    /// @brief return type of the absolute value
    /// @tparam T type to check
    template <typename T>
    struct AbsReturnType
    {
        using type = T;
    };

    /// @brief specialization for std::complex<T>
    /// @tparam T type of the complex number
    template <typename T>
    struct AbsReturnType<std::complex<T>>
    {
        using type = T;
    };

    /// @brief alias for the return type of the absolute value
    /// @tparam T type to check
    template <typename T>
    using AbsReturnType_t = typename AbsReturnType<T>::type;

    /// @brief concept to check if the type is a number that supports addition, multiplication and absolute value
    /// @tparam T type to check
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

    /// @brief specialization for ColMajor
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