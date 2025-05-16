/**
 * @file proxy.hpp
 * @author 
 * @brief Defines the Proxy class for sparse matrix element access and modification.
 * 
 * This file contains the declaration and implementation of the algebra::Proxy class template,
 * which provides a proxy object for accessing and modifying elements in a sparse matrix.
 * The Proxy enforces sparse storage rules by ensuring that zero values are not explicitly stored
 * in the underlying uncompressed storage format. Assignments and arithmetic operations are
 * intercepted to insert, update, or erase elements as appropriate.
 * 
 * @copyright
 * Copyright (c) 2024
 * 
 * @see algebra::Proxy
 * @see algebra::UncompressedStorage
 */
#ifndef PROXY_HPP
#define PROXY_HPP

#include "storage.hpp"

namespace algebra
{
    /**
     * @brief Proxy class for matrix elements that enforces sparse storage rules.
     * 
     * This proxy wraps access to a matrix element and ensures that zero values are not explicitly stored
     * in the underlying uncompressed storage. Assignments and arithmetic operations are intercepted:
     * - Assigning zero removes the element from storage.
     * - Assigning a nonzero value inserts or updates the element.
     * - Addition and subtraction update or erase the element as appropriate.
     * 
     * @tparam T Type of the matrix elements.
     * @tparam S Storage order of the matrix (RowMajor or ColumnMajor).
     * 
     * @note Used internally to prevent explicit storage of zero values in sparse matrices.
     */
    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class Proxy
    {
    private:
        UncompressedStorage<T, S> &uncompressed_format; /// reference to the uncompressed storage
        size_t row; /// row index
        size_t col; /// column index

    public:

        /// @brief Constructor
        /// @param uncompressed_format
        /// @param row 
        /// @param col 
        Proxy(UncompressedStorage<T, S> &uncompressed_format, size_t row, size_t col)
            : uncompressed_format(uncompressed_format), row(row), col(col)
        {
        }
        
        /// @brief implicit conversion operator
        /// @return the value of the matrix at (row, col)
        operator T() const
        {
            // find the value in the uncompressed format
            auto it = uncompressed_format.find({row, col});
            if (it == uncompressed_format.end())
            {
                // if the value is not found, return 0
                return T(0);
            }
            // if the value is found, return it
            return uncompressed_format[{row, col}];
        }

        /// @brief assignment operator
        /// @param val the value to assign
        /// @return the proxy
        /// @note this operator is used to assign a value to the matrix at (row, col)
        /// @note if the value is 0, the value is erased from the matrix
        /// @note if the value is not 0, the value is set in the matrix
        Proxy &operator=(T const &val)
        {
            if (val == T(0))
            {
                // erase the value
                uncompressed_format.erase({row, col});
            }
            else
            {
                // set the value
                uncompressed_format[{row, col}] = val;
            }
            return *this;
        }

        
        /// @brief addition operator
        /// @param val the value to add
        /// @return the proxy
        /// @note this operator is used to add a value to the matrix at (row, col)
        /// @note if the value is 0, the value is erased from the matrix
        /// @note if the value is not 0, the value is set in the matrix
        Proxy &operator+=(T const &val)
        {
            T newv = uncompressed_format[{row, col}] + val;
            if (newv == T(0))
            {
                // erase the value
                uncompressed_format.erase({row, col});
            }
            else
            {
                // set the value
                uncompressed_format[{row, col}] = newv;
            }
            return *this;
        }

        /// @brief subtraction operator
        /// @param val the value to subtract
        /// @return the proxy
        /// @note this operator is used to subtract a value from the matrix at (row, col)
        /// @note if the value is 0, the value is erased from the matrix
        /// @note if the value is not 0, the value is set in the matrix
        Proxy &operator-=(T const &val)
        {
            T newv = uncompressed_format[{row, col}] - val;
            if (newv == T(0))
            {
                // erase the value
                uncompressed_format.erase({row, col});
            }
            else
            {
                // set the value
                uncompressed_format[{row, col}] = newv;
            }
            return *this;
        }
    };
};

#endif // PROXY_HPP