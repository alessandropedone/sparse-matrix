#ifndef PROXY_HPP
#define PROXY_HPP

#include "storage.hpp"

namespace algebra
{
    /// @brief A little proxy that “wraps” a T& and checks assignments and sums to make sure that the value is not zero
    /// @tparam T type of the matrix elements
    /// @tparam S storage order of the matrix (RowMajor or ColumnMajor)
    /// @note this is used to prevent zero values in the matrix
    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class Proxy
    {
    private:
        UncompressedStorage<T, S> &uncompressed_format;
        size_t row;
        size_t col;

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

        /// @brief friend abs function
        /// @param proxy the proxy object
        /// @return the absolute value of the matrix at (row, col)
        friend AbsReturnType_t<T> abs(const Proxy &proxy)
        {
            // find the value in the uncompressed format
            auto it = proxy.uncompressed_format.find({proxy.row, proxy.col});
            if (it == proxy.uncompressed_format.end())
            {
                // if the value is not found, return 0
                return AbsReturnType_t<T>(0);
            }
            // if the value is found, return it
            return std::abs(proxy.uncompressed_format[{proxy.row, proxy.col}]);
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
}

#endif // PROXY_HPP