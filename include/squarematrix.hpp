#ifndef SQUAREMATRIX_HPP
#define SQUAREMATRIX_HPP

#include "matrix.hpp"
#include "storage.hpp"
#include "proxy.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <execution>
#include <cmath>

namespace algebra
{
    /// @brief Square Matrix class, child of Matrix class
    /// @tparam T type of the matrix elements
    /// @tparam S storage order of the matrix (RowMajor or ColumnMajor)
    template <AddMulType T, StorageOrder S = StorageOrder::RowMajor>
    class SquareMatrix : public Matrix<T, S>{
        public:
            // delete default constructor
            SquareMatrix() = delete;

            /// @brief constructor with size
            /// @param size number of rows and columns
            SquareMatrix(int size) : Matrix<T, S>(size, size) {
                this->rows = size;
                this->cols = size;
                this->compressed = false;
            };

        private:
            ModifiedCompressedStorage<T> mod_comp_format; // MSR or MSC format
        };
        
};

#include "squarematrix.tpp"

#endif // SQUAREMATRIX_HPP