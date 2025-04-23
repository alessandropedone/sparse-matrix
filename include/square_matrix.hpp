#ifndef SQUARE_MATRIX_HPP
#define SQUARE_MATRIX_HPP

#include "storage.hpp"
#include "matrix.hpp"
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
            SquareMatrix(int size) : Matrix<T, S>(size, size) : rows(size), cols(size) {
                this->compressed = false;
                this->modified = false;
            };

            /// @brief default destructor
            virtual ~SquareMatrix() override = default;

            /// @brief check if the matrix is in a compressed format
            /// @return true if the matrix is compressed, false otherwise            bool is_compressed() const { return compressed; };
            bool is_modified() const { return modified; };

            /// @brief compress the matrix in modified format
            void compress_mod();

            /// @brief set an element in the matrix (dynamic construction of the matrix)
            /// @param row row index
            /// @param col column index
            /// @param value value to set
            virtual void set(size_t row, size_t col, const T &value) override;

            /// @brief compress the matrix if it is in an uncompressed format
            virtual void compress() override;

            /// @brief compress the matrix in parallel if it is in an uncompressed format
            virtual void compress_parallel() override;

            /// @brief uncompress the matrix if it is in a compressed format
            virtual void uncompress() override;

            /// @brief uncompress the matrix in parallel if it is in a compressed format
            virtual void uncompress_parallel() override;

            /// @brief call operator() const version
            /// @param row row index
            /// @param col column index
            /// @return element at (row, col)
            virtual T operator()(size_t row, size_t col) const override;


            /// @brief call operator() non-const version
            /// @param row row index
            /// @param col column index
            /// @return reference to the element at (row, col) with proxy (to avoid setting zero values)
          //  virtual Proxy<T> operator()(size_t row, size_t col) override;

            /// @brief resize the matrix
            /// @param rows number of rows
            /// @param cols number of columns
            // Overload of Matrix<T, S>::resize_and_clear
            void resize_and_clear(size_t dim);

            /// @brief get the number of non-zero elements
            /// @return number of non-zero elements
            virtual size_t get_nnz() override const{
                if (modified)
                {
                    return mod_comp_format.values.size(); //this doesn't account for zeros in the diagonal
                }
                else
                {
                    return Matrix<T, S>::get_nnz();
                }
            };

        private:
            bool modified = false; // flag to check if the matrix is in modified compressed format

            // storage for the matrix
            ModifiedCompressedStorage<T> compressed_format_mod; // MSR or MSC format
        };
        
};

#include "square_matrix.tpp"

#endif // SQUARE_MATRIX_HPP