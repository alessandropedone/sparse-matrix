#ifndef SQUAREMATRIX_TPP
#define SQUAREMATRIX_TPP

#include "squarematrix.hpp"

namespace algebra
{
    
    /// @brief compress the matrix in modified format
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::mod_compress(){
        if(modified)
            return;
        
        //clear the modified compressed matrix
        mod_comp_format.values.clear();
        mod_comp_format.bind.clear();
        
        //reserve space for modified compressed structure
        const size_t nnz = this->get_nnz();
        mod_comp_format.values.resize(nnz);
        mod_comp_format.bind.resize(nnz);
        std::fill(mod_comp_format.values.begin(), mod_comp_format.values.end(), 0);
        std::fill(mod_comp_format.bind.begin(), mod_comp_format.bind.end(), 0);

        //to store correctly the pointers in the bind vector (that don't account for the diagonal elements)
        size_t off_diag_idx = 0;

        if(compressed){
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                // iterate over columns of m
                for (size_t col_idx = 0; col_idx < this->cols; col_idx++)
                {
                    //set col pointer
                    mod_comp_format.bind[col_idx] = off_diag_idx + this->rows;

                    // iterate over rows of m that are non-zero in the column "col" of m
                    size_t start = compressed_format.inner[col_idx];
                    size_t end = compressed_format.inner[col_idx + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        // get the row index of the non-zero element
                        size_t row_idx = compressed_format.outer[j];
    
                        // distinguish diagonal and off-diagonal elements
                        if(row_idx == col_idx){// diagonal element
                            mod_comp_format.values[row_idx] = compressed_format.values[j];
                        }
                        else{// off-diagonal element
                            mod_comp_format.values[this->rows + off_diag_idx] = compressed_format.values[j];
                            mod_comp_format.bind[this->rows + off_diag_idx] = row_idx;
                            ++off_diag_idx;
                        }
                    }
                }
            }
            else
            {
                // iterate over rows of m
                for (size_t row_idx = 0; row_idx < this->rows; row_idx++)
                {
                    //set row pointer
                    mod_comp_format.bind[row_idx] = off_diag_idx + this->rows;
                    
                    // iterate over columns of m that are non-zero in the row "row" of m
                    size_t start = compressed_format.inner[row_idx];
                    size_t end = compressed_format.inner[row_idx + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        // get the col index of the non-zero element
                        size_t col_idx = compressed_format.outer[j];

                        // distinguish diagonal and off-diagonal elements
                        if(row_idx == col_idx){// diagonal element
                            mod_comp_format.values[row_idx] = compressed_format.values[j];
                        }
                        else{// off-diagonal element
                            mod_comp_format.values[this->rows + off_diag_idx] = compressed_format.values[j];
                            mod_comp_format.bind[this->rows + off_diag_idx] = col_idx;
                            ++off_diag_idx;
                        }
                    }
                }
            }
    
            // clear the compressed matrix
            compressed_format.inner.clear();
            compressed_format.outer.clear();
            compressed_format.values.clear();
        }
        else{
            for (const auto &it : this->uncompressed_format)
            {
                if constexpr (S == StorageOrder::ColumnMajor){
                    size_t current_col = 0;
                    // set column pointers
                    if(current_col == it.first.col){
                        mod_comp_format.bind[it.first.col] = off_diag_idx + this->rows;
                        ++current_col;
                    }

                    if(it.first.col == it.first.row){// diagonal element
                        mod_comp_format.values[it.first.row] = it.second;
                    }
                    else{// off-diagonal element
                        mod_comp_format.values[this->rows + off_diag_idx] = it.second;
                        mod_comp_format.bind[this->rows + off_diag_idx] = it.first.row;
                        ++off_diag_idx;
                    }
                }
                else{
                    size_t current_row = 0;
                    // set row pointers
                    if(current_row == it.first.row){
                        mod_comp_format.bind[it.first.row] = off_diag_idx + this->rows;
                        ++current_row;
                    }
                    if(it.first.col == it.first.row){// diagonal element
                        mod_comp_format.values[it.first.row] = it.second;
                    }
                    else{// off-diagonal element
                        mod_comp_format.values[this->rows + off_diag_idx] = it.second;
                        mod_comp_format.bind[this->rows + off_diag_idx] = it.first.col;
                        ++off_diag_idx; }
                }
            // clear the uncompressed matrix
            uncompressed_format.clear();
            }
        }
        
        // update flags
        this->compressed = false;
        this->modified = true;
        return;
    };

    /// @brief set an element in the matrix (dynamic construction of the matrix)
    /// @param row row index
    /// @param col column index
    /// @param value value to set
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::set(size_t row, size_t col, const T &value){

    };

    /// @brief compress the matrix if it is in an uncompressed format
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::compress(){

    };

    /// @brief compress the matrix in parallel if it is in an uncompressed format
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::compress_parallel(){

    };

    /// @brief uncompress the matrix if it is in a compressed format
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::uncompress(){

    };

    /// @brief uncompress the matrix in parallel if it is in a compressed format
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::uncompress_parallel(){

    };

    /// @brief call operator() const version
    /// @param row row index
    /// @param col column index
    /// @return element at (row, col)
    template <AddMulType T, StorageOrder S>
    T SquareMatrix<T, S>::operator()(size_t row, size_t col) const{

    };


    /// @brief call operator() non-const version
    /// @param row row index
    /// @param col column index
    /// @return reference to the element at (row, col) with proxy (to avoid setting zero values)
   // template <AddMulType T, StorageOrder S>
   // Proxy<T> SquareMatrix<T, S>::operator()(size_t row, size_t col){

   // };

    /// @brief resize the matrix
    /// @param rows number of rows
    /// @param cols number of columns
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::resize_and_clear(size_t dim){
        this->rows = dim;
        this->cols = dim;
        this->compressed = false;
        this->modified = false
        this->uncompressed_format.clear();
        this->compressed_format.inner.clear();
        this->compressed_format.outer.clear();
        this->compressed_format.values.clear();
        this->mod_comp_format.values.clear();
        this->mod_comp_format.bind.clear();
    };

};

#endif // SQUAREMATRIX_TPP