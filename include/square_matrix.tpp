#ifndef SQUARE_MATRIX_TPP
#define SQUARE_MATRIX_TPP

#include "square_matrix.hpp"

namespace algebra
{
    
    /// @brief compress the matrix in modified format
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::compress_mod(){
        if(modified)
            return;
        
        //clear the modified compressed matrix
        compressed_format_mod.values.clear();
        compressed_format_mod.bind.clear();
        
        //reserve space for modified compressed structure
        const size_t nnz = this->get_nnz();
        compressed_format_mod.values.resize(nnz);
        compressed_format_mod.bind.resize(nnz);
        std::fill(compressed_format_mod.values.begin(), compressed_format_mod.values.end(), 0);
        std::fill(compressed_format_mod.bind.begin(), compressed_format_mod.bind.end(), 0);

        //to store correctly the pointers in the bind vector (that don't account for the diagonal elements)
        size_t off_diag_idx = 0;

        if(compressed){
            if constexpr (S == StorageOrder::ColumnMajor)
            {
                // iterate over columns of m
                for (size_t col_idx = 0; col_idx < this->cols; col_idx++)
                {
                    //set col pointer
                    compressed_format_mod.bind[col_idx] = off_diag_idx + this->rows;

                    // iterate over rows of m that are non-zero in the column "col" of m
                    size_t start = compressed_format.inner[col_idx];
                    size_t end = compressed_format.inner[col_idx + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        // get the row index of the non-zero element
                        size_t row_idx = compressed_format.outer[j];
    
                        // distinguish diagonal and off-diagonal elements
                        if(row_idx == col_idx){// diagonal element
                            compressed_format_mod.values[row_idx] = compressed_format.values[j];
                        }
                        else{// off-diagonal element
                            compressed_format_mod.values[this->rows + off_diag_idx] = compressed_format.values[j];
                            compressed_format_mod.bind[this->rows + off_diag_idx] = row_idx;
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
                    compressed_format_mod.bind[row_idx] = off_diag_idx + this->rows;
                    
                    // iterate over columns of m that are non-zero in the row "row" of m
                    size_t start = compressed_format.inner[row_idx];
                    size_t end = compressed_format.inner[row_idx + 1];
                    for (size_t j = start; j < end; j++)
                    {
                        // get the col index of the non-zero element
                        size_t col_idx = compressed_format.outer[j];

                        // distinguish diagonal and off-diagonal elements
                        if(row_idx == col_idx){// diagonal element
                            compressed_format_mod.values[row_idx] = compressed_format.values[j];
                        }
                        else{// off-diagonal element
                            compressed_format_mod.values[this->rows + off_diag_idx] = compressed_format.values[j];
                            compressed_format_mod.bind[this->rows + off_diag_idx] = col_idx;
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
                        compressed_format_mod.bind[it.first.col] = off_diag_idx + this->rows;
                        ++current_col;
                    }

                    if(it.first.col == it.first.row){// diagonal element
                        compressed_format_mod.values[it.first.row] = it.second;
                    }
                    else{// off-diagonal element
                        compressed_format_mod.values[this->rows + off_diag_idx] = it.second;
                        compressed_format_mod.bind[this->rows + off_diag_idx] = it.first.row;
                        ++off_diag_idx;
                    }
                }
                else{
                    size_t current_row = 0;
                    // set row pointers
                    if(current_row == it.first.row){
                        compressed_format_mod.bind[it.first.row] = off_diag_idx + this->rows;
                        ++current_row;
                    }
                    if(it.first.col == it.first.row){// diagonal element
                        compressed_format_mod.values[it.first.row] = it.second;
                    }
                    else{// off-diagonal element
                        compressed_format_mod.values[this->rows + off_diag_idx] = it.second;
                        compressed_format_mod.bind[this->rows + off_diag_idx] = it.first.col;
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
        if (modified)
        {
            std::cout << "Matrix is in modified compressed format, uncompressing..." << std::endl;
            // uncompress the matrix
            uncompress();
        }
        Matrix<T, S>::set(row, col, value);
    };

    /// @brief compress the matrix if it is in an uncompressed format
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::compress(){
        if(compressed)
            return;
        if(modified){
            //clear the compressed matrix
            this->compressed_format.inner.clear();
            this->compressed_format.outer.clear();
            this->compressed_format.values.clear();

            // reserve space for the compressed matrix
            if constexpr (S == StorageOrder::ColumnMajor){
                compressed_format.inner.resize(cols + 1);
            }
            else{
                compressed_format.inner.resize(rows + 1);
            }
            std::fill(compressed_format.inner.begin(), compressed_format.inner.end(), 0);
            
            // fill the compressed matrix
            size_t index = 0; //keeps track of nnz elements
            if constexpr (S == StorageOrder::ColumnMajor){
                for(size_t i = 0; i < this->rows; ++i){
                    bool flag = 0; // flag to check if the diagonal element has been inserted
                    size_t start = compressed_format_mod.bind[i];
                    size_t end = compressed_format_mod.bind[i + 1];
                    // handle last column
                    if (i+1 == this->rows){
                        end = compressed_format_mod.values.size() - 1;
                    }
                    for(size_t j = start; j < end; ++j){
                        size_t rowidx = compressed_format_mod.bind[j];
                        if(rowidx < i){
                            compressed_format_values.push_back(compressed_format_mod.values[j]);
                            compressed_format.outer.push_back(rowidx);
                        }
                        else{
                            if(!flag && compressed_format_mod.values[i]!= 0){
                                compressed_format.values.push_back(compressed_format_mod.values[i]);
                                compressed_format.outer.push_back(rowidx);
                                ++index;
                                flag = true;
                            }
                        compressed_format.values.push_back(compressed_format_mod.values[j]);
                        compressed_format.outer.push_back(rowidx);
                        }
                    }
                    index += end - start;
                    compressed_format.inner[i + 1] = index;
                }
            }
            else{
                for(size_t i = 0; i < this->rows; ++i){
                    bool flag = 0; // flag to check if the diagonal element has been inserted
                    size_t start = compressed_format_mod.bind[i];
                    size_t end = compressed_format_mod.bind[i + 1];
                    // handle last row
                    if (i+1 == this->rows){
                        end = compressed_format_mod.values.size() - 1;
                    }
                    for(size_t j = start; j < end; ++j){
                        size_t colidx = compressed_format_mod.bind[j];
                        if(colidx < i){
                            compressed_format_values.push_back(compressed_format_mod.values[j]);
                            compressed_format.outer.push_back(colidx);
                        }
                        else{
                            if(!flag && compressed_format_mod.values[i]!= 0){
                                compressed_format.values.push_back(compressed_format_mod.values[i]);
                                compressed_format.outer.push_back(colidx);
                                ++index;
                                flag = true;
                            }
                        compressed_format.values.push_back(compressed_format_mod.values[j]);
                        compressed_format.outer.push_back(colidx);
                        }
                    }
                    index += end - start;
                    compressed_format.inner[i + 1] = index;
                }
            }

            // clear the modified compressed matrix
            compressed_format_mod.values.clear();
            compressed_format_mod.bind.clear();

            // update the compressed flag
            this->compressed = true;
            return;
        }
        Matrix<T, S>::compress();
        return;
    };

    /// @brief compress the matrix in parallel if it is in an uncompressed format
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::compress_parallel(){

    };

    /// @brief uncompress the matrix if it is in a compressed format
    template <AddMulType T, StorageOrder S>
    void SquareMatrix<T, S>::uncompress(){
        if(modified){
            // clear the uncompressed format
            this->uncompressed_format.clear();

            // fill the uncompressed matrix
            if constexpr(S == StorageOrder::ColumnMajor){
                for (size_t col_idx = 0; col_idx < this->cols; ++col_idx){
                    // add diagonal element
                    if (mod_comp_format.values[col_idx] != 0){
                        this->uncompressed_format[{col_idx, col_idx}] = mod_comp_format.values[col_idx];
                    }
                    size_t start = mod_comp_format.bind[col_idx];
                    size_t end = mod_comp_format.bind[col_idx + 1];
                    if (col_idx + 1 == this->cols)
                        end = mod_comp_format.values.size() - 1;
                    for(size_t j = start; j < end; ++j){
                        size_t row_idx = mod_comp_format.bind[j];
                        this->uncompressed_format[{row_idx, col_idx}] = mod_comp_format.values[j];
                    }
                }
            }
            else{
                for (size_t row_idx = 0; row_idx < this->rows; ++row_idx){
                    // add diagonal element
                    if (mod_comp_format.values[row_idx] != 0){
                        this->uncompressed_format[{row_idx, row_idx}] = mod_comp_format.values[row_idx];
                    }
                    size_t start = mod_comp_format.bind[row_idx];
                    size_t end = mod_comp_format.bind[row_idx + 1];
                    if (row_idx + 1 == this->rows)
                        end = mod_comp_format.values.size() - 1;
                    for(size_t j = start; j < end; ++j){
                        size_t col_idx = mod_comp_format.bind[j];
                        this->uncompressed_format[{row_idx, col_idx}] = mod_comp_format.values[j];
                    }
                }
            }

            return;
        }
        Matrix<T, S>::uncompress();
        return;
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
        this->compressed_format_mod.values.clear();
        this->compressed_format_mod.bind.clear();
    };

};

#endif // SQUARE_MATRIX_TPP