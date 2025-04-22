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
        size_t nnz = this->get_nnz();
        mod_comp_format.values.resize(nnz);
        mod_comp_format.bind.resize(nnz);
        std::fill(mod_comp_format.values.begin(), mod_comp_format.values.end(), 0);

        if(compressed){
            
            // clear the compressed matrix
            compressed_format.inner.clear();
            compressed_format.outer.clear();
            compressed_format.values.clear();
        }
        else{

            // clear the uncompressed matrix
            uncompressed_format.clear();
        }

        
        // update flags
        compressed = false;
        modified = true;
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
    template <AddMulType T, StorageOrder S>
    Proxy<T> SquareMatrix<T, S>::operator()(size_t row, size_t col){

    };

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