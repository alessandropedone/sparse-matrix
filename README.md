[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/HlQKP7Zu)

# Sparse matrix 

test also the product with modified compressed format for SquareMatrix

how do i manage to compute products (views and SquareMatrix)

Parallelization with std algorithms (remove parallel implementations)

Modified compressed techniques (MSR, MSC)
- norm method issue (matrix pointer to squareMatrix object would call matrix::norm):
    - cannot be virtual
    - try to overcome the problem with dynamic_cast
    [See here for more details on the matter and proposed solution](https://chatgpt.com/share/680cae04-c850-800c-b63c-dece2a3d7728)
    - we could also try CTRP with abstract class
- add all methods to views
- SquareMatrix has a transpose and diagonal view
  
README.md (
    need of TBB, 
    we didn't implement COO since COOmap is clearly more efficient, 
    explain why parallel implementation of compress is slower, put a test about this)