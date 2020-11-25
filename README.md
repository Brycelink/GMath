# GMath
A simple Math Library for Matrices, vectors and Quartenions.
GMath uses templates on its classes to allow use of different types for its Matrices and vectors.
It's default type however, is floating point values.

# Usage
GMath only requires you to utilize the GMath.h header file to access all other resources of GMath.
As noted in the GMath.h header file, one may modify the GMath.h file to suite the naming convention
of the type definitions. The existing format specifies the original concept when the library was being
built.

## Recommended Modifications
### GMATH_TYPE_PRECISION (GMath/GMath_defs.h)
Can be modified to change the default type of all Matrices, Quartenions and Vectors. It is by default set to float.
### MatrixC, MatrixSqC (GMath.h)
Refers to Column Major Matrices
### MatrixR, MatrixSqR (GMath.h)
Refers to Row Major Matrices
### Matrix, MatrixSq (GMath.h)
Selects between Row Major and Column Major as default. NOTE: Ensure both Matrix and MatrixSq use same Major.
### Vector, Vector4, Vector3, Vector2 (GMath.h)
Refers to Vectors and their type
### Quartenion (GMath.h)
Refers to Quartenions and their type

#Example
```c
GMath::MatrixSq<4> m = {1, 2, 3, 4, 5, 6, 7, 8, 553, 13, 22, 16, 17, 9, 76, 33};
//GMath::Matrix< 4, 4> m = {1, 2, 3, 4, 5, 6, 7, 8, 553, 13, 22, 16, 17, 9, 76, 33};

GMath::Matrix<4> n(1);
GMath::Matrix<4> result;

result = m.Multiply(n);
//result = m*n;

result = m.Add(n);
//result = m+n;

result = m.Subtract(n);
//result = m-n;
```

# License
GMath is licensed under the MIT open-source License [https://opensource.org/licenses/MIT]
