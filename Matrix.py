from random import randint

class Matrix:
    """Matrix objects. Compatible with all classes that have defined addition,
      substraction, multipliaciton and division methods.
    """

    def __init__(self, coeffs: list):
        correct_dimension = all(len(row)==len(coeffs[0]) for row in coeffs)
        if not correct_dimension:
            raise IndexError("Matrix dimension not valid")
        self.coeffs = coeffs
        self.rows = len(coeffs)
        self.columns = len(coeffs[0])

    def __str__(self):
        rows_str = []
        for row in self.coeffs:
            rows_str.append("[" +" ".join(f'{item:5}' for item in row) +
                            "]")
        return "\n".join(rows_str)
    
    def __repr__(self):
        return f'Mat({self.rows}x{self.columns})'

    def __iter__(self):
        return iter(self.coeffs)
    
    def __getitem__(self, i):
        return self.coeffs[i]

    def __add__(self, other):
        if other == 0:
            return self
        _compare_matrices(self, other)
        tmp_matrix = []
        for i in range(self.rows):
            tmp_matrix.append([self[i][j] + other[i][j] for j in
                               range(self.columns)])
        return Matrix(tmp_matrix)
    
    def __radd__(self, other):
        return self.__add__(other)
    
    def __sub__(self, other):
        if other == 0:
            return self
        _compare_matrices(self, other)
        tmp_matrix = []
        for i in range(self.rows):
            tmp_matrix.append([self[i][j] - other[i][j] for j in
                               range(self.columns)])
        return Matrix(tmp_matrix)
    
    def __neg__(self):
        tmp_matrix = []
        for row in self:
            tmp_matrix.append([-i for i in row])
        return Matrix(tmp_matrix)

    def __mul__(self, other):
        if other == 1:
            return self
        if not isinstance(other, Matrix):
            raise TypeError("Other not Matrix type")
        if self.columns != other.rows:
            raise IndexError("Matrices are incompatible")
        tmp_matrix = []
        for i in range(self.rows):
            row = []
            for j in range(other.columns):
                row.append(sum([self[i][k]*other[k][j] for k in
                                range(self.columns)]))
            tmp_matrix.append(row)
        return Matrix(tmp_matrix)
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def __eq__(self, other):
        if other == 0:
            is_row_zero = [all(i == 0 for i in row) for row in self]
            return all(is_row_zero)
        if type(other) != Matrix:
            return False
        return self - other == 0

    def transpose(self):
        """Returns the transpose of the matrix.

        Matrix remains Unaffected.
        """
        tmp_matrix = []
        for j in range(self.columns):
            tmp_matrix.append([self[i][j] for i in range(self.rows)])
        return Matrix(tmp_matrix)
    
    def copy(self):
        """Makes a copy of the matrix
        """
        return Matrix(self.coeffs)

    def det(self):
        """Computes the determinant of the matrix.
        
        Needs to be square matrix.
        """
        if self.rows != self.columns:
            raise ValueError("Must be square matrix")
        if self.rows == 1:
            return self.coeffs[0][0]

        column1 = [row[0] for row in self.coeffs]
        tmp_matrix = self.remove_column(0)
        sum = 0

        for i in range(len(column1)):
            reduced_matrix = tmp_matrix.remove_row(i)
            if i%2 == 0:
                sum += column1[i] * reduced_matrix.det()
            else:
                sum -= column1[i] * reduced_matrix.det()
        
        return sum

    def remove_row(self, i=None):
        """Returns a matrix with row i removed
        """
        rows = self.coeffs
        if i != None:
            rows = rows[:i] + rows[i+1:]
        return Matrix(rows)
    
    def remove_column(self, j=None):
        """Returns a matrix with column j removed
        """
        rows = self.coeffs
        if j != None:
            rows = [row[:j] + row[j+1:] for row in rows]
        return Matrix(rows)

    def row_swap(self, i, j):
        """Swaps row i and j of the matrix
        """
        row_i, row_j = self.coeffs[i], self.coeffs[j]
        self.coeffs[i], self.coeffs[j] = row_j, row_i
        return self

    def row_divide(self, i, scale):
        """Divides all items on row i by 'scale'
        """
        row_i = self.coeffs[i]
        self.coeffs[i] = [j/scale for j in row_i]
        return self

    def row_add(self, i, vector: list):
        """Adds elements in vector to their respective place in row i
        """
        if len(vector) != self.columns:
            raise ValueError(f"Vector needs to be of size {self.columns}")
        row_i = self.coeffs[i]
        self.coeffs[i] = [vector[j] + row_i[j] for j in range(len(vector))]
        return self

    def solve(self, b: list):
        """Solves Ax = b by Gaussian elimination
        """
        A = self.transpose()
        coeffs = A.coeffs
        A = Matrix(coeffs + [b]).transpose()

        row_reduction(A)
        
        for i in range(A.rows // 2):
            A.row_swap(i, A.rows-1 - i)
        A = A.transpose()

        for i in range((A.rows-1) // 2):
            A.row_swap(i, A.rows-2 - i)
        A = A.transpose()
        
        row_reduction(A)

        return A.transpose()[-1][::-1]
    
def identity_matrix(n: int) -> Matrix:
    """Returns n by n identity matrix
    """
    tmp_matrix = []
    for i in range(n):
        row = [0]*n
        row[i] = 1
        tmp_matrix.append(row)
    return Matrix(tmp_matrix)

def zero_matrix(n: int, m: int) -> Matrix:
    """Returns n by m zero matrix
    """
    tmp_matrix = []
    for i in range(n):
        tmp_matrix.append([0]*m)
    return Matrix(tmp_matrix)

def random_matrix(n: int, m: int, a=-5,b=5) -> Matrix:
    """Generates a random n*m matrix with integer coefficients
    """
    tmp_matrix = []
    for i in range(n):
        row = [randint(a,b) for j in range(m)]
        tmp_matrix.append(row)
    return Matrix(tmp_matrix)

def row_reduction(A: Matrix):
    """Preforms Gaussian elimination on the matrix
    """
    n = min(A.rows, A.columns)
    for i in range(n):
        scaling_coeff = A[i][i]
        if scaling_coeff == 0:
            column = A.transpose()[i]
            k = max([column.index(l) for l in column if l != 0])
                #finds the largest non-zero element in column i
            A.row_swap(i, k)
            scaling_coeff = A[i][i]
        A.row_divide(i, scaling_coeff)
        row_to_add = A[i]

        for j in range(i + 1, n):
            leading_coeff = A[j][i]
            A.row_add(j, [-leading_coeff*k for k in row_to_add])

def _compare_matrices(A: Matrix, B: Matrix):
    """Checks if matrices A and B are of the same size
    """
    if A.rows != B.rows or A.columns != B.columns:
        raise ValueError("Both matrcies must be of same size")