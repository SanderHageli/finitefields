from finitefield_functions import ALWAYS_REDUCE, isPrime
from IntegersModP import IntegersMod
from irred_poly_finder import modulo_method

class FiniteField:
    """The finite field with p^n elements, where p is a prime
    """

    def __init__(self, base_prime, degree):
        """
        Creates the finite field with base_prime**degree elements
        """
        if not isPrime(base_prime):
            raise ValueError(f"Base number {base_prime} must be prime")
        if degree != 0 and not isinstance(degree, int):
            raise ValueError(f"Degree {degree} must be postive integer")
        self.characteristic = base_prime
        self.degree = degree

    def __str__(self):
        return f"F({self.characteristic}^{self.degree})"

    def irred_poly(self):
        """Returns the polynomial used to define the relation between
        field elements
        """
        prime = self.characteristic
        irred_polys = modulo_method(self.degree, IntegersMod(prime))
        return irred_polys[-1]

    def size(self):
        """Returns the size of the field"""
        return self.characteristic ** self.degree
    
    def __call__(self, vector: list):
        return FieldElement(self.characteristic, self.degree, vector)
    
    def identity(self):
        """Returns the identity element of the field
        """
        return self([1])
    
    def zero(self):
        """Retruns the zero element of the field
        """
        return self([0])
        

class FieldElement(FiniteField):
    """
    Field elements are regarded as vectors over the base field F_{base_prime}
    """

    def __init__(self, base_prime, degree, vector):
        super().__init__(base_prime, degree)

        if isinstance(vector[0], int):
            base_field = IntegersMod(base_prime)
            vector = [base_field(i) for i in vector]
        self.vector = vector

    def __str__(self):
        self.reduce_element()
        def to_superscript(n):
            superscript_digits = '⁰¹²³⁴⁵⁶⁷⁸⁹'
            return ''.join(superscript_digits[int(digit)] for digit in str(n))
        
        if self == super().zero():
            return "0"
        
        coeffs = [i.value if i!=0 else 0 for i in self.vector]
        root_terms = []
        for degree in range(len(coeffs)):
            coef = coeffs[degree]
            if coef == 0:
                continue
            
            coef_str = "" if abs(coef) == 1 and \
                degree != 0 else str(abs(coef))
            variable_str = "" if degree == 0 else "ω"
            exponent_str = "" if degree <= 1 else to_superscript(degree)
            
            term = f"{coef_str}{variable_str}{exponent_str}"
            if coef < 0:
                term = f"-{term}"
            elif root_terms:
                term = f"+{term}"
            
            root_terms.append(term)
        
        return "".join(root_terms)

    def __call__(self, vector):
        return super().__call__(vector)
    
    def __getitem__(self, i):
        if i >= self.len():
            return 0
        return self.vector[i]

    def __add__(self, other):
        if other == 0:
            return self
        _compare_class(self, other)
        length = max(self.len(), other.len())
        return FieldElement(self.characteristic, self.degree,
                [(self[i] + other[i]) for i in range(length)]).reduce_element()
    
    def __radd__(self, other):
        return self.__add__(other)
    
    def __sub__(self, other):
        if other == 0:
            return self
        _compare_class(self, other)
        length = max(self.len(), other.len())
        return FieldElement(self.characteristic, self.degree,
                [(self[i] - other[i]) for i in range(length)]).reduce_element()
    
    def __mul__(self, other):
        if other == 0:
            return super().zero()
        if other == 1:
            return self
        _compare_class(self, other)
        tmp_vector = [0]*(self.len()+other.len())
        self.reduce_element()
        other.reduce_element()
        for i in range(self.len()):
            for j in range(other.len()):
                tmp_vector[i+j] += self[i]*other[j]
        return FieldElement(self.characteristic, self.degree,
                            tmp_vector).reduce_element()
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def inverse(self):
        """Returns the inverse of the field element
        """
        if self == 0:
            raise ZeroDivisionError("Zero has no inverse")
        pass

    def __truediv__(self, other):
        _compare_class(self, other)
        return self * other.inverse()

    def __neg__(self):
        return FieldElement(self.characteristic, self.degree,
                            [-i for i in self.vector])

    def __eq__(self, other):
        if other == 0:
            return all(i == 0 for i in self.vector)
        if not isinstance(other, FieldElement):
            return False
        return (self - other).reduce_element() == 0
    
    def len(self):
        """Returns the length of the vector representation
        """
        return len(self.vector)
    
    def reduce_element(self, reduce_status = ALWAYS_REDUCE):
        """Reduces the element to a representaion of length less than the
        degree of the field if the variable ALWAYS_REDUCE is True
        """
        if not reduce_status:
            return self
        
        vector_repr = self.vector
        deg = self.degree
        irred_poly = super().irred_poly()

        algebraic_element = [-i for i in irred_poly]
        del algebraic_element[-1]

        while len(vector_repr) > deg:
            for i in range(deg):
                vector_repr[len(vector_repr)-1-deg+i] += \
                            vector_repr[-1]*algebraic_element[i]
            del vector_repr[-1]
        return FieldElement(self.characteristic, deg, vector_repr)


def _compare_class(field_element_1: FieldElement,
                   field_element_2: FieldElement):
    if not (isinstance(field_element_1, FieldElement) and
            isinstance(field_element_2, FieldElement)):
        raise TypeError("Both must be of FieldElement type")
    if field_element_1.characteristic != field_element_2.characteristic \
        or field_element_1.degree != field_element_2.degree:
        raise ValueError("Elements are from different fields")