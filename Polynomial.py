from IntegersModP import IntegersMod
from finitefield_functions import largest_index


class Polynomial:
    """Polynomial objects. Compatible with all classes that have defined
    addition, substraction, multipliaciton and division methods.
    """

    def __init__(self, coeffs: list):
        """
        Creates a polynomial of the form
        a_0 + a_1 x + a_2 x^2 + ...+ a_n x^n, 
        where indices are same as in list
        """
        self.coeffs = coeffs
        self.degree = largest_index(self.coeffs)

    def __str__(self):
        def to_superscript(n):
            superscript_digits = '⁰¹²³⁴⁵⁶⁷⁸⁹'
            return ''.join(superscript_digits[int(digit)] for digit in str(n))
        
        if all(coef == 0 for coef in self.coeffs):
            return "0"
        
        if isinstance(self.coeffs[0], IntegersMod):
            coeffs = [coef.value for coef in self.coeffs]
            mod_str = f" (mod {self.coeffs[0].prime})"
        else:
            coeffs = self.coeffs
            mod_str = ""
        
        polynomial_terms = []
        for degree in range(len(coeffs) - 1, -1, -1):
            coef = coeffs[degree]
            if coef == 0:
                continue
            
            coef_str = "" if abs(coef) == 1 and \
                degree != 0 else str(abs(coef))
            variable_str = "" if degree == 0 else "x"
            exponent_str = "" if degree <= 1 else to_superscript(degree)
            
            term = f"{coef_str}{variable_str}{exponent_str}"
            if coef < 0:
                term = f"-{term}"
            elif polynomial_terms:
                term = f"+{term}"
            
            polynomial_terms.append(term)
        
        return "".join(polynomial_terms) + mod_str
    
    def __repr__(self):
        return str(self)
    
    def __iter__(self):
        return iter(self.coeffs)

    def __add__(self, other):
        if other == 0:
            return self
        if not isinstance(other, Polynomial):
            raise TypeError("Both need to be polynomials")
        n = min(self.degree, other.degree) + 1 
        coeffs1 = self.coeffs
        coeffs2 = other.coeffs
        return Polynomial([coeffs1[i] + coeffs2[i] for i in range(n)]
                          + coeffs1[n:] + coeffs2[n:])
    
    def __radd__(self, other):
        return self.__add__(other)
    
    def __sub__(self, other):
        if other == 0:
            return self
        if not isinstance(other, Polynomial):
            raise TypeError("Both need to be polynomials")
        n = min(self.degree, other.degree) + 1 
        coeffs1 = self.coeffs
        coeffs2 = other.coeffs
        return Polynomial([coeffs1[i] - coeffs2[i] for i in range(n)]
                          + coeffs1[n:] + [-i for i in coeffs2[n:]])
    
    def __neg__(self):
        return Polynomial([-i for i in self.coeffs])

    def __mul__(self, other):
        if other == 0:
            return Polynomial([0])
        if other == 1:
            return self
        if not isinstance(other, Polynomial):
            raise TypeError("Both need to be polynomials")
        deg_self = self.degree + 1
        deg_other = other.degree + 1
        tmp_poly = [0]*(deg_self + deg_other)

        for i in range(deg_self):
            for j in range(deg_other):
                tmp_poly[i+j] += self.coeffs[i]*other.coeffs[j]
        return Polynomial(tmp_poly[:largest_index(tmp_poly) + 1])
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def division(self, other):
        """Euclidian algorithm for polynomials
        """

        if not isinstance(other, Polynomial):
            raise TypeError("Both need to be polynomials")
        if other == 0:
            raise ValueError("The divisor polynomial cannot be zero.")
        
        dividend = self.coeffs[::-1]
        divisor = other.coeffs[::-1]

        quotient = [0] * (self.degree - other.degree + 1)
        remainder = dividend[:]
        
        for i in range(self.degree - other.degree + 1):
            leading_coeff = remainder[i] / divisor[0]
            quotient[i] = leading_coeff
            
            for j in range(len(divisor)):
                remainder[i + j] -= leading_coeff * divisor[j]
        return Polynomial(quotient[::-1]), \
               Polynomial(remainder[:largest_index(remainder)+1][::-1])

    def __truediv__(self, other):
        quotient, remainder = self.division(other)
        if remainder != 0:
            raise ValueError("Other does not divide self")
        return quotient
    
    def __floordiv__(self, other):
        return self.division(other)[0]
    
    def __mod__(self, other):
        return self.division(other)[1]

    def __eq__(self, other):
        if other == 0:
            return all([i == 0 for i in self.coeffs])
        if not isinstance(other, Polynomial):
            return False
        return self - other == 0