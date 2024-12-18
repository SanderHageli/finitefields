from finitefield_functions import isPrime

class IntegersMod:
    """Implementation of the field of integers modulo p, for some prime p
    """

    def __init__(self, prime):
        if not isPrime(prime):
            raise ValueError(f"Base number {prime} must be prime")
        self.characteristic = prime
    
    def __call__(self,value):
        return IntegersModElement(value, self.characteristic)

    def elements(self):
        """Returns a list of all field elements
        """
        return [IntegersModElement(i, self.characteristic) for i in
                                    range(self.characteristic)]
    
    def identity(self):
        """Returns the identity element of the field
        """
        return self(1)
    
    def zero(self):
        """Retruns the zero element of the field
        """
        return self(0)

class IntegersModElement(IntegersMod):
    """Field elements of the field of integers modulo p
    """

    def __init__(self, value, prime):
        super().__init__(prime)
        self.value = value % prime
        self.characteristic = prime

    def __str__(self):
        return f'{self.value} (mod {self.characteristic})'
    
    def __repr__(self):
        return str(self)
    
    def __format__(self, format_spec):
        if isinstance(format_spec, int):
            l = len(str(self.value))

            return " "*(format_spec-l) + str(self)
        return str(self)

    def __add__(self, other):
        if other == 0:
            return self
        _compare_class(self, other)
        return IntegersModElement((self.value+other.value) 
                                  % self.characteristic, self.characteristic)
    
    def __radd__(self, other):
        return self.__add__(other)
    
    def __sub__(self, other):
        if other == 0:
            return self
        _compare_class(self, other)
        return IntegersModElement((self.value-other.value)
                                  % self.characteristic, self.characteristic)
    
    def __neg__(self):
        return IntegersModElement(-self.value % self.characteristic,
                                  self.characteristic)
    
    def __mul__(self, other):
        if other == 0:
            return IntegersModElement(0, self.characteristic)
        if other == 1:
            return self
        _compare_class(self, other)
        return IntegersModElement((self.value*other.value)
                                  % self.characteristic, self.characteristic)
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        if not isinstance(other, int):
            raise TypeError("Other needs to be int")
        return IntegersModElement(pow(self.value, other, self.characteristic),
                                  self.characteristic)

    def inverse(self):
        """Computes the inverse of the field element

        This is done by computing a^(p-2)
        """
        if self == 0:
            raise ZeroDivisionError("Zero has no inverse")
        return self ** (self.characteristic-2)
    
    def __truediv__(self, other):
        _compare_class(self, other)
        return self * other.inverse()

    def __eq__(self, other):
        if other == 0:
            return self.value == 0
        if other == 1:
            return self.value == 1
        if not isinstance(other, IntegersModElement):
            return False
        return self - other == 0

def _compare_class(element1: IntegersModElement,
                   element2: IntegersModElement):
    if not (isinstance(element1, IntegersModElement) and
            isinstance(element2, IntegersModElement)):
        raise TypeError(f"Both must be IntegersModPElement type")
    if element1.characteristic != element2.characteristic:
        raise ValueError("Primes must be equal")