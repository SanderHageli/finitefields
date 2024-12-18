from math import sqrt
import os

from IntegersModP import IntegersMod
from Polynomial import Polynomial
from finitefield_functions import IRREDUCIBLE_POLYS_PATH


def gen_all_polys(degree, field):
    """
    Generates all monic polynomials of degree 'degree' over the given field
    """
    all_polys = []
    field_elements = field.elements()

    def set_coeffs(poly, deg):
        if deg == degree:
            all_polys.append(poly+[field.identity()])
            return
        for i in field_elements:
            set_coeffs(poly+[i], deg+1)

    set_coeffs([],0)
    return [Polynomial(p) for p in all_polys]

def write_poly_to_file(polys, file_folder = IRREDUCIBLE_POLYS_PATH):
    prime = (polys[0].coeffs)[0].characteristic
    poly_deg = polys[0].degree
    path = f'./{file_folder}/irred_polys_{prime}.txt'
    poly_line = f'{poly_deg}:'+";".join([",".join([str(i.value) for i in p.coeffs]) for p in polys])+"\n"

    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path, "r") as file:
            for line in file:
                pass
            last_written_deg = int(line[0])

        if poly_deg > last_written_deg:
            with open(path, "a") as file:
                file.write(poly_line)
    else:
        with open(path, "w") as file:
            file.write(poly_line)

def read_poly_from_file(prime, degree, file_folder = IRREDUCIBLE_POLYS_PATH):
    path = f'./{file_folder}/irred_polys_{prime}.txt'
    irred_polys = []
    F = IntegersMod(prime)
    i = 0

    if os.path.exists(path):
        with open(path, "r") as file:
            for line in file:
                i += 1
                if i == degree:                 
                    deg_k_polys_str = line.split(":")[1]
                    deg_k_polys_list = deg_k_polys_str.split(";")
                    irred_polys = [Polynomial([F(int(j)) for j in p.split(",")]) for p in deg_k_polys_list]
                if i > degree:
                    break

    return irred_polys, i
    

def modulo_method(degree: int, field):
    """
    Returns a list of all irreducible polynomials up to degree 'degree',
    by checking, for all polynomials p(x) of degree 'degree', if
    q(x) divides p(x) for any irreducible polynomial q(x) of lower degree
    """

    saved_polys = read_poly_from_file(field.characteristic, degree)[0]
    irreducible_polys = saved_polys
    lower_deg_irred_polys = []

    if degree > 1:
        lower_deg_irred_polys = modulo_method(degree-1, field)

    if saved_polys == []: #If there are no polys of spesified degree in file
        deg_k_polys = gen_all_polys(degree, field)
        if degree == 1:
            write_poly_to_file(deg_k_polys)
            return deg_k_polys
        
        lower_deg_irred_polys = modulo_method(degree-1, field) #kanskje noe redundencies her med tanke på cutoff-en

        cutoff = int(sqrt(degree))
        polys_to_check = [p for p in lower_deg_irred_polys if p.degree <= cutoff]
        
        irreducible_polys = []
        for p in deg_k_polys:
            if all(p % q != 0 for q in polys_to_check): 
                irreducible_polys.append(p)
        
        write_poly_to_file(irreducible_polys)

    return lower_deg_irred_polys + irreducible_polys

def sieve_element_method(degree: int, field):
    # Bare om det går an å sile med hensyn til indeksering. Tar jo O(n) å fjerne p fra deg_k_polys

    deg_k_polys = gen_all_polys(degree, field)
    if degree == 1:
        return deg_k_polys
    
    lower_deg_irred_polys = sieve_element_method(degree-1, field)

    cutoff = int(sqrt(degree))
    polys_to_check = [p for p in lower_deg_irred_polys if p.degree <= cutoff]
    
    i = 0
    k = len(deg_k_polys)
    while i < k:
        p = deg_k_polys[i]
        is_irred = True
        for q in polys_to_check:
            if p%q==0:
                is_irred = False
                while True:
                    try:
                        deg_k_polys.remove(p)
                        k -= 1
                        p += q
                    except:
                        break                    
        if is_irred:
            i+=1
        
    return lower_deg_irred_polys + deg_k_polys