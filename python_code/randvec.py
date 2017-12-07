import math 
import numpy as np
from random import random
from random import gauss

def random_vector(dims):
    """Return unifrom random vectors

    dims: dimensions
    """
    vec = [gauss(0,1) for i in range(dims)]
    mag = sum(x**2 for x in vec)**.5
    return [x/mag for x in vec]

def random_vector_inside_sphere():
    """Return uniform distribution inside a unit ball

    """
    y = random()**(1.0/3.0)   
    vec = random_vector(3) 
    return [x*y for x in vec]

def random_vector_inside_shell(inner,outer):
    """Return uniform distribution inside a shell

    """
    a = 0 
    b = 0 
    c = 0 
    r = 0
    while ( r < inner or r > outer ):
        a = (random()-0.5)*2*outer
        b = (random()-0.5)*2*outer
        c = (random()-0.5)*2*outer
        r = a**2 + b**2 + c**2
    return [a,b,c]
    
