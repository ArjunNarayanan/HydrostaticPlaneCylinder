import sympy
from sympy import symbols
sympy.init_printing()

R,lamda,mu,b,theta = symbols("R lamda mu b theta")
A1c,A1s,A2s,B = symbols("A1c A1s A2s B")

M = sympy.Matrix([[R,-R,-1/R,0],
                  [2*(lamda+mu),-2*(lamda+mu),2*mu/R**2,0],
                  [0,2*(lamda+mu),-2*mu/b**2,lamda],
                  [2*lamda*R**2,2*lamda*(b**2 - R**2),0,(lamda+2*mu)*b**2]])
K = lamda + sympy.Rational(2,3)*mu
R = sympy.Matrix([0,-K*theta,K*theta,K*theta*(b**2-R**2)])

from sympy import linsolve
coeffs = linsolve((M,R),A1c,A1s,A2s,B)

coefficients = coeffs.args[0]
