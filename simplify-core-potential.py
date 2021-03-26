import sympy
from sympy import symbols
sympy.init_printing()

g,l,m,b,t = symbols("gamma lambda mu b theta")
K = l + sympy.Rational(2,3)*m

C1 = (l-2*m)/(l+2*m)

t1 = sympy.simplify(-l/3+(l+m)/3*C1)

t2 = -(l-2*m)/18 + (l+m)*C1/9 - (l+m)*(C1**2)/18
t2 = sympy.simplify(t2)

t3 = (l-2*m)/9 + (l+m)*(C1**2)/9
t3 = sympy.simplify(t3)