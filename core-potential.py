import sympy
from sympy import symbols
sympy.init_printing()

g,l,m,b,t = symbols("gamma lambda mu b theta")
K = l + sympy.Rational(2,3)*m

C1 = (l-2*m)/(l+2*m)
C2 = K/(2*(l+2*m))

A1c = C1*t/6*((g**2) - 1)
A1s = C1*t/6*((g**2) + 2/C1)
A2s = -C2*t*(b**2)*(g**2)
B = t/3*(1-(g**2))

cerr = A1c
cett = A1c
cezz = B

cekk = cerr + cett + cezz

csrr = 2*(l+m)*A1c + l*B
cstt = 2*(l+m)*A1c + l*B
cszz = 2*l*A1c + (l+2*m)*B

cse = (cerr*csrr + cett*cstt + cezz*cszz)/2
ccw = (1+cekk)*csrr

corepotential = cse - ccw
corepotential = sympy.expand(corepotential)
