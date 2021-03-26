corepotential = sympy.simplify(corepotential)

serr = A1s - A2s/(b**2*g**2) - t/3
sett = A1s + A2s/(b**2*g**2) - t/3
sezz = B

ssrr = 2*(l+m)*A1s - 2*m*A2s/(b**2*g**2) + l*B - K*t
sstt = 2*(l+m)*A1s + 2*m*A2s/(b**2*g**2) + l*B - K*t
sszz = 2*l*A1s + (l+2*m)*B - K*t

sse = serr*ssrr + 