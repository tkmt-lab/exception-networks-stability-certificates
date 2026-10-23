"""Independent models from the allowed matrices; no producer-code imports."""
import json
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parent
x,y,z=S.symbols('x y z',positive=True)
a,b,c,d,e,f=S.symbols('a b c d e f',positive=True)
K=(a,b,c,d,e,f);X=(x,y,z)

def model(index):
    data=next(n for n in json.loads((ROOT/'input/networks.json').read_text())['networks'] if n['candidate_index']==index)
    Y=S.Matrix(data['Y']);N=S.Matrix(data['N'])
    monomials=S.Matrix([S.prod(X[i]**Y[i,j] for i in range(3)) for j in range(6)])
    F=N*S.diag(*K)*monomials
    if index==33449:
        expected=S.Matrix([a*z+b*y-d*x,-b*y-c*y*z+e*x*z,-a*z+b*y-c*y*z+f*x*y])
        D=d*b+(d*c-b*e)*z
        phi={x:a*z*(b+c*z)/D,y:a*e*z**2/D}
        P=S.expand(D**2-e*z*(b-c*z)*D-a*e*f*z**2*(b+c*z))
        power=1
    elif index==239089:
        expected=S.Matrix([a*z+c*y*z-d*x,a*z-b*y-c*y*z+e*x*z,-c*y*z-e*x*z+f*x*y])
        D=d*(b+c*z)-c*e*z**2
        phi={x:a*z*(b+2*c*z)/D,y:a*z*(d+e*z)/D}
        P=S.expand(D*(c*d+e*b+3*c*e*z)-a*f*(b+2*c*z)*(d+e*z))
        power=2
    else:raise ValueError(index)
    assert F==expected
    return {'data':data,'Y':Y,'N':N,'monomials':monomials,'F':F,'J':F.jacobian(X),'D':D,'phi':phi,'P':P,'power':power}
