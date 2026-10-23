"""Exact symbolic verification of the independently derived identities."""
import itertools,json,platform,resource,time
import sympy as S
from models import ROOT,X,K,x,y,z,a,b,c,d,e,f,model
start=time.monotonic();output={'sympy':S.__version__,'python':platform.python_version(),'networks':{}}
p,q,r=S.symbols('p q r',positive=True)
for index in [33449,239089]:
    m=model(index);Y,N,F,J,D,phi,P,power=(m[k] for k in ['Y','N','F','J','D','phi','P','power'])
    assert N.rank()==3 and not N.T.nullspace()
    assert all(t in [0,1] for t in Y) and all(t in [0,1] for t in N+Y)
    assert max(sum(Y[:,j]) for j in range(6))==2
    assert S.cancel(F[0].subs(phi))==0 and S.cancel(F[1].subs(phi))==0
    residual=S.cancel(F[2].subs(phi))
    assert S.cancel(residual+a*z**power*P/D**2)==0
    minor=J.extract([0,1],[0,1]).det()
    assert S.expand(minor-D)==0
    determinant=S.cancel(J.det().subs(phi))
    assert S.cancel(determinant-D*S.diff(residual,z))==0
    difference=S.cancel(determinant+a*z**power*S.diff(P,z)/D)
    assert S.cancel(difference/P-a*z**(power-1)*(2*z*S.diff(D,z)-power*D)/D**2)==0
    poly=S.Poly(P,z)
    if index==33449:
        assert S.expand(poly.LC()-c*e*(c*d-b*e-a*f))==0
        assert poly.TC()==b**2*d**2
        assert S.expand(P.subs({a:1,b:1,c:1,d:2,e:1,f:1})-(z**2+2*z+4))==0
        flux=S.Matrix([p,q,r,p+q,q+r,p-q+r])
        expected_minor=[p*(q+r),q*(2*p+r),r*(p+q+r)]
    else:
        assert poly.LC()==-3*c**2*e**2
        assert S.expand(poly.TC()-b*d*(c*d+b*e-a*f))==0
        flux=S.Matrix([p,p-q+r,q,p+q,r,q+r])
        expected_minor=[p*(p+q+r),r*(p+q),q*(p+2*r)]
    assert N*flux==S.zeros(3,1)
    Jflux=J.subs({k:flux[j]/m['monomials'][j] for j,k in enumerate(K)},simultaneous=True).applyfunc(S.simplify)
    minors=[]
    for ids,expected in zip([(0,1),(0,2),(1,2)],expected_minor):
        actual=S.factor(Jflux.extract(ids,ids).det())
        assert S.cancel(actual-expected/(X[ids[0]]*X[ids[1]]))==0
        minors.append(str(actual))
    coeff=Jflux.charpoly().all_coeffs();b1,b2,b3=coeff[1:]
    A,B,C=[-Jflux[i,i] for i in range(3)]
    M12,M13,M23=[Jflux.extract(ids,ids).det() for ids in [(0,1),(0,2),(1,2)]]
    cycles=Jflux[0,1]*Jflux[1,2]*Jflux[2,0]+Jflux[0,2]*Jflux[1,0]*Jflux[2,1]
    positive_H2=(A+B)*M12+(A+C)*M13+(B+C)*M23+2*A*B*C+cycles
    assert S.cancel(b1*b2-b3-positive_H2)==0
    rank_minor=next((list(ids),int(N[:,ids].det())) for ids in itertools.combinations(range(6),3) if N[:,ids].det())
    output['networks'][str(index)]={'mass_action':[str(t) for t in F],'rank':3,'full_rank_minor_columns_zero_based':rank_minor,'conservation_left_kernel':[], 'D':str(D),'admissibility':'z>0 and D(z)>0','x_of_z':str(phi[x]),'y_of_z':str(phi[y]),'P':str(P),'P_coefficients_descending':[str(t) for t in poly.all_coeffs()],'residual_formula':f'-a*z**{power}*P/D**2','determinant_at_root':f'-a*z**{power}*Pprime/D','fluxes':[str(t) for t in flux],'flux_condition':'p>0,q>0,r>0,p-q+r>0','flux_jacobian':[[str(t) for t in row] for row in Jflux.tolist()],'positive_principal_minors':minors,'H2_positive_identity':True,'origin_charpoly':str(S.factor(J.subs({x:0,y:0,z:0}).charpoly().as_expr())),'symbolic_checks':'PASS'}
output['wall_seconds']=time.monotonic()-start;output['ru_maxrss_bytes_macos']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
(ROOT/'independent-identities.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({k:{'P':v['P'],'det_at_root':v['determinant_at_root'],'minors':v['positive_principal_minors']} for k,v in output['networks'].items()},indent=2))
print('PASS: matrix-derived ODE, Cramer elimination, determinant chain rule, flux Jacobians and Hurwitz identity')
