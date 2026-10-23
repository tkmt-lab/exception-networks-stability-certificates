"""Package reconciliation adapted from the post-freeze R023 comparison.
This packaged run is not a new independent verification.
"""
import hashlib,json,resource,time
from fractions import Fraction as Q
import sympy as S
from models import ROOT,K,X,x,y,z,a,b,c,d,e,f,model
import rational_polynomials as R
start=time.monotonic()
t=S.Symbol('t',positive=True);u=S.Symbol('u',positive=True)
L=e*a*t+e*b-c*d;DA=t*L
QA=-e*a*a*t**3+a*c*d*t*t+b*(e*b+f*a-2*c*d)*t+b*b*f
Aprod={x:b*(b+a*t)/DA,y:b*d/DA,z:b*d/L}
AA=a+c*u;U=f*AA-c*d;DB=d*e*AA
QB=-b*d*e*AA+d*(a-c*u)*U+u*U**2
Bprod={x:u*U/(d*e),y:u,z:u*U/(e*AA)}
forms={33449:(t,QA,Aprod),239089:(u,QB,Bprod)}
identities={}
for index in [33449,239089]:
    m=model(index);variable,poly,point=forms[index]
    F=m['F'].subs(point,simultaneous=True).applyfunc(S.cancel)
    det=S.cancel(m['J'].det().subs(point,simultaneous=True))
    if index==33449:
        assert F[0]==F[1]==0 and S.cancel(F[2]-b*d*QA/DA**2)==0
        expected=-b*d*(S.diff(QA,t)*DA-2*QA*S.diff(DA,t))/DA**2
        tz=m['D']/(a*e*z)
        assert S.cancel(QA.subs(t,tz)+b*d*m['P']/(a*e**2*z**3))==0
        assert S.cancel(S.diff(tz,z)+b*d/(a*e*z*z))==0
        identities[str(index)]={'producer_residual_identity':True,'producer_determinant_identity':True,'coordinate':'t=z/y = D_A(z)/(a*e*z)','orientation':'strictly decreasing in z','polynomial_relation':'QA(t(z)) = -b*d*P_A(z)/(a*e**2*z**3)'}
    else:
        assert F[0]==F[2]==0 and S.cancel(F[1]-u*QB/DB)==0
        expected=u*U*((QB+u*S.diff(QB,u))*DB-u*QB*S.diff(DB,u))/DB**2
        yz=m['phi'][y]
        assert S.cancel(S.diff(yz,z)-a*d*(b*d+2*b*e*z+2*c*e*z*z)/m['D']**2)==0
        identities[str(index)]={'producer_residual_identity':True,'producer_determinant_identity':True,'coordinate':'u=y = a*z*(d+e*z)/D_B(z)','orientation':'strictly increasing in z','point_equivalence':'Both parametrizations solve the same independently generated ODE on their positive domains; the Cramer parametrization is unique for a fixed z.'}
    assert S.cancel(det-expected)==0
producer=json.loads((ROOT/'pilot-exact.json').read_text());ours=json.loads((ROOT/'independent-samples.json').read_text());sample_results=[]
assert producer['inputs_sha256']==ours['input_sha256']
for pr in producer['results']:
    own=next(v for v in ours['samples'] if v['id']==pr['sample'])
    index=own['candidate_index'];variable,poly,point=forms[index]
    rates=dict(zip(K,map(S.Rational,own['kappa'])))
    poly_num=S.Poly(poly.subs(rates),variable)
    qc=[Q(str(poly_num.nth(i))) for i in range(4)]
    assert qc==list(map(Q,pr['q_ascending']))
    sturm=R.sturm(qc)
    assert len(sturm[-1])==1
    assert R.count(sturm,Q(pr['admissible_lower']),Q(pr['strict_cauchy_upper']))==pr['distinct_positive_equilibria']==own['positive_equilibria']
    assert pr['squarefree'] and pr['stable']==own['exponentially_stable_equilibria']
    matches=[]
    for number,root in enumerate(pr['roots'],1):
        lo,hi=map(Q,root['interval']);assert lo>Q(pr['admissible_lower']) and R.count(sturm,lo,hi)==1
        image={}
        for axis in X:
            expr=point[axis].subs(rates)
            vals=[Q(str(S.cancel(expr.subs(variable,S.Rational(endpoint.numerator,endpoint.denominator))))) for endpoint in [lo,hi]]
            image[str(axis)]=(min(vals),max(vals))
        # Monotonicity on the positive domains makes these exact image bounds.
        for axis,reported in zip(['x','y','z'],root['coordinates']):
            reported=tuple(map(Q,reported));assert reported[0]<=image[axis][0]<=image[axis][1]<=reported[1]
        matching=[(i+1,r) for i,r in enumerate(own['roots']) if Q(r['z_interval'][0])<=image['z'][0] and image['z'][1]<=Q(r['z_interval'][1])]
        assert len(matching)==1;own_number,own_root=matching[0]
        for axis in ['x','y']:
            box=tuple(map(Q,own_root['value_intervals'][axis]));assert box[0]<=image[axis][0]<=image[axis][1]<=box[1]
        rh={}
        for key,bounds in zip(['b1','b2','b3','H2'],root['characteristic_b1_b2_b3_H2']):
            exact_bounds=tuple(map(Q,bounds));sgn=R.interval_sign(exact_bounds)
            assert sgn==own_root['signs'][key]
            independent_bounds=tuple(map(Q,own_root['value_intervals'][key]));assert max(exact_bounds[0],independent_bounds[0])<=min(exact_bounds[1],independent_bounds[1])
            rh[key]=sgn
        stable=root['verdict']=='stable_Routh_Hurwitz';assert stable==own_root['exponentially_stable']
        matches.append({'producer_root_order':number,'independent_z_order':own_number,'producer_Q_interval_certified':True,'mapped_coordinate_bounds_inside_reported_bounds':True,'mapped_coordinate_bounds_inside_independent_bounds':True,'RH_signs':rh,'exponentially_stable':stable,'status':'MATCH'})
    assert sorted(v['independent_z_order'] for v in matches)==[1,2,3]
    sample_results.append({'sample':own['id'],'positive_roots':own['positive_equilibria'],'stable':own['exponentially_stable_equilibria'],'matching':matches,'status':'MATCH'})
output={'independent_freeze_commit':'9c1dde75aa476ba15bc1bd97c847c2ad2658ad3d','run_kind':'package reconciliation; both original implementations are now available','symbolic_comparison':identities,'samples':sample_results,'mathematical_discrepancies':[],'presentation_differences':['A: independent variable z yields degree at most 3; producer t=z/y yields a cubic with nonzero leading coefficient for all positive rates. These are equivalent, not a contradiction.','A root order reverses between t and z; B order agrees between y and z.','Independent stability proof directly establishes RH positivity for these two networks; producer cites Tang-Zhang Lemma 3.'],'status':'MATCH','wall_seconds':time.monotonic()-start,'ru_maxrss_bytes_macos':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(ROOT/'comparison-check.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({'status':output['status'],'symbolic':identities,'samples':[{k:v for k,v in s.items() if k!='matching'} for s in sample_results]},indent=2))
