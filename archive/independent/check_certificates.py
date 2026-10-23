"""Second arithmetic route: IVT+cubic degree and SymPy exact root counts for signs."""
import json,time,resource
import sympy as S
from models import ROOT,K,x,y,z,model
start=time.monotonic();data=json.loads((ROOT/'independent-samples.json').read_text());results=[]
for rec in data['samples']:
    m=model(rec['candidate_index']);rates=dict(zip(K,map(S.Rational,rec['kappa'])))
    P=S.Poly(m['P'].subs(rates),z);D=S.Poly(m['D'].subs(rates),z)
    roots=rec['roots'];assert len(roots)==3 and P.degree()==3
    assert D.nth(2)<=0  # concave or affine, hence interval minimum at an endpoint
    previous=S.Rational(0);sign_checks=0
    for root in roots:
        lo,hi=map(S.Rational,root['z_interval']);assert previous<lo<hi;previous=hi
        assert P.eval(lo)*P.eval(hi)<0
        assert D.eval(lo)>0 and D.eval(hi)>0
        assert P.diff().count_roots(lo,hi)==0
        assert S.sign(P.diff().eval((lo+hi)/2))==root['Pprime_sign']
        # Zero-free numerator and denominator imply a constant exact sign.
        for key,value in rec['rational_expressions_mod_P'].items():
            numerator=S.Poly.from_list(list(reversed(list(map(S.Rational,value['numerator'])))),z)
            denominator=S.Poly.from_list(list(reversed(list(map(S.Rational,value['denominator'])))),z)
            assert numerator.count_roots(lo,hi)==0 and denominator.count_roots(lo,hi)==0
            exact_sign=int(S.sign(numerator.eval((lo+hi)/2)/denominator.eval((lo+hi)/2)))
            assert exact_sign==root['signs'][key];sign_checks+=1
        assert root['nondegenerate'] and root['admissible']
    # Each interval has a real root by IVT; three disjoint intervals exhaust a cubic.
    # Consequently the roots are all simple. This count does not use our Sturm code.
    assert S.gcd(P,P.diff()).degree()==0
    assert P.count_roots(0,S.oo)==3
    results.append({'id':rec['id'],'ivt_disjoint_intervals':3,'polynomial_degree':3,'gcd_degree':0,'separate_sympy_positive_root_count':3,'rational_sign_checks':sign_checks,'status':'PASS'})
result={'method':'Independent count by three sign-changing intervals + cubic degree; separate SymPy Sturm root counts for numerator/denominator signs; no imports from rational_polynomials.py','samples':results,'wall_seconds':time.monotonic()-start,'ru_maxrss_bytes_macos':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'status':'PASS'}
(ROOT/'certificate-recheck.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
