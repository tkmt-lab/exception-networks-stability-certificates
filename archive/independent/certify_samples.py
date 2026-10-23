"""Exact fixed-sample certification; root decisions use only Fraction arithmetic."""
import hashlib,json,platform,resource,time
from fractions import Fraction as Q
import sympy as S
from models import ROOT,K,X,x,y,z,model
import rational_polynomials as R
start=time.monotonic()

def coefficients(expr):
    poly=S.Poly(expr,z,domain=S.QQ)
    return R.trim([Q(str(poly.nth(i))) for i in range(poly.degree()+1)] if not poly.is_zero else [0])

def rational_data(expr,P):
    num,den=S.cancel(expr).as_numer_denom()
    n=coefficients(num);d=coefficients(den)
    reduced=R.divide(n,P)[1]
    return reduced,d

output={'status':'PASS','sympy':S.__version__,'python':platform.python_version(),'method':'N,Y to mass-action ODE; Cramer parametrization; independent Fraction Sturm sequence and rational interval Horner evaluation; no floating-point root or sign decisions','samples':[]}
for sample in json.loads((ROOT/'input/samples.json').read_text())['samples']:
    m=model(sample['candidate_index']);rates={k:S.Rational(v) for k,v in zip(K,sample['kappa'])}
    P=coefficients(m['P'].subs(rates));D=coefficients(m['D'].subs(rates));phi={k:S.cancel(v.subs(rates)) for k,v in m['phi'].items()}
    seq,bound,intervals=R.isolate_positive(P)
    assert len(seq[-1])==1, 'Polynomial is not square-free; a separate degenerate analysis is required'
    J=m['J'].subs(rates).subs(phi)
    b1=S.cancel(-S.trace(J));b2=S.cancel(sum(J.extract(ids,ids).det() for ids in [(0,1),(0,2),(1,2)]));b3=S.cancel(-J.det());H2=S.cancel(b1*b2-b3)
    cp=J.charpoly().all_coeffs()
    assert cp[0]==1 and all(S.cancel(t-u)==0 for t,u in zip(cp[1:],[b1,b2,b3]))
    remainders=[]
    for expr in m['F'].subs(rates).subs(phi):
        num,den=S.cancel(expr).as_numer_denom();rem=R.divide(coefficients(num),P)[1];assert R.zero(rem);remainders.append(R.strings(rem))
    expressions={'x':phi[x],'y':phi[y],'b1':b1,'b2':b2,'b3':b3,'H2':H2,'detJ':-b3}
    data={key:rational_data(expr,P) for key,expr in expressions.items()}
    record={'id':sample['id'],'candidate_index':sample['candidate_index'],'kappa':sample['kappa'],'P_coefficients_ascending':R.strings(P),'D_coefficients_ascending':R.strings(D),'sturm_chain':list(map(R.strings,seq)),'positive_root_bound':str(bound),'sturm_at_zero':R.variation_at(seq,Q(0)),'sturm_at_positive_infinity':R.variation_infinity(seq),'distinct_positive_polynomial_roots':len(intervals),'square_free_gcd':'1','ode_numerator_remainders':remainders,'rational_expressions_mod_P':{key:{'numerator':R.strings(v[0]),'denominator':R.strings(v[1])} for key,v in data.items()},'roots':[]}
    for interval in intervals:
        for step in range(256):
            d_interval=R.ipoly(D,interval);dp_interval=R.ipoly(R.derivative(P),interval)
            dsign=R.interval_sign(d_interval)
            if dsign<0:break
            evaluated={key:R.irational(*value,interval) for key,value in data.items()}
            signs={key:R.interval_sign(v) for key,v in evaluated.items()}
            if dsign>0 and R.interval_sign(dp_interval) and all(signs.values()):break
            interval=R.refine(P,seq,interval)
        else:raise AssertionError('Rational sign isolation did not finish')
        lo,hi=interval;assert R.count(seq,lo,hi)==1
        root={'z_interval':R.strings(interval),'sturm_variations':[R.variation_at(seq,lo),R.variation_at(seq,hi)],'admissible':dsign>0,'D_interval':R.strings(d_interval),'Pprime_interval':R.strings(dp_interval),'Pprime_sign':R.interval_sign(dp_interval)}
        if dsign>0:
            assert signs['x']==signs['y']==1
            assert signs['detJ']==-R.interval_sign(dp_interval)
            assert signs['b1']==signs['b2']==signs['H2']==1
            first_column=[1,signs['b1'],signs['H2']*signs['b1'],signs['b3']]
            stable=all(signs[key]>0 for key in ['b1','b2','b3','H2'])
            root.update({'value_intervals':{key:R.strings(v) for key,v in evaluated.items()},'signs':signs,'routh_first_column_signs':first_column,'right_half_plane_eigenvalues':R.variations(first_column),'exponentially_stable':stable,'nondegenerate':signs['detJ']!=0})
        record['roots'].append(root)
    positive=[r for r in record['roots'] if r['admissible']]
    record['positive_equilibria']=len(positive);record['exponentially_stable_equilibria']=sum(r['exponentially_stable'] for r in positive)
    output['samples'].append(record)
    print(sample['id'], 'positive equilibria',len(positive),'stable',record['exponentially_stable_equilibria'])
    for root in positive:print('  z:',root['z_interval'],'det:',root['signs']['detJ'],'Routh:',root['routh_first_column_signs'])
output['wall_seconds']=time.monotonic()-start;output['ru_maxrss_bytes_macos']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
output['input_sha256']={name:hashlib.sha256((ROOT/'input'/name).read_bytes()).hexdigest() for name in ['networks.json','samples.json']}
(ROOT/'independent-samples.json').write_text(json.dumps(output,indent=2)+'\n')
