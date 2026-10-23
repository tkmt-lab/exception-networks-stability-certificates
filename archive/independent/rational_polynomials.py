"""Small independent Fraction polynomial/Sturm/interval routines, ascending order."""
from fractions import Fraction as Q
from functools import reduce

def trim(p):
    p=list(map(Q,p))
    while len(p)>1 and p[-1]==0:p.pop()
    return p

def zero(p):return trim(p)==[0]
def degree(p):return len(trim(p))-1
def evaluate(p,x):
    ans=Q(0)
    for c in reversed(p):ans=ans*x+c
    return ans

def derivative(p):return trim([i*p[i] for i in range(1,len(p))] or [0])
def divide(p,q):
    p=trim(p);q=trim(q);assert not zero(q)
    result=[Q(0)]*max(1,len(p)-len(q)+1)
    while not zero(p) and len(p)>=len(q):
        shift=len(p)-len(q);factor=p[-1]/q[-1];result[shift]+=factor
        for i,c in enumerate(q):p[i+shift]-=factor*c
        p=trim(p)
    return trim(result),p

def positive_normalize(p):
    p=trim(p)
    return p if zero(p) else [c/abs(p[-1]) for c in p]
def sturm(p):
    seq=[positive_normalize(p),positive_normalize(derivative(p))]
    assert not zero(seq[-1])
    while True:
        rem=divide(seq[-2],seq[-1])[1]
        if zero(rem):break
        seq.append(positive_normalize([-c for c in rem]))
    return seq

def sign(x):return (x>0)-(x<0)
def variations(signs):
    signs=[s for s in signs if s]
    return sum(a!=b for a,b in zip(signs,signs[1:]))
def variation_at(seq,x):return variations([sign(evaluate(p,x)) for p in seq])
def variation_infinity(seq):return variations([sign(p[-1]) for p in seq])
def count(seq,lo,hi):return variation_at(seq,lo)-variation_at(seq,hi)

def split(p,lo,hi):
    # Never use a root as a subdivision endpoint. More candidates than degree.
    for n in range(2,len(p)+3):
        m=((n-1)*lo+hi)/n
        if evaluate(p,m):return m
    raise AssertionError('Every subdivision candidate is a root')

def isolate_positive(p):
    p=trim(p);seq=sturm(p);assert evaluate(p,Q(0))!=0
    bound=1+max(abs(c/p[-1]) for c in p[:-1]);B=Q(bound.numerator//bound.denominator+1)
    assert evaluate(p,B)!=0
    total=variation_at(seq,Q(0))-variation_infinity(seq)
    assert count(seq,Q(0),B)==total
    stack=[(Q(0),B,total)];out=[]
    while stack:
        lo,hi,n=stack.pop()
        if n==0:continue
        if n==1 and hi-lo<Q(1,2**24):out.append((lo,hi));continue
        mid=split(p,lo,hi);left=count(seq,lo,mid);right=count(seq,mid,hi)
        assert left+right==n
        stack.extend([(mid,hi,right),(lo,mid,left)])
    return seq,B,sorted(out)

def refine(p,seq,interval):
    lo,hi=interval;mid=split(p,lo,hi)
    if count(seq,lo,mid)==1:return lo,mid
    assert count(seq,mid,hi)==1
    return mid,hi

def imul(a,b):
    products=[x*y for x in a for y in b];return min(products),max(products)
def iadd(a,b):return a[0]+b[0],a[1]+b[1]
def ipoly(p,interval):
    ans=(Q(0),Q(0))
    for c in reversed(p):ans=iadd(imul(ans,interval),(c,c))
    return ans

def irational(num,den,interval):
    n=ipoly(num,interval);d=ipoly(den,interval)
    if d[0]<=0<=d[1]:return None
    return imul(n,(1/d[1],1/d[0]))
def interval_sign(interval):
    if interval is None:return 0
    if interval[0]>0:return 1
    if interval[1]<0:return -1
    return 0

def strings(p):return [str(c) for c in p]
