import json
from challenge import N,MASK,PM
W=1024
OUTW=2048
ALL=(1<<W)-1

class CNF:
    def __init__(self):
        self.nvars=0; self.clauses=[]; self.names={}
    def var(self,name=None):
        self.nvars+=1
        if name is not None: self.names[name]=self.nvars
        return self.nvars
    def add(self,*lits):
        out=[]; seen=set()
        for lit in lits:
            if lit is True: return
            if lit is False or lit==0: continue
            lit=int(lit)
            if -lit in seen: return
            if lit not in seen: seen.add(lit); out.append(lit)
        if not out: raise RuntimeError('empty clause generated')
        self.clauses.append(out)
    @staticmethod
    def neg(a):
        if a is True:return False
        if a is False:return True
        return -a
    def and_gate(self,a,b):
        if a is False or b is False:return False
        if a is True:return b
        if b is True:return a
        if a==b:return a
        if a==-b:return False
        z=self.var(); self.add(-a,-b,z); self.add(a,-z); self.add(b,-z); return z
    def or_gate(self,a,b):
        if a is True or b is True:return True
        if a is False:return b
        if b is False:return a
        if a==b:return a
        if a==-b:return True
        z=self.var(); self.add(a,b,-z); self.add(-a,z); self.add(-b,z); return z
    def xor_gate(self,a,b):
        if a is False:return b
        if b is False:return a
        if a is True:return self.neg(b)
        if b is True:return self.neg(a)
        if a==b:return False
        if a==-b:return True
        z=self.var(); self.add(-a,-b,-z); self.add(a,b,-z); self.add(a,-b,z); self.add(-a,b,z); return z
    def half_adder(self,a,b): return self.xor_gate(a,b),self.and_gate(a,b)
    def full_adder(self,a,b,c):
        vals=[a,b,c]; ones=sum(v is True for v in vals); xs=[v for v in vals if v is not True and v is not False]
        if not xs:return bool(ones&1),bool(ones>>1)
        if len(xs)==1:
            u=xs[0]
            if ones==0:return u,False
            if ones==1:return -u,u
            return u,True
        if len(xs)==2:
            u,v=xs
            if ones==0:return self.half_adder(u,v)
            return self.neg(self.xor_gate(u,v)),self.or_gate(u,v)
        a,b,c=xs; s=self.var(); d=self.var()
        self.add(a,b,-d); self.add(d,s,-c); self.add(c,-d,-s); self.add(d,-a,-b)
        self.add(a,b,c,-s); self.add(a,c,s,-b); self.add(b,c,s,-a)
        self.add(a,-b,-c,-s); self.add(b,-a,-c,-s); self.add(s,-a,-b,-c)
        return s,d
    def force(self,sig,value):
        if sig is True:
            if not value: raise RuntimeError('forcing true to false')
            return
        if sig is False:
            if value: raise RuntimeError('forcing false to true')
            return
        self.add(sig if value else -sig)
    def write(self,path):
        with open(path,'w',buffering=1024*1024) as f:
            f.write('p cnf %d %d\n'%(self.nvars,len(self.clauses)))
            for c in self.clauses:f.write(' '.join(map(str,c))+' 0\n')

def bit(v,i):return (v>>i)&1

def q_low_values():
    vals=[]; mod=1<<265
    for a in range(16):
        plow=(PM|(a<<150))&(mod-1); vals.append((N*pow(plow,-1,mod))%mod)
    return vals

def q_high_prefix_for_d(d):
    known=PM|(d<<920); fixed=MASK|(0xF<<920); pmin=known; pmax=known|(ALL^fixed)
    qmin=(N+pmax-1)//pmax; qmax=N//pmin; diff=qmin^qmax
    common=W-diff.bit_length() if diff else W
    return qmin,qmax,common,qmin>>(W-common)

def make_inputs(cnf,d):
    p=[None]*W
    for i in range(W):
        if (MASK>>i)&1:p[i]=bool((PM>>i)&1)
        elif 920<=i<924:p[i]=bool((d>>(i-920))&1)
        else:p[i]=cnf.var('p%d'%i)
    qvals=q_low_values(); q=[None]*W; alits=[p[150+k] for k in range(4)]
    assert all(isinstance(z,int) and not isinstance(z,bool) for z in alits)
    for j in range(265):
        vals=[bit(v,j) for v in qvals]
        if all(z==vals[0] for z in vals):q[j]=bool(vals[0])
        else:
            qj=cnf.var('q%d'%j); q[j]=qj
            for av in range(16):
                clause=[]
                for k,alit in enumerate(alits):clause.append(alit if bit(av,k)==0 else -alit)
                clause.append(qj if vals[av] else -qj); cnf.add(*clause)
    qmin,qmax,common,prefix=q_high_prefix_for_d(d)
    for j in range(W-common,W):q[j]=bool(bit(qmin,j))
    for j in range(W):
        if q[j] is None:q[j]=cnf.var('q%d'%j)
    return p,q,(qmin,qmax,common,prefix)

def compress_product(cnf,p,q):
    cols=[[] for _ in range(OUTW+12)]; const=0; ands=0; terms=0
    for i,pi in enumerate(p):
        if pi is False:continue
        for j,qj in enumerate(q):
            if qj is False:continue
            k=i+j
            if pi is True and qj is True:const+=1<<k; terms+=1
            elif pi is True:cols[k].append(qj);terms+=1
            elif qj is True:cols[k].append(pi);terms+=1
            else:cols[k].append(cnf.and_gate(pi,qj));ands+=1;terms+=1
    k=0
    while const:
        if const&1:cols[k].append(True)
        const>>=1;k+=1
    full=0
    for k in range(len(cols)-1):
        work=cols[k]
        while len(work)>=3:
            a=work.pop();b=work.pop();c=work.pop();s,carry=cnf.full_adder(a,b,c)
            if s is not False:work.append(s)
            if carry is not False:cols[k+1].append(carry)
            full+=1
    if len(cols[-1])>=3:raise RuntimeError('insufficient carry columns')
    carry=False
    for k,work in enumerate(cols):
        if len(work)>2:raise RuntimeError('column %d still has %d signals'%(k,len(work)))
        a=work[0] if len(work)>=1 else False; b=work[1] if len(work)>=2 else False
        s,carry=cnf.full_adder(a,b,carry); cnf.force(s,bit(N,k) if k<OUTW else 0)
    cnf.force(carry,False)
    return {'terms':terms,'ands':ands,'full_adders':full,'ripple_columns':len(cols)}

def build_instance(d,cnf_path,map_path):
    cnf=CNF();p,q,qrange=make_inputs(cnf,d);stats=compress_product(cnf,p,q);cnf.write(cnf_path)
    isvar=lambda z:isinstance(z,int) and not isinstance(z,bool)
    meta={'d':d,'nvars':cnf.nvars,'clauses':len(cnf.clauses),
          'p_ids':[int(z) if isvar(z) else (1 if z else 0) for z in p],
          'p_const':[None if isvar(z) else int(bool(z)) for z in p],
          'q_ids':[int(z) if isvar(z) else (1 if z else 0) for z in q],
          'q_const':[None if isvar(z) else int(bool(z)) for z in q],
          'qrange':list(qrange[:3]),'stats':stats}
    with open(map_path,'w') as f:json.dump(meta,f)
    return meta
