#!/usr/bin/env python3
import argparse,json,os,subprocess,sys,time
from challenge import N,E,CT,MASK,PM
from satlib import build_instance

def parse_model(text):
    sat=None; values={}
    for line in text.splitlines():
        if line.startswith('s '):
            if 'SATISFIABLE' in line and 'UNSATISFIABLE' not in line:sat=True
            if 'UNSATISFIABLE' in line:sat=False
        if line.startswith('v '):
            for z in line.split()[1:]:
                n=int(z)
                if n:values[abs(n)]=n>0
    return sat,values

def recover(meta,values):
    def get(ids,consts):
        out=0
        for i,(vid,c) in enumerate(zip(ids,consts)):
            out|=(c if c is not None else int(values.get(vid,False)))<<i
        return out
    return get(meta['p_ids'],meta['p_const']),get(meta['q_ids'],meta['q_const'])

def decrypt(p,q):
    phi=(p-1)*(q-1);priv=pow(E,-1,phi);m=pow(CT,priv,N)
    return m.to_bytes((m.bit_length()+7)//8,'big')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--d',type=int,required=True);ap.add_argument('--solver',required=True)
    ap.add_argument('--work',default='sat-work');ap.add_argument('--timeout',type=int,default=19000);args=ap.parse_args()
    assert 0<=args.d<16;os.makedirs(args.work,exist_ok=True)
    cnf=os.path.join(args.work,'d%02d.cnf'%args.d);mapping=os.path.join(args.work,'d%02d.json'%args.d)
    t=time.time();meta=build_instance(args.d,cnf,mapping)
    print('CNF_META',json.dumps(meta['stats'],sort_keys=True),'vars',meta['nvars'],'clauses',meta['clauses'],
          'q_common_high',meta['qrange'][2],'gen_sec',round(time.time()-t,2),flush=True)
    t=time.time();cp=subprocess.run(['timeout','%ds'%args.timeout,args.solver,cnf],stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT,text=True,check=False)
    print(cp.stdout,flush=True);print('SOLVER_EXIT',cp.returncode,'elapsed',round(time.time()-t,2),flush=True)
    sat,values=parse_model(cp.stdout)
    if sat is not True:
        print('NO_SAT_MODEL d=%d sat=%r'%(args.d,sat),flush=True);return 1
    p,q=recover(meta,values);print('MODEL_P=0x%x'%p,flush=True);print('MODEL_Q=0x%x'%q,flush=True)
    if p*q!=N or (p&MASK)!=PM:
        print('MODEL_FAILED_VERIFY product=%s mask=%s'%(p*q==N,(p&MASK)==PM),flush=True);return 2
    raw=decrypt(p,q)
    print('FOUND_P=0x%x'%p,flush=True);print('FOUND_Q=0x%x'%q,flush=True)
    print('FOUND_P_DEC=%d'%p,flush=True);print('FOUND_Q_DEC=%d'%q,flush=True)
    print('PLAINTEXT_HEX=%s'%raw.hex(),flush=True);print('PLAINTEXT_REPR=%r'%raw,flush=True)
    with open(os.path.join(args.work,'FOUND.txt'),'w') as f:
        f.write('p=0x%x\nq=0x%x\nplaintext_hex=%s\nplaintext_repr=%r\n'%(p,q,raw.hex(),raw))
    return 0

if __name__=='__main__':sys.exit(main())
