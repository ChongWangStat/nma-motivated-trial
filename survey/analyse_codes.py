## Run from the survey/ folder: python3 analyse_codes.py
import json, math, csv, glob
A={r['pmcid']:r for f in sorted(glob.glob('codes/A_batch_*.json')) for r in json.load(open(f))}
B={r['pmcid']:r for f in sorted(glob.glob('codes/B_rel_*.json')) for r in json.load(open(f))}
def kappa(a,b):
    n=len(a); po=sum(x==y for x,y in zip(a,b))/n; cats=set(a)|set(b)
    pe=sum((a.count(c)/n)*(b.count(c)/n) for c in cats); return po,(po-pe)/(1-pe) if pe<1 else float('nan')
ids2=[p for p in B if p in A]
for f in ['report_type','related','related_nma','motivation','promising_signal','design_use','contextualised','combined']:
    po,k=kappa([str(A[p][f]).lower() for p in ids2],[str(B[p][f]).lower() for p in ids2]); print(f'agreement {f:18s} {po:.2f} kappa {k:.2f}')
## adjudication of the 10 reports on which the two codings disagreed on a key field (lower level when ambiguous)
adj={'PMC11569494':{'motivation':'1','promising_signal':'no','motivating_ma_type':'none'},
     'PMC11644361':{'related':'no','related_nma':'no'},
     'PMC11465541':{'motivation':'1','promising_signal':'no','motivating_ma_type':'none'},
     'PMC11348930':{'motivation':'1','promising_signal':'no','motivating_ma_type':'none'},
     'PMC11660996':{'related':'no','related_nma':'no'},
     'PMC11543726':{'report_type':'not_rct'},
     'PMC11287987':{'related':'no'},
     'PMC11256146':{'related':'no','related_nma':'no'},
     'PMC11492218':{'related':'no','motivation':'1','promising_signal':'no','motivating_ma_type':'none'},
     'PMC10785538':{'motivation':'1','promising_signal':'no','motivating_ma_type':'none'}}
F={}
for p,r in A.items():
    d={k:str(v).strip().lower() for k,v in r.items() if k not in ('evidence','note','pmcid')}
    if p in B:
        b={k:str(v).strip().lower() for k,v in B[p].items() if k not in ('evidence','note','pmcid')}
        for k in d:
            if d[k]!=b[k] and k not in adj.get(p,{}):
                # unresolved minor field disagreement: take the lower/negative level
                order=['0','1','2']; 
                if k=='motivation': d[k]=min(d[k],b[k])
                elif d[k] in ('yes','no') and b[k] in ('yes','no'): d[k]='no'
        d.update(adj.get(p,{}))
    d['pmcid']=p; d['evidence']=A[p]['evidence']; d['double_coded']= p in B
    F[p]=d
S=json.load(open('sample_ids.json')); nma=S['nma_open_access_coded']; ma=S['ma_only_open_access_coded']
for p in nma: F[p]['stratum']='NMA'
for p in ma: F[p]['stratum']='MA_only'
with open('survey_codes_final.csv','w',newline='') as f:
    keys=['pmcid','stratum','double_coded','report_type','related','related_nma','motivation','promising_signal','design_use','contextualised','combined','motivating_ma_type','evidence']
    w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); [w.writerow({k:F[p].get(k) for k in keys}) for p in nma+ma]
def wilson(x,n,z=1.96):
    if n==0: return (float('nan'),)*3
    p=x/n; den=1+z*z/n; c=(p+z*z/(2*n))/den; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return p, max(0,c-h), min(1,c+h)
def pr(lbl,x,n): p,l,u=wilson(x,n); print(f'  {lbl:55s} {x:4d}/{n:<4d} {100*p:5.1f}% ({100*l:.1f}-{100*u:.1f})'); return p
res={}
for s,ids in [('NMA',nma),('MA_only',ma)]:
    R=[F[p] for p in ids]; n=len(R)
    print('==',s,'n coded',n, {t:sum(r['report_type']==t for r in R) for t in ['primary','secondary_analysis','protocol','not_rct']})
    P=[r for r in R if r['report_type']=='primary']; k=len(P)
    out={}
    out['related']=pr('related (MA same condition & intervention)',sum(r['related']=='yes' for r in P),k)
    out['related_nma']=pr('related NMA',sum(r['related_nma']=='yes' for r in P),k)
    out['mot2']=pr('motivation=2 (MA given as reason for trial)',sum(r['motivation']=='2' for r in P),k)
    out['mot12']=pr('motivation>=1 (MA cited in intro/rationale)',sum(r['motivation'] in('1','2') for r in P),k)
    m2=[r for r in P if r['motivation']=='2']
    pr('  promising signal among motivated',sum(r['promising_signal']=='yes' for r in m2),len(m2))
    out['prom']=pr('promising signal (of all primary)',sum(r['promising_signal']=='yes' and r['motivation']=='2' for r in P),k)
    pr('  motivating MA type = nma/both among motivated',sum(r['motivating_ma_type'] in ('nma','both') for r in m2),len(m2))
    out['mot_nma']=pr('motivated by an NMA (of all primary)',sum(r['motivation']=='2' and r['motivating_ma_type'] in ('nma','both') for r in P),k)
    out['design']=pr('design use',sum(r['design_use']=='yes' for r in P),k)
    out['context']=pr('contextualised in discussion',sum(r['contextualised']=='yes' for r in P),k)
    out['comb']=pr('combined with MA',sum(r['combined']=='yes' for r in P),k)
    rel=[r for r in P if r['related']=='yes']
    pr('  motivation=2 among related',sum(r['motivation']=='2' for r in rel),len(rel))
    pr('  contextualised among related',sum(r['contextualised']=='yes' for r in rel),len(rel))
    res[s]=(out,k)
# weighted overall over trials citing any MA (OA frame): NMA 655, MA-only 5396
W={'NMA':655,'MA_only':5396}
print('== weighted (trials citing any MA, n frame 6051)')
for key in ['related','mot2','prom','design','context','comb']:
    v=sum(W[s]*res[s][0][key] for s in W)/sum(W.values()); print(f'  {key:10s} {100*v:.1f}%')
