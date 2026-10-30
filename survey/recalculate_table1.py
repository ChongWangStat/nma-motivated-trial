"""Recompute the manuscript's sampled survey percentages; does not verify LLM labels."""
import csv,math
from pathlib import Path
p=Path(__file__).with_name('survey_codes_numeric.csv')
records=list(csv.DictReader(p.open()))
assert len(records)==256
for stratum in ('NMA','MA_only'):
 R=[r for r in records if r['stratum']==stratum]
 P=[r for r in R if r['report_type']=='primary']; n=len(P)
 print(stratum,'sampled',len(R),'primary',n)
 for label,fn in [
   ('related meta-analysis',lambda x:x['related']=='yes'),
   ('related NMA',lambda x:x['related_nma']=='yes'),
   ('motivating',lambda x:x['motivation']=='2'),
   ('promising',lambda x:x['motivation']=='2' and x['promising_signal']=='yes'),
   ('design',lambda x:x['design_use']=='yes'),
   ('discussion',lambda x:x['contextualised']=='yes'),
   ('combined',lambda x:x['combined']=='yes')]:
  k=sum(fn(x) for x in P);print('  %-25s %3d/%d = %5.1f%%'%(label,k,n,100*k/n))
