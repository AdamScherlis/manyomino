"""Cells in each ring of the made pictures against the average real ring over many roots.

analysis/e2e_test.py pairs each made picture with the real cells round the same root.  For the count of cells in
a ring that pairing buys little, since the real ring is only loosely tied to the 32 cells in the middle that the
two share, and the real rings of a few hundred roots are themselves a per cent or two off the average.  This sets
the made rings (and the real ones of the same roots) against 100 roots from each sample.

  python3 analysis/ring_pop.py [TAG=final ...]       # reads zoomtest/e2e/TAG_pics.npz; writes analysis/ring_pop.json
"""
import numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e2e_test import read, ROOT
edges=[16,32,64,128,256,512]
def rings_from_pics(tag,which):
    d=np.load(ROOT+'/zoomtest/e2e/%s_pics.npz'%tag); p=np.unpackbits(d[which],axis=-1).astype(bool)
    n,G=p.shape[0],p.shape[1]; R=G//2; yy,xx=np.mgrid[0:G,0:G]; a=np.where(xx<R,R-1-xx,xx-R); b=np.where(yy<R,R-1-yy,yy-R); rad=np.maximum(a,b)
    c=np.array([np.bincount(rad[q],minlength=R)[:R] for q in p],float); return np.array([c[:,lo:hi].sum(1) for lo,hi in zip(edges[:-1],edges[1:]) if hi<=R]).T
rng=np.random.default_rng(5); pop=[]
snaps=sorted(f for f in os.listdir(ROOT+'/prod') if f.startswith('snap_') and 'n1024000' in f and f.endswith('.bin'))
for s in snaps:
    T=read(ROOT+'/prod/'+s)
    for _ in range(100):
        c=T[rng.integers(len(T))]; d=T-c; a=np.where(d[:,0]<0,-d[:,0]-1,d[:,0]); b=np.where(d[:,1]<0,-d[:,1]-1,d[:,1]); r=np.maximum(a,b); r=r[r<512]
        h=np.bincount(r,minlength=512); pop.append([h[lo:hi].sum() for lo,hi in zip(edges[:-1],edges[1:])])
pop=np.array(pop,float); pm=pop.mean(0); ps=pop.std(0,ddof=1)/np.sqrt(len(pop))
print('real, %d roots of %d samples: rings %s (s.e. %s %%)'%(len(pop),len(snaps),' '.join('%.0f'%v for v in pm),' '.join('%.1f'%(100*s/m) for s,m in zip(ps,pm))))
out={'real roots':len(pop),'rings':[float(v) for v in pm],'s.e. %':[float(100*a/b) for a,b in zip(ps,pm)],'runs':{}}
for tag in (sys.argv[1:] or ['final']):
    for which in ['real','made']:
        r=rings_from_pics(tag,which); k=r.shape[1]; m=r.mean(0); se=r.std(0,ddof=1)/np.sqrt(len(r))
        out['runs'][tag+' '+which]={'roots':len(r),'against all real roots, %':[float(100*(m[i]/pm[i]-1)) for i in range(k)],'in s.e.':[float((m[i]-pm[i])/np.hypot(se[i],ps[i])) for i in range(k)]}
        print('%-8s %-4s %d roots: against all real roots %s'%(tag,which,len(r),'  '.join('%+.1f%% (%.1f s.e.)'%(100*(m[i]/pm[i]-1),(m[i]-pm[i])/np.hypot(se[i],ps[i])) for i in range(k))))
import json; json.dump(out,open(ROOT+'/analysis/ring_pop.json','w'),indent=1)
