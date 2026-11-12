"""Robustness analyses for the Rv1830 (ResR) correspondence.

Rebuilds the 40,276-isolate analysis set from the three source files (same rule
as verify_rv1830.py) and runs four robustness analyses:
  1. Confounding by phylogenetic population structure (lineage adjustment).
  2. Between-collection heterogeneity (random-effects meta-analysis).
  3. MIC comparison stratified by target-gene status.
  4. Treatment-outcome analysis (Collection 3, adjusted).
Writes robustness_results.json and prints a summary.

Usage:
  python robustness_rv1830.py "Data set 1.xlsx" "Data set 2.csv" "Data set 3.csv"
Requires: pandas, numpy, scipy, statsmodels, openpyxl.
"""
import re, sys, math, json, collections, warnings
import numpy as np, pandas as pd
from scipy import stats as ss
import statsmodels.formula.api as smf
warnings.filterwarnings("ignore")

DS1, DS2, DS3 = sys.argv[1], sys.argv[2], sys.argv[3]
D1 = pd.ExcelFile(DS1).parse('Combined_17719_isolates')
D2 = pd.read_csv(DS2, low_memory=False)
D3 = pd.read_csv(DS3, low_memory=False)

def tokclass(t):
    t = t.strip()
    if not t or t in ('-', 'WT'): return None
    if re.match(r'^[acgt]-\d+', t) or t.startswith('-'): return 'upstream'
    m = re.fullmatch(r'([A-Z!])(\d+)([A-Z!*])', t)
    if m: return 'syn' if m.group(1) == m.group(3) else 'protein'
    if re.fullmatch(r'\d+_(ins|del)_[acgt]+', t): return 'indel_cds'
    return 'unknown'
def prot(s):
    if str(s).strip() in ('-', 'nan'): return None
    return [t.strip() for t in re.split(r'[;,]', str(s)) if tokclass(t) in ('protein', 'indel_cds')]
def s1(s, m):
    return [t.strip() for t in re.split(r'[;,]', str(s)) if t.strip() and t.strip() != '-'] if m == 1 else []
def majL(x):
    x = str(x)
    if ';' in x or x.lower() == 'mixed': return 'Mixed'
    x = x.replace('lineage', 'L').replace('Lineage', 'L')
    if x.startswith('L') and len(x) > 1 and x[1].isdigit(): return 'L' + x[1]
    if 'bovis' in x.lower() or x.startswith('La1'): return 'M.bovis'
    if 'orygis' in x.lower(): return 'M.orygis'
    return 'Other'

m1 = {'Sensitive':0,'HR-TB':1,'RR-TB':1,'MDR-TB':2,'Pre-XDR-TB':3,'XDR-TB':4}
m2 = {'no_resistance_called':0,'resistant_other':1,'MDR':2,'MDR+FQ(pre-XDR)':3,'MDR+FQ+inj(XDR-like)':4}
m3 = {'Sensitive':0,'Mono DR':1,'Poly DR':1,'MDR non XDR':2,'Pre-XDR':3,'XDR':4}

p2 = [prot(s) for s in D2['Rv1830'].fillna('-')]
k2 = np.array([t is not None for t in p2]) & (D2['resistance_class'] != 'untested').values
k1 = (D1['Primary_endpoint'] != 'Excluded (Other)').values
k3 = (D3['rv1830_status'] == 'ok').values

A = pd.DataFrame({'coll':'C1','sev':D1['Phenotype'].map(m1),'lin':D1['Major_lineage'].map(majL),
     'subs':[s1(s,m) for s,m in zip(D1['Rv1830_substitutions'].fillna('-'), D1['Rv1830_mutant'])]})[k1]
B = pd.DataFrame({'coll':'C2','sev':D2['resistance_class'].map(m2),'lin':D2['LINEAGE'].map(majL),
     'subs':[t if t is not None else [] for t in p2]})[k2]
Cc = pd.DataFrame({'coll':'C3','sev':D3['type_of_resistance'].map(m3),'lin':D3['main_lineage'].map(majL),
     'outcome':D3['outcome'].astype(str),
     'subs':[[t.strip() for t in re.split(r'[;,]', str(s)) if t.strip() and t.strip() != '-'] if pd.notna(s) and str(s) != '-' else [] for s in D3['rv1830_aa_changes']]})[k3]
ALL = pd.concat([A, B, Cc], ignore_index=True)
ALL = ALL[ALL['sev'].notna()].copy(); ALL['sev'] = ALL['sev'].astype(int)
ALL['carrier'] = np.array([len(L) > 0 for L in ALL['subs']]).astype(int)
ALL['resistant'] = (ALL['sev'] > 0).astype(int)
ALL['D144N'] = np.array([any(t.startswith('D144N') for t in L) for L in ALL['subs']]).astype(int)
ALL['nonD'] = ((ALL['carrier'] == 1) & (ALL['D144N'] == 0)).astype(int)
vc = ALL['lin'].value_counts(); ALL['linC'] = ALL['lin'].where(ALL['lin'].isin(vc[vc >= 50].index), 'Other')

def orrow(res, term):
    b = res.params[term]; ci = res.conf_int().loc[term]
    return round(math.exp(b), 2), round(math.exp(ci[0]), 2), round(math.exp(ci[1]), 2)
def OR2(a, b, c, d):
    o = (a*d)/(b*c); se = math.sqrt(1/a+1/b+1/c+1/d)
    return round(o, 2), round(o*math.exp(-1.96*se), 2), round(o*math.exp(1.96*se), 2)

R = {'n_total': int(len(ALL)), 'n_carriers': int(ALL['carrier'].sum())}

# ---- 1. population structure ----
crude = smf.logit("resistant ~ carrier", data=ALL).fit(disp=0)
adj   = smf.logit("resistant ~ carrier + C(linC) + C(coll)", data=ALL).fit(disp=0)
adj2  = smf.logit("resistant ~ nonD + D144N + C(linC) + C(coll)", data=ALL).fit(disp=0)
l4 = ALL[ALL['lin'] == 'L4']
R['lineage'] = {
    'carrier_OR_crude': orrow(crude, 'carrier'),
    'carrier_OR_adj_lineage_collection': orrow(adj, 'carrier'),
    'nonD144N_OR_adj': orrow(adj2, 'nonD'),
    'D144N_OR_adj': orrow(adj2, 'D144N'),
    'nonD144N_within_L4_OR': OR2(int(((l4['nonD']==1)&(l4['resistant']==1)).sum()),
                                 int(((l4['nonD']==1)&(l4['resistant']==0)).sum()),
                                 int(((l4['carrier']==0)&(l4['resistant']==1)).sum()),
                                 int(((l4['carrier']==0)&(l4['resistant']==0)).sum())),
    'D144N_lineage': ALL[ALL['D144N']==1]['lin'].value_counts().to_dict(),
}

# ---- 2. heterogeneity (DerSimonian-Laird) ----
logors, ses = [], []
for k, g in ALL.groupby('coll'):
    a=int(((g['carrier']==1)&(g['resistant']==1)).sum()); b=int(((g['carrier']==1)&(g['resistant']==0)).sum())
    c=int(((g['carrier']==0)&(g['resistant']==1)).sum()); d=int(((g['carrier']==0)&(g['resistant']==0)).sum())
    logors.append(math.log((a*d)/(b*c))); ses.append(math.sqrt(1/a+1/b+1/c+1/d))
logors=np.array(logors); ses=np.array(ses); w=1/ses**2
fe=np.sum(w*logors)/np.sum(w); Q=float(np.sum(w*(logors-fe)**2)); dfQ=len(logors)-1
C_=np.sum(w)-np.sum(w**2)/np.sum(w); tau2=max(0,(Q-dfQ)/C_); I2=max(0,100*(Q-dfQ)/Q)
wr=1/(ses**2+tau2); re=np.sum(wr*logors)/np.sum(wr); sere=math.sqrt(1/np.sum(wr))
R['heterogeneity'] = {'random_effects_OR': [round(math.exp(re),2), round(math.exp(re-1.96*sere),2), round(math.exp(re+1.96*sere),2)],
                      'I2_percent': round(I2,0), 'tau2': round(tau2,3), 'Q_p': round(float(ss.chi2.sf(Q,dfQ)),4)}

# ---- 3. MIC stratified ----
def tomic(v):
    if pd.isna(v): return np.nan
    try: return float(str(v).replace('<=','').replace('>=','').replace('>','').replace('<',''))
    except ValueError: return np.nan
sub = D2[k2].copy(); sub['carrier'] = np.array([len(t)>0 for t in np.array(p2,dtype=object)[k2]])
def mic_test(mask, drug):
    m = sub.loc[mask, drug+'_MIC'].map(tomic); ok = m.notna()
    cm = mask & ok & sub['carrier']; nm = mask & ok & (~sub['carrier'])
    x = np.log2(sub.loc[cm, drug+'_MIC'].map(tomic)); y = np.log2(sub.loc[nm, drug+'_MIC'].map(tomic))
    if len(x) < 5 or len(y) < 5: return None
    return {'nC':int(len(x)),'nN':int(len(y)),'medC':round(float(2**np.median(x)),3),
            'medN':round(float(2**np.median(y)),3),'p':round(float(ss.mannwhitneyu(x,y,alternative='two-sided').pvalue),4)}
inhwt=(sub['katG'].astype(str)=='WT')&(sub['inhA'].astype(str)=='WT')&(sub['fabG1'].astype(str)=='WT')
R['mic_stratified'] = {
    'RIF_rpoB_WT': mic_test(sub['rpoB'].astype(str)=='WT','RIF'),
    'RIF_rpoB_mut': mic_test(sub['rpoB'].astype(str)!='WT','RIF'),
    'INH_target_WT': mic_test(inhwt,'INH'),
    'INH_target_mut': mic_test(~inhwt,'INH'),
}

# ---- 4. outcomes (C3) ----
c3 = ALL[ALL['coll']=='C3'].copy()
def outcome(unfav_set, name):
    keep = c3[c3['outcome'].isin({'Cured','Completed'} | unfav_set)].copy()
    keep['unfav'] = keep['outcome'].isin(unfav_set).astype(int)
    keep['linC'] = keep['lin'].where(keep['lin'].isin(['L2','L4']),'Other')
    a=int(((keep['carrier']==1)&(keep['unfav']==1)).sum()); b=int(((keep['carrier']==1)&(keep['unfav']==0)).sum())
    cc=int(((keep['carrier']==0)&(keep['unfav']==1)).sum()); d=int(((keep['carrier']==0)&(keep['unfav']==0)).sum())
    mod = smf.logit("unfav ~ carrier + C(sev) + C(linC)", data=keep).fit(disp=0)
    return {'definition':name,'n':int(len(keep)),'carrier_unfav_pct':round(100*a/(a+b),1),
            'noncarrier_unfav_pct':round(100*cc/(cc+d),1),'crude_OR':OR2(a,b,cc,d),
            'adj_OR':orrow(mod,'carrier'),'adj_p':round(float(mod.pvalues['carrier']),3)}
R['outcomes'] = {
    'broad_died_fail_ltfu': outcome({'Died','Failure','Lost to follow up'},'Died+Failure+LTFU'),
    'strict_died_fail': outcome({'Died','Failure'},'Died+Failure'),
}

json.dump(R, open("robustness_results.json","w"), indent=1, default=str)
print(json.dumps(R, indent=1, default=str))
print("\nrobustness_results.json written")
