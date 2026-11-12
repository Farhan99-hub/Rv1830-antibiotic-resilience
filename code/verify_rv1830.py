
"""Independent re-verification of the Rv1830 (ResR) correspondence.

Reloads the three source files from scratch, rebuilds the analysis set,
recomputes every statistic quoted in the manuscript, then reads the .docx
back and checks each quoted value against the recomputed one.
"""
import re, sys, collections, json
import numpy as np, pandas as pd
from scipy import stats as ss
from statsmodels.stats.contingency_tables import Table2x2, StratifiedTable, Table
from docx import Document

DS1, DS2, DS3, DOCX = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

# ---------- load ----------
D1 = pd.ExcelFile(DS1).parse('Combined_17719_isolates')
D2 = pd.read_csv(DS2, low_memory=False)
D3 = pd.read_csv(DS3, low_memory=False)

# ---------- no shared isolates ----------
s1, s2, s3 = (set(D1['Isolate_ID'].astype(str)), set(D2['UNIQUEID'].astype(str)),
              set(D3['sra_id'].astype(str)))
assert not (s1 & s2) and not (s1 & s3) and not (s2 & s3), "collections overlap"

# ---------- variant classification ----------
def tokclass(t):
    t = t.strip()
    if not t or t in ('-', 'WT'): return None
    if re.match(r'^[acgt]-\d+', t) or t.startswith('-'): return 'upstream'
    m = re.fullmatch(r'([A-Z!])(\d+)([A-Z!*])', t)
    if m: return 'syn' if m.group(1) == m.group(3) else 'protein'
    if re.fullmatch(r'\d+_(ins|del)_[acgt]+', t): return 'indel_cds'
    return 'unknown:' + t

def prot_tokens(s):
    if str(s).strip() in ('-', 'nan'): return None          # no call
    return [t.strip() for t in re.split(r'[;,]', str(s))
            if tokclass(t) in ('protein', 'indel_cds')]

SEV = {0: 'Susceptible', 1: 'Other resistance', 2: 'MDR', 3: 'Pre-XDR', 4: 'XDR'}
m1 = {'Sensitive':0,'HR-TB':1,'RR-TB':1,'MDR-TB':2,'Pre-XDR-TB':3,'XDR-TB':4}
m2 = {'no_resistance_called':0,'resistant_other':1,'MDR':2,
      'MDR+FQ(pre-XDR)':3,'MDR+FQ+inj(XDR-like)':4}
m3 = {'Sensitive':0,'Mono DR':1,'Poly DR':1,'MDR non XDR':2,'Pre-XDR':3,'XDR':4}

p2 = [prot_tokens(s) for s in D2['Rv1830'].fillna('-')]
keep2 = np.array([t is not None for t in p2]) & (D2['resistance_class'] != 'untested').values
keep1 = (D1['Primary_endpoint'] != 'Excluded (Other)').values
keep3 = (D3['rv1830_status'] == 'ok').values

def subs1(s, m):
    return [t.strip() for t in re.split(r'[;,]', str(s)) if t.strip() and t.strip() != '-'] if m == 1 else []

A = pd.DataFrame({'coll':'C1','sev':D1['Phenotype'].map(m1),
                  'subs':[subs1(s,m) for s,m in zip(D1['Rv1830_substitutions'].fillna('-'), D1['Rv1830_mutant'])]})[keep1]
B = pd.DataFrame({'coll':'C2','sev':D2['resistance_class'].map(m2),
                  'subs':[t if t is not None else [] for t in p2]})[keep2]
C = pd.DataFrame({'coll':'C3','sev':D3['type_of_resistance'].map(m3),
                  'subs':[[t.strip() for t in re.split(r'[;,]', str(s)) if t.strip() and t.strip()!='-']
                          if pd.notna(s) and str(s)!='-' else [] for s in D3['rv1830_aa_changes']]})[keep3]
ALL = pd.concat([A,B,C], ignore_index=True)
ALL = ALL[ALL['sev'].notna()].copy()
ALL['sev'] = ALL['sev'].astype(int)
ALL['carrier'] = [len(L) > 0 for L in ALL['subs']]
ALL['resistant'] = ALL['sev'] > 0
ALL['D144N'] = [any(t.startswith('D144N') for t in L) for L in ALL['subs']]

R = {}
R['n_total'] = len(ALL)
R['n_by_coll'] = ALL['coll'].value_counts().to_dict()
R['n_carriers'] = int(ALL['carrier'].sum())
R['pct_carriers'] = round(100*ALL['carrier'].mean(), 2)
R['pct_by_coll'] = {k: round(100*v, 2) for k, v in ALL.groupby('coll')['carrier'].mean().items()}

def tab(flag, df):
    x = pd.crosstab(flag, df['resistant'])
    return np.array([[x.loc[True,True], x.loc[True,False]],
                     [x.loc[False,True], x.loc[False,False]]], float)

def orci(arr):
    o = Table2x2(arr)
    return (round(o.oddsratio,2), round(o.oddsratio_confint()[0],2),
            round(o.oddsratio_confint()[1],2), ss.fisher_exact(arr)[1])

A0 = tab(ALL['carrier'], ALL)
R['pct_res_carrier'] = round(100*A0[0,0]/A0[0].sum(), 1)
R['pct_res_noncarrier'] = round(100*A0[1,0]/A0[1].sum(), 1)
R['OR_crude'] = orci(A0)

strat = []
for k, g in ALL.groupby('coll'):
    arr = tab(g['carrier'], g); strat.append(arr); R[f'OR_{k}'] = orci(arr)
st = StratifiedTable(np.dstack(strat))
R['OR_MH'] = (round(st.oddsratio_pooled,2), round(st.oddsratio_pooled_confint()[0],2),
              round(st.oddsratio_pooled_confint()[1],2), st.test_null_odds().pvalue)
R['breslow_day_p'] = st.test_equal_odds().pvalue

gg = pd.crosstab(ALL['sev'], ALL['carrier'])
R['gradient'] = [round(100*gg.loc[i,True]/(gg.loc[i,True]+gg.loc[i,False]), 2) for i in sorted(gg.index)]
R['gradient_counts'] = [(int(gg.loc[i,True]), int(gg.loc[i,True]+gg.loc[i,False])) for i in sorted(gg.index)]
ordt = Table(gg.values).test_ordinal_association()
R['trend_z'], R['trend_p'] = round(ordt.zscore,2), ordt.pvalue

sD = ALL[ALL['D144N'] | ~ALL['carrier']]
R['n_D144N'] = int(ALL['D144N'].sum())
aD = tab(sD['D144N'], sD); R['pct_res_D144N'] = round(100*aD[0,0]/aD[0].sum(),1); R['OR_D144N'] = orci(aD)
other = ALL['carrier'] & ~ALL['D144N']
sO = ALL[other | ~ALL['carrier']]
R['n_other'] = int(other.sum())
aO = tab(other.loc[sO.index], sO); R['pct_res_other'] = round(100*aO[0,0]/aO[0].sum(),1); R['OR_other'] = orci(aO)

cnt = collections.Counter(t for L in ALL['subs'] for t in L)
R['n_distinct_subs'] = len(cnt)
def pos(t):
    m = re.match(r'^[A-Z!](\d+)', t); return int(m.group(1)) if m else None
tot = sum(n for t,n in cnt.items() if pos(t))
inr = sum(n for t,n in cnt.items() if pos(t) and 65 <= pos(t) <= 137)
R['changes_total'], R['changes_in_region'] = tot, inr
R['pct_in_region'] = round(100*inr/tot, 1)
R['distinct_in_region'] = (sum(1 for t in cnt if pos(t) and 65 <= pos(t) <= 137),
                           sum(1 for t in cnt if pos(t)))
R['top14'] = cnt.most_common(14)

# c-87t prevalence in the analysable MIC collection
s2v = D2.loc[keep2, 'Rv1830'].fillna('-')
R['c87t'] = (int(s2v.map(lambda v: 'c-87t' in [t.strip() for t in re.split(r'[;,]', str(v))]).sum()),
             int(keep2.sum()))

# MIC comparison
sub = D2[keep2].copy(); sub['carrier'] = [len(t) > 0 for t in np.array(p2, dtype=object)[keep2]]
def tomic(v):
    if pd.isna(v): return np.nan
    try: return float(str(v).replace('<=','').replace('>=','').replace('>','').replace('<',''))
    except ValueError: return np.nan
R['mic'] = {}
for drug in ['RIF','INH','EMB','LEV','MXF','AMI','KAN','ETH']:
    m = sub[f'{drug}_MIC'].map(tomic); ok = m.notna()
    x, y = np.log2(m[ok & sub['carrier']]), np.log2(m[ok & ~sub['carrier']])
    u = ss.mannwhitneyu(x, y, alternative='two-sided')
    R['mic'][drug] = dict(median_carrier=float(2**np.median(x)), median_noncarrier=float(2**np.median(y)),
                          diff_log2=round(float(x.mean()-y.mean()),3), p=round(float(u.pvalue),4),
                          sig_bonferroni=bool(u.pvalue < 0.05/8))

# ---------- check the manuscript against R ----------
doc = Document(DOCX)
txt = "\n".join(p.text for p in doc.paragraphs)
checks = [
    ("isolates analysed",        f"{R['n_total']:,}"),
    ("carriers",                 str(R['n_carriers'])),
    ("resistant in carriers",    f"{R['pct_res_carrier']} percent"),
    ("resistant in non-carriers", f"{R['pct_res_noncarrier']} percent"),
    ("crude OR",                 f"{R['OR_crude'][0]:.2f}-fold"),
    ("crude OR CI",              f"{R['OR_crude'][1]:.2f} to {R['OR_crude'][2]:.2f}"),
    ("Mantel-Haenszel OR",       f"{R['OR_MH'][0]:.2f} ({R['OR_MH'][1]:.2f} to {R['OR_MH'][2]:.2f}"),
    ("per-collection ORs",       f"{R['OR_C1'][0]:.2f}, {R['OR_C2'][0]:.2f} and {R['OR_C3'][0]:.2f}"),
    ("gradient susceptible",     f"{R['gradient'][0]:.2f} percent"),
    ("gradient other",           f"{R['gradient'][1]:.2f} percent"),
    ("gradient MDR",             f"{R['gradient'][2]:.2f} percent"),
    ("gradient pre-XDR",         f"{R['gradient'][3]:.2f} percent"),
    ("gradient XDR",             f"{R['gradient'][4]:.2f} percent"),
    ("trend z",                  f"z = {R['trend_z']:.2f}"),
    ("D144N carriers",           f"by {R['n_D144N']} isolates"),
    ("D144N resistant pct",      f"{R['pct_res_D144N']} percent of its carriers"),
    ("D144N OR",                 f"{R['OR_D144N'][0]:.2f} ({R['OR_D144N'][1]:.2f} to {R['OR_D144N'][2]:.2f}"),
    ("remaining carriers",       f"remaining {R['n_other']} carriers"),
    ("OR excluding D144N",       f"{R['OR_other'][0]:.2f} ({R['OR_other'][1]:.2f} to {R['OR_other'][2]:.2f}"),
    ("resistant pct excl D144N", f"{R['pct_res_other']} percent of them resistant"),
    ("distinct substitutions",   f"{R['n_distinct_subs']} distinct"),
    ("changes in region",        f"{R['changes_in_region']} of {R['changes_total']}"),
    ("pct in region",            f"{R['pct_in_region']} percent"),
    ("distinct in region",       f"{R['distinct_in_region'][0]} of the {R['distinct_in_region'][1]}"),
    ("c-87t count",              f"{R['c87t'][0]:,} of the {R['c87t'][1]:,}"),
    ("RIF MIC p",                f"p = {R['mic']['RIF']['p']:.2f}"),
    ("RIF medians",              f"{R['mic']['RIF']['median_carrier']:g} mg/L"),
    ("INH nominal p",            f"isoniazid (p = {R['mic']['INH']['p']:.3f})"),
    ("ETH nominal p",            f"ethionamide (p = {R['mic']['ETH']['p']:.3f})"),
]
fails = []
print(f"{'check':30s} {'expected string in manuscript':48s} ok")
for name, expect in checks:
    ok = expect in txt
    if not ok: fails.append((name, expect))
    print(f"{name:30s} {expect:48s} {'PASS' if ok else 'FAIL'}")

stale = [v for v in ['39,733','9,884','835','1.75','1.88','684','2.25','215'] if v in txt]
print("\nstale pre-correction values still present:", stale or "none")
print("no drug significant after Bonferroni:", not any(d['sig_bonferroni'] for d in R['mic'].values()))
json.dump({k: (list(v) if isinstance(v, tuple) else v) for k, v in R.items()},
          open("verification_results.json", "w"), indent=1, default=str)
print(f"\n{len(checks)-len(fails)}/{len(checks)} checks passed")
sys.exit(1 if (fails or stale) else 0)
