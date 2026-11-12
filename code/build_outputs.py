"""Rebuild the analysis set, Table 1 and the substitution list from the raw files.

Reconstructs the 40,276-isolate analysis set (same rule as verify_rv1830.py) and
writes:
  analysis_set_40276.csv                      one row per isolate
  Table1_main.xlsx                            carriage by resistance category
  Supplementary_Table_S2_substitutions.xlsx   all 202 distinct substitutions

Usage:
  python build_outputs.py "Data set 1.xlsx" "Data set 2.csv" "Data set 3.csv" [out_dir]
Requires: pandas, numpy, openpyxl.
"""
import re, sys, math, collections
import numpy as np, pandas as pd, openpyxl
from openpyxl.styles import Font, PatternFill, Alignment

DS1, DS2, DS3 = sys.argv[1], sys.argv[2], sys.argv[3]
OUT = sys.argv[4] if len(sys.argv) > 4 else "."
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

SEV = {0:'Susceptible',1:'Other resistance',2:'MDR',3:'Pre-XDR',4:'XDR'}
m1 = {'Sensitive':0,'HR-TB':1,'RR-TB':1,'MDR-TB':2,'Pre-XDR-TB':3,'XDR-TB':4}
m2 = {'no_resistance_called':0,'resistant_other':1,'MDR':2,'MDR+FQ(pre-XDR)':3,'MDR+FQ+inj(XDR-like)':4}
m3 = {'Sensitive':0,'Mono DR':1,'Poly DR':1,'MDR non XDR':2,'Pre-XDR':3,'XDR':4}

p2 = [prot(s) for s in D2['Rv1830'].fillna('-')]
k2 = np.array([t is not None for t in p2]) & (D2['resistance_class'] != 'untested').values
k1 = (D1['Primary_endpoint'] != 'Excluded (Other)').values
k3 = (D3['rv1830_status'] == 'ok').values
A = pd.DataFrame({'isolate_id':D1['Isolate_ID'].astype(str),'coll':'C1','sev':D1['Phenotype'].map(m1),
     'subs':[s1(s,m) for s,m in zip(D1['Rv1830_substitutions'].fillna('-'), D1['Rv1830_mutant'])]})[k1]
B = pd.DataFrame({'isolate_id':D2['UNIQUEID'].astype(str),'coll':'C2','sev':D2['resistance_class'].map(m2),
     'subs':[t if t is not None else [] for t in p2]})[k2]
Cc = pd.DataFrame({'isolate_id':D3['sra_id'].astype(str),'coll':'C3','sev':D3['type_of_resistance'].map(m3),
     'subs':[[t.strip() for t in re.split(r'[;,]', str(s)) if t.strip() and t.strip() != '-'] if pd.notna(s) and str(s) != '-' else [] for s in D3['rv1830_aa_changes']]})[k3]
ALL = pd.concat([A, B, Cc], ignore_index=True)
ALL = ALL[ALL['sev'].notna()].copy(); ALL['sev'] = ALL['sev'].astype(int)
ALL['carrier'] = [len(L) > 0 for L in ALL['subs']]
ALL['resistant'] = (ALL['sev'] > 0).astype(int)
ALL['D144N'] = [any(t.startswith('D144N') for t in L) for L in ALL['subs']]
assert len(ALL) == 40276, f"expected 40276, got {len(ALL)}"

# ---- analysis set ----
pd.DataFrame({'isolate_id':ALL['isolate_id'],'collection':ALL['coll'],'severity_level':ALL['sev'],
    'severity_label':ALL['sev'].map(SEV),'resistant':ALL['resistant'],
    'rv1830_carrier':ALL['carrier'].astype(int),
    'rv1830_protein_substitutions':[';'.join(L) for L in ALL['subs']],
    'D144N_carrier':[int(x) for x in ALL['D144N']]}).to_csv(f"{OUT}/analysis_set_40276.csv", index=False)

HDR = Font(bold=True, color="FFFFFF"); FILL = PatternFill("solid", fgColor="2C6FBB")
def style(ws, widths):
    for c in ws[1]:
        c.font = HDR; c.fill = FILL; c.alignment = Alignment(wrap_text=True, vertical="center", horizontal="center")
    for i, w in enumerate(widths, 1): ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = w

# ---- Table 1 ----
def orc(a,b,c,d):
    o=(a*d)/(b*c); se=math.sqrt(1/a+1/b+1/c+1/d); return o, o*math.exp(-1.96*se), o*math.exp(1.96*se)
base = ALL[ALL['sev']==0]; bc=int(base['carrier'].sum()); bn=len(base)-bc
wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Table1"
ws.append(["Clinical resistance category","Total isolates (N)","All ResR carriers n (%)","Non-D144N carriers n (%)","D144N carriers n (%)","Crude OR for carriage vs susceptible (95% CI)"])
tot = {'N':0,'a':0,'nd':0,'d':0}
for lv in [0,1,2,3,4]:
    g=ALL[ALL['sev']==lv]; N=len(g); a=int(g['carrier'].sum()); dd=int((g['carrier']&g['D144N']).sum()); nd=a-dd
    tot['N']+=N; tot['a']+=a; tot['nd']+=nd; tot['d']+=dd
    orstr = "1.00 (reference)" if lv==0 else "%.2f (%.2f–%.2f)"%orc(a,N-a,bc,bn)
    ws.append([SEV[lv],f"{N:,}",f"{a} ({100*a/N:.2f}%)",f"{nd} ({100*nd/N:.2f}%)",f"{dd} ({100*dd/N:.2f}%)",orstr])
ws.append(["Total cohort",f"{tot['N']:,}",f"{tot['a']} ({100*tot['a']/tot['N']:.2f}%)",
           f"{tot['nd']} ({100*tot['nd']/tot['N']:.2f}%)",f"{tot['d']} ({100*tot['d']/tot['N']:.2f}%)","—"])
style(ws,[32,16,20,22,18,30])
wb.save(f"{OUT}/Table1_main.xlsx")

# ---- S2 substitutions ----
cnt = collections.Counter(t for L in ALL['subs'] for t in L)
def pos(t):
    m = re.match(r'^[A-Z!](\d+)', t); return int(m.group(1)) if m else None
wb2 = openpyxl.Workbook(); ws = wb2.active; ws.title = "S2_substitutions"
ws.append(["Substitution","Position","Isolates","Region","In_DNA_binding_65_137"])
for subn, n in cnt.most_common():
    p = pos(subn)
    region = "DNA-binding HTH (65-137)" if (p and 65<=p<=137) else ("N-terminal (1-64)" if (p and p<65) else ("C-terminal (138-225)" if p else "unassigned"))
    ws.append([subn, p, n, region, "Yes" if (p and 65<=p<=137) else "No"])
style(ws,[14,10,10,26,22])
wb2.save(f"{OUT}/Supplementary_Table_S2_substitutions.xlsx")

print(f"OK: {len(ALL)} isolates, {int(ALL['carrier'].sum())} carriers, {len(cnt)} distinct substitutions")
print(f"wrote analysis_set_40276.csv, Table1_main.xlsx, Supplementary_Table_S2_substitutions.xlsx to {OUT}")
