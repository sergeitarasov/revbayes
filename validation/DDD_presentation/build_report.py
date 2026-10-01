#!/usr/bin/env python3
"""Create a landscape, presentation-ready PDF from saved numerical results."""
import csv
import hashlib
import json
from pathlib import Path
import re
from xml.sax.saxutils import escape
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph, Table, TableStyle
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.colors import HexColor, white
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'validation/DDD_presentation'
DATA = BASE/'results'
FIG = BASE/'figures'
OUT = ROOT/'output/pdf/ddd_simulation_report.pdf'
OUT.parent.mkdir(parents=True, exist_ok=True)
fonts = Path('/Users/taravser/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/libreoffice-headless/libreoffice/LibreOfficeDev.app/Contents/Resources/fonts/truetype')
# On other systems, set DDD_FONT_DIR to a directory with DejaVuSans*.ttf.
import os
fonts = Path(os.environ.get('DDD_FONT_DIR', fonts))
pdfmetrics.registerFont(TTFont('DV', str(fonts/'DejaVuSans.ttf')))
pdfmetrics.registerFont(TTFont('DV-Bold', str(fonts/'DejaVuSans-Bold.ttf')))
pdfmetrics.registerFont(TTFont('DV-Mono', str(fonts/'DejaVuSansMono.ttf')))
pdfmetrics.registerFontFamily('DV',normal='DV',bold='DV-Bold')

def read(name):
    return list(csv.DictReader((DATA/name).open(), delimiter='\t'))
sims=read('simulations.tsv'); comp=read('likelihood_comparisons.tsv'); fits=read('fits.tsv')
groups=read('group_summary.tsv'); meta=read('illustrative_metadata.tsv')[0]
selected=next(r for r in fits if r['id']==meta['id'])
verification=(DATA/'Rev_verification.txt').read_text()
match=re.search(r'SIMULATION_REV_CHECKS_PASSED count=(\d+) max_error=([-+.\deE]+)',verification)
if not match or re.search(r'\b(Error|Exception):',verification):
    raise RuntimeError('Actual Rev verification did not pass')
assert len(sims)==100 and len(comp)==300 and len(fits)==100 and int(match[1])==481
max_fixed=max(abs(float(r['error'])) for r in comp)
max_fit=max(abs(float(r['lnL_difference'])) for r in fits)
max_lambda=max(abs(float(r['lambda_difference'])) for r in fits)
max_cap=max(abs(float(r['cutoff_delta'])) for r in fits)
lower=sum(r['boundary']=='lower' for r in fits); upper=sum(r['boundary']=='upper' for r in fits)
assert max_fixed<1e-7 and max_fit<1e-7 and max_cap<1e-7

def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
manifest=dict(status='passed',date='2026-10-01',trees=len(sims),parameter_point_checks=len(comp),
    actual_Rev_checks=int(match[1]),actual_Rev_max_error=float(match[2]),max_parameter_point_error=max_fixed,
    max_fitted_lnL_difference=max_fit,max_fitted_lambda_difference=max_lambda,max_fit_cutoff_change=max_cap,
    lower_bound_fits=lower,upper_bound_fits=upper,reference='DDD 5.2.5',
    kernel_sha256=digest(ROOT/'src/core/functions/phylogenetics/DiversityDependentFbdLikelihood.cpp'),
    rb_sha256=digest(ROOT/'.local-build/revbayes-build/rb'),
    inputs_sha256={p.name:digest(p) for p in DATA.glob('*.tsv')},
    scripts_sha256={p.name:digest(p) for p in BASE.iterdir() if p.suffix in ('.R','.py')},
    scope='Crown-conditioned complete extant trees. Only lambda0 fitted; mu and K fixed. Pilot, not joint recovery or FBD validation.')
(DATA/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

W,H=960,540
INK=HexColor('#183047'); MUTED=HexColor('#526779'); BLUE=HexColor('#126A91'); PALE=HexColor('#ECF3F7')
ORANGE=HexColor('#C56724'); GREEN=HexColor('#298168')
c=canvas.Canvas(str(OUT),pagesize=(W,H))
c.setTitle('Diversity dependence in RevBayes: DDD benchmark and simulation pilot')
c.setAuthor('RevBayes diversity-dependent birth-death development')
c.setSubject('100 simulations, 481 Rev checks, conditional one-parameter fitting, numerical convergence')
page=0

def para(text,x,y,width,size=13,leading=None,colour=INK,font='DV'):
    style=ParagraphStyle('p',fontName=font,fontSize=size,leading=leading or size*1.35,textColor=colour)
    p=Paragraph(text,style)
    _,height=p.wrap(width,1000)
    if y-height<25: raise ValueError(f'Text extends into footer on page {page}: {text[:70]}')
    p.drawOn(c,x,y-height)
    return y-height

def start(title,subtitle=''):
    global page
    page+=1
    c.setFillColor(white);c.rect(0,0,W,H,fill=1,stroke=0)
    c.setFillColor(BLUE);c.rect(0,H-7,W,7,fill=1,stroke=0)
    para('REVBAYES  /  DIVERSITY-DEPENDENT BIRTH-DEATH',36,514,760,size=9,colour=MUTED,font='DV-Bold')
    para(title,36,487,890,size=27,leading=32,font='DV-Bold')
    if subtitle: para(subtitle,36,447,890,size=12.5,colour=MUTED)

def finish():
    c.setStrokeColor(HexColor('#DDE6EC'));c.line(36,28,924,28)
    c.setFont('DV',8);c.setFillColor(MUTED)
    c.drawString(36,15,'Simulation pilot | 1 October 2026 | DDD 5.2.5 | Complete extant sampling')
    c.drawRightString(924,15,f'{page:02d}')
    c.showPage()

def figure(name,y=100,h=325):
    path=FIG/(name+'.png')
    with Image.open(path) as im: iw,ih=im.size
    scale=min(888/iw,h/ih); dw,dh=iw*scale,ih*scale
    c.drawImage(str(path),36+(888-dw)/2,y+(h-dh)/2,width=dw,height=dh,mask='auto')

def caption(text): para(text,42,85,876,size=11.5,leading=15)

start('Diversity dependence in RevBayes','A tree-only likelihood benchmark against DDD, followed by a simulation pilot')
metrics=[('100','simulated trees'),('481','checks in the Rev executable'),('5.2 × 10⁻¹²','largest Rev-versus-DDD difference')]
for i,(number,label) in enumerate(metrics):
    x=36+i*300
    c.setFillColor(PALE);c.roundRect(x,300,286,113,9,fill=1,stroke=0)
    para(number,x+17,389,255,size=31,font='DV-Bold',colour=BLUE)
    para(label,x+17,343,255,size=12,colour=MUTED)
para('<b>Result:</b> the implementations agree numerically under matched rate laws, sampling and conditioning. Parameter estimates still vary substantially among simulated trees.',40,276,872,size=18,leading=25)
blocks=[('GENERATE','Two rate models; crown ages 3 and 8; 25 trees per scenario. λ₀ = 0.8, μ = 0.2, K = 15.'),
 ('COMPARE','Evaluate identical parameter points and fitted points. Check the actual Rev function, not just a standalone prototype.'),
 ('ESTIMATE','Fit λ₀ only, holding μ and K at their generating values. Retain all trees and all boundary estimates.')]
for i,(title,body) in enumerate(blocks):
    x=40+i*300;para(title,x,188,275,size=10,font='DV-Bold',colour=BLUE);para(body,x,164,266,size=12.5)
para('Scope: a numerical benchmark and a small estimation experiment. This does not establish joint parameter recovery or validate fossil-data MCMC.',40,68,868,size=11,colour=MUTED)
finish()

start('Speciation responds to total living diversity','N(t) includes both reconstructed lineages and lineages with no living sampled descendants.')
figure('01_rate_laws')
caption('The simulations use DDD’s linear model and its power-law model. The latter is called “exponential” in DDD documentation, but it is not exp(−αN). Both simulated curves satisfy λ(K) = μ; K is not a hard cap on diversity.')
finish()

start('The observed tree contains only part of the history','Extinct lineages affect speciation rates even when they leave no tip in the reconstructed tree.')
figure('02_tree_and_hidden_diversity')
caption(f'This history contains {meta["total_species"]} species in total, including {meta["extinct_species"]} extinct species and {meta["extant_tips"]} living tips. The likelihood receives only the extant tree. Illustration selected as the tree closest to the median tip count in the power-law, age-8 group; all 100 trees remain in the analyses.')
finish()

start('The same model gives the same likelihood','No arbitrary normalization offset was fitted or subtracted.')
figure('03_likelihood_agreement')
caption(f'Three parameter points per tree give a maximum absolute difference of {max_fixed:.2g}. The full Rev check includes these 300 points, 100 fitted points and an 81-point profile: 481 checks, maximum difference {float(match[2]):.2g}. All use crown survival, complete sampling and DDD’s phylogeny-density convention.')
finish()

start('Independent fits trace the same likelihood profile','Common scalar search and parameter bounds; extinction and K held fixed.')
figure('04_likelihood_profile')
caption(f'For this illustrative tree, λ̂₀ = {float(selected["Rev_lambda"]):.3f}, compared with the generating value 0.8. Across all 100 fits, the largest fitted log-likelihood difference is {max_fit:.2g}; the largest λ̂₀ difference is {max_lambda:.2g}. Fitting calls either DDD or the shared C++ kernel; final points are also checked in Rev.')
finish()

start('Estimates vary even when the software agrees','A 25-replicate-per-scenario pilot, with all accepted trees and search-bound estimates retained.')
figure('05_estimation_variation',y=101,h=326)
caption(f'{lower} fits reached the lower bound (0.20001) and {upper} reached the upper bound (8). Triangles identify these restricted-search results. They are not evidence that an unconstrained optimum lies at exactly that value. With μ and K known, this is an easier inference problem than estimating all three parameters.')
finish()

start('Pilot estimates remain variable across scenarios','Generating λ₀ = 0.8. Time units are arbitrary; rates are per lineage per unit time.')
headers=['Rate model','Crown\nage','Trees','Extant tips\nmedian [range]','λ̂₀\nmedian [IQR]','Bound hits\nlower / upper']
rows=[headers]
for g in sorted(groups,key=lambda x:(x['model'],float(x['age']))):
    rows.append(['Linear' if g['model']=='DDDlinear' else 'Power law',g['age'],g['n'],
      f"{float(g['median_tips']):g} [{g['min_tips']}-{g['max_tips']}]",
      f"{float(g['median_lambda']):.3f} [{float(g['q25']):.3f}-{float(g['q75']):.3f}]",
      f"{g['lower']} / {g['upper']}"])
t=Table(rows,colWidths=[138,80,66,185,242,165],rowHeights=[57,47,47,47,47])
t.setStyle(TableStyle([('FONTNAME',(0,0),(-1,-1),'DV'),('FONTNAME',(0,0),(-1,0),'DV-Bold'),
 ('FONTSIZE',(0,0),(-1,-1),12),('LEADING',(0,0),(-1,-1),16),('TEXTCOLOR',(0,0),(-1,-1),INK),
 ('BACKGROUND',(0,0),(-1,0),PALE),('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#F7F9FA')]),
 ('VALIGN',(0,0),(-1,-1),'MIDDLE'),('LEFTPADDING',(0,0),(-1,-1),11),('LINEBELOW',(0,0),(-1,0),1,HexColor('#CBD9E2'))]))
t.wrapOn(c,880,300);t.drawOn(c,42,171)
para('Interpretation',42,145,870,size=14,font='DV-Bold')
para('The power-law runs show wide variation and search-bound estimates at both ages. This pilot does not support a claim of reliable joint recovery, or a general claim that older trees solve the estimation problem. IQRs describe variation among simulated point estimates; they are not confidence intervals.',42,121,868,size=13,leading=18)
finish()

start('Hidden-lineage truncation must be checked','Solver tolerance controls propagation accuracy; it does not control the hidden-state cutoff error.')
figure('06_cutoff_convergence')
cut512=float((DATA/'stress_H512.txt').read_text())
cutrows=list(csv.DictReader((ROOT/'tests/test_DDD_benchmark/results/cutoff_convergence.tsv').open(),delimiter='\t'))
cut256=float(cutrows[-1]['kernel'])
caption(f'This earlier stress case deliberately uses insufficient cutoffs. Rev kills probability escaping the upper boundary; DDD’s matrix backend modifies its top diagonal. Their low-cutoff values need not match. Rev’s H=256 to H=512 change is {abs(cut512-cut256):.2g}. In the 100 new fits, H=128 to H=256 changes log likelihood by at most {max_cap:.2g}.')
finish()

start('Run the bundled example in RevBayes','A self-contained six-tip teaching example; its K = 8 differs from the simulation study’s K = 15.')
code='''tree = readTrees("examples/DDD/data/simulated.tre")[1]
lnL := fnDiversityDependentLogLikelihood(
    tree, lambda0=0.8, mu=0.2, K=8,
    rateModel="DDDpower", start="crown",
    condition="survival", density="DDDphylogeny",
    maxHiddenLineages=128, numericalTolerance=1e-13)
print(lnL)'''
c.setFillColor(PALE);c.roundRect(40,222,880,199,8,fill=1,stroke=0)
y=395
for line in code.splitlines():
    c.setFont('DV-Mono',13);c.setFillColor(INK);c.drawString(57,y,line);y-=23
para('Expected log likelihood: <b>−7.9263828301537</b>. The complete script also checks H=256, evaluates the other two rate laws and writes a λ₀ likelihood profile.',42,201,870,size=14,leading=20)
para('From the repository root:',42,137,870,size=11,colour=MUTED)
para('<font face="DV-Mono">.local-build/revbayes-build/rb examples/DDD/likelihood.Rev</font>',42,116,870,size=12)
para('Independent R reference: <font face="DV-Mono">Rscript examples/DDD/compare_in_R.R</font>',42,89,870,size=12)
para('This is a function in the development build. A deterministic lnL variable alone does not add a likelihood factor to an MCMC model.',42,57,870,size=10.5,colour=MUTED)
finish()

start('Model and likelihood conventions','Per-lineage speciation depends on total living diversity N; per-lineage extinction μ is constant.')
para('RATE LAWS',42,412,400,size=11,font='DV-Bold',colour=BLUE)
para('Linear (DDD model 1)',42,382,400,size=15,font='DV-Bold')
para('λ<sub>N</sub> = max[0, λ₀ − (λ₀ − μ)N/K]',42,352,408,size=17)
para('Power law (DDD model 2)',42,304,400,size=15,font='DV-Bold')
para('λ<sub>N</sub> = λ₀(N + 1)<super>−a</super><br/>a = log(λ₀/μ) / log(K + 1)',42,274,408,size=17,leading=27)
para('Original exponential form',42,204,400,size=15,font='DV-Bold')
para('λ<sub>N</sub> = λ₀ exp[−α(N − 1)]',42,174,408,size=17)
para('Only the first two forms are simulated here. K denotes the diversity at which λ equals μ; it is not the maximum possible number of species.',42,123,400,size=12,leading=17)
para('HIDDEN-LINEAGE CALCULATION',497,412,418,size=11,font='DV-Bold',colour=BLUE)
para('Let k be reconstructed lineages and h hidden living lineages, so N = k + h.',497,384,416,size=13,leading=18)
para('A[h, h] = −(k + h)(λ<sub>k+h</sub> + μ)<br/>A[h + 1, h] = (h + 2k)λ<sub>k+h</sub><br/>A[h − 1, h] = hμ',497,331,416,size=17,leading=32)
para('Propagate the weighted tree-embedding vector between nodes. At an observed branching, multiply by λ<sub>k+h</sub> and increase k. At the present, select h = 0 for complete extant sampling.',497,216,416,size=12.5,leading=18)
para('A crown start has k = 2 and omits the root birth factor. Condition on both crown sides surviving in the coupled diversity process. A is an embedding operator, not the unconditional diversity-count generator.',497,132,416,size=12.5,leading=18)
finish()

start('Methods, limitations and reproducibility','Saved trees, seeds, numerical results, executable checks and plotting code accompany the report.')
y=417
methods=[
 ('Simulation','DDD::dd_sim; models 1 and 2; λ₀ = 0.8, μ = 0.2, K = 15; crown ages 3 and 8; 25 replicates each. DDD retains histories in which both crown sides survive. No additional filtering by tree size. Seeds: 20261001 + model × 10000 + age × 100 + replicate.'),
 ('Likelihood and fitting','Complete extant sampling; crown survival condition; btorph=1. H=128, uniformization tolerance 10⁻¹³. Three parameter points per tree: linear λ₀ = 0.4, 0.6, 0.8; power λ₀ = 0.6, 0.8, 1.2. Fit λ₀ within [0.20001, 8], with μ and K fixed, using a 13-point logarithmic scan, local scalar optimization and explicit endpoint checks.'),
 ('Limits','Only extant trees and one fitted parameter are examined here. No fossils, character data, incomplete sampling, joint λ₀/μ/K recovery, confidence-interval coverage or fossil-model MCMC is validated. Twenty-five trees per group constitute a pilot; bound hits and finite search ranges constrain interpretation.'),
 ('Reproduce','Rscript validation/DDD_presentation/simulate_and_compare.R; run results/verify_in_Rev.Rev with rb; Rscript validation/DDD_presentation/make_figures.R; python3 validation/DDD_presentation/build_report.py. See the accompanying README for commands and dependencies.')]
for title,body in methods:
    para(title,42,y,158,size=13,font='DV-Bold',colour=BLUE)
    bottom=para(escape(body),210,y,706,size=11.5,leading=15.2)
    y=bottom-19
para('References',42,y,155,size=12,font='DV-Bold',colour=BLUE)
refs='Etienne et al. (2012), Proc. R. Soc. B 279:1300-1309. <link href="https://doi.org/10.1098/rspb.2011.1439" color="#126A91">doi:10.1098/rspb.2011.1439</link><br/>DDD package and source: <link href="https://rsetienne.github.io/DDD/" color="#126A91">rsetienne.github.io/DDD</link>. Installed reference: DDD 5.2.5; session details saved in results/session.txt.'
y=para(refs,210,y,706,size=10.5,leading=14)-14
para(f'Tested binary SHA-256: {manifest["rb_sha256"][:20]}…  |  Full hashes and numerical checks: results/manifest.json',42,y,874,size=9,colour=MUTED)
finish()
c.save()
manifest['report_pages']=page
manifest['report_sha256']=digest(OUT)
(DATA/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(OUT)
print(f'{page} pages; {len(sims)} trees; {match[1]} Rev checks')
