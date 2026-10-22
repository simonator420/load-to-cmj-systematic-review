#!/usr/bin/env python3
"""Post hoc exploratory sensitivity analysis: coarser strata for the three-cohort gate (not prespecified).
Run from the repository root:  python3 scripts/run_exploratory_sensitivity.py   (writes into results/sensitivity/)."""
import json, csv, math, os, sys
import numpy as np
from scipy import stats, optimize

HERE = os.path.dirname(os.path.abspath(__file__))
V14 = os.path.join(os.path.abspath(os.path.join(HERE, "..")), "data")
OUTD = os.path.join(os.path.abspath(os.path.join(HERE, "..")), "results", "sensitivity")
os.makedirs(OUTD, exist_ok=True)
ROWS = os.path.join(V14, "08_classified_correlation_rows.json")
LINK = os.path.join(V14, "05_cohort_linkage.csv")

# ---------- data ----------
rows = [r for r in json.load(open(ROWS)) if r.get("poolable_candidate") == "yes"]
assert len(rows) == 71, len(rows)
link = {r["report_id"]: r for r in csv.DictReader(open(LINK))}
def sens_cohort(r):
    l = link.get(r["record_id"])
    if l and l["linkage_status"] == "possible_overlap_unresolved":
        return f'{r["cohort"]}__{r["record_id"]}'   # split pair: one token per report
    return r["cohort"]
def lag(r):
    v = r["lag_hours"]
    return None if v in ("", None) else float(v)
def is_ext(r): return r["exposure_domain"].startswith("external")
for r in rows:
    r["sens_cohort"] = sens_cohort(r); r["lag"] = lag(r)
    r["r"] = float(r["r_numeric"]); r["n"] = float(r["n_numeric"])

# ---------- strata (see sa_definitions.md) ----------
def inl(r, lo, hi, hi_incl=True):
    l = r["lag"]; return l is not None and l >= lo and (l <= hi if hi_incl else l < hi)
STRATA = {
 "S1":  ("any external x jump height (change|absolute) x 24-72 h [user-requested]", "moderate", lambda r: is_ext(r) and inl(r,24,72)),
 "S1c": ("any external x jump_height_change x 24-72 h", "strong", lambda r: is_ext(r) and r["outcome_class"]=="jump_height_change" and inl(r,24,72)),
 "S2":  ("any load (ext+int) x jump height x 24-72 h", "weak", lambda r: inl(r,24,72)),
 "S3":  ("any external x jump height x 0-<30 h", "weak", lambda r: is_ext(r) and inl(r,0,30,False)),
 "S4":  ("any external x jump height x all known lags <=168 h", "weak", lambda r: is_ext(r) and inl(r,0,168)),
 "S5":  ("any load x jump height x all known lags <=168 h", "weak", lambda r: inl(r,0,168)),
}
MID = {"S1":48,"S1c":48,"S2":48,"S3":15,"S4":84,"S5":84}
PRIO = ["total_distance","player_load","high_speed_running","sprint_distance","sprint_number","sprint_count"]
def cprio(c):
    for i,p in enumerate(PRIO):
        if c==p: return i
    return len(PRIO) + (0 if c.startswith("accel_decel") else 1)
def tier(r):
    d = r["dependence_status"]
    if d == "within_athlete_modelled": return 1
    if d == "independent_athletes" and r["outcome_class"]=="jump_height_change": return 3
    return 4

# ---------- meta-analysis ----------
def reml(y, v, tol=1e-12, it=10000):
    k=len(y); t2=max(np.var(y,ddof=1)-np.mean(v),0.0) if k>1 else 0.0
    for _ in range(it):
        w=1/(v+t2); mu=np.sum(w*y)/np.sum(w)
        new=max(0.0, np.sum(w**2*((y-mu)**2 - v))/np.sum(w**2) + 1/np.sum(w))  # REML fixed-point (Raudenbush)
        # standard REML update: sum w^2((y-mu)^2 - v)/sum w^2 + 1/sum w
        if abs(new-t2)<tol: t2=new; break
        t2=new
    return t2
def meta(y, v):
    y=np.asarray(y,float); v=np.asarray(v,float); k=len(y)
    out={"k":k}
    if k<2: return out
    t2=reml(y,v); w=1/(v+t2); mu=np.sum(w*y)/np.sum(w)
    q_re=np.sum(w*(y-mu)**2)
    se_hk=math.sqrt(max(q_re/(k-1),0)/np.sum(w))
    tq=stats.t.ppf(0.975,k-1)
    wf=1/v; muf=np.sum(wf*y)/np.sum(wf); Q=np.sum(wf*(y-muf)**2)
    I2=max(0.0,(Q-(k-1))/Q)*100 if Q>0 else 0.0
    out.update(tau2=t2, mu_z=mu, se_hk=se_hk, ci_z=(mu-tq*se_hk, mu+tq*se_hk), Q=Q, I2=I2)
    if k>=3:
        se_re=math.sqrt(1/np.sum(w)); tp=stats.t.ppf(0.975,k-2); s=math.sqrt(t2+se_re**2)
        out["pi_z"]=(mu-tp*s, mu+tp*s)
    else: out["pi_z"]=None
    out["r"]=math.tanh(mu); out["ci_r"]=tuple(math.tanh(x) for x in out["ci_z"])
    out["pi_r"]=tuple(math.tanh(x) for x in out["pi_z"]) if out["pi_z"] else None
    return out

def crosscheck():
    """Independent checks: (1) REML tau2 by direct bounded optimisation of REML log-lik;
    (2) Hartung-Knapp SE = statsmodels WLS bse (WLS residual scale = sum w e^2/(k-1)) and t CI."""
    import statsmodels.api as sm
    res=[]
    toys={"toy_k3":([-0.2,0.35,0.1],[0.02,0.05,0.03]),"toy_k5":([0.1,-0.4,0.3,-0.1,0.25],[0.03,0.04,0.02,0.05,0.06]),
          "toy_k4_data_like":([0.04,0.32,-0.42,-0.22],[1/7,1/38,1/31,1/9])}
    for name,(y,v) in toys.items():
        y=np.array(y);v=np.array(v);k=len(y)
        def nll(t2):
            w=1/(v+t2); mu=np.sum(w*y)/np.sum(w)
            return 0.5*(np.sum(np.log(v+t2))+np.log(np.sum(w))+np.sum(w*(y-mu)**2))
        o=optimize.minimize_scalar(nll,bounds=(0,10),method="bounded",options={"xatol":1e-12})
        t2_direct=o.x if nll(o.x)<nll(0) else 0.0
        m=meta(y,v); w=1/(v+m["tau2"])
        fit=sm.WLS(y,np.ones(k),weights=w).fit(use_t=True)
        ci=fit.conf_int()[0]
        res.append(dict(case=name,tau2_mine=m["tau2"],tau2_direct_REML=float(t2_direct),
                        mu_mine=m["mu_z"],mu_sm=float(fit.params[0]),
                        se_hk_mine=m["se_hk"],se_hk_statsmodels_WLS=float(fit.bse[0]),
                        ci_mine=list(m["ci_z"]),ci_sm=[float(ci[0]),float(ci[1])],
                        ok=bool(abs(m["tau2"]-t2_direct)<1e-5 and abs(m["se_hk"]-fit.bse[0])<1e-8 and abs(m["ci_z"][0]-ci[0])<1e-8)))
    return res

def select(cohort_rows, mid):
    best=min(tier(r) for r in cohort_rows)
    cand=[r for r in cohort_rows if tier(r)==best]
    key=lambda r:(cprio(r["exposure_construct"]), abs((r["lag"] or 0)-mid), r["row_id"] if isinstance(r["row_id"],int) else 10**9)
    cand.sort(key=key)
    return cand[0], cand, best

def run(coding):
    ck="cohort" if coding=="primary" else "sens_cohort"
    out={}
    for sid,(label,defens,f) in STRATA.items():
        rs=[r for r in rows if f(r)]
        coh={}
        for r in rs: coh.setdefault(r[ck],[]).append(r)
        d=dict(label=label,defensibility=defens,n_rows=len(rs),n_reports=len({r["record_id"] for r in rs}),
               n_cohorts=len(coh),cohorts=sorted(coh),gate_met=len(coh)>=3)
        if d["gate_met"] and coding=="primary":
            sel=[];excl=[]
            for c,cr in sorted(coh.items()):
                r,cand,t=select(cr,MID[sid])
                rec=dict(cohort=c,record_id=r["record_id"],row_id=r["row_id"],construct=r["exposure_construct"],domain=r["exposure_domain"],
                         outcome=r["outcome_class"],lag_hours=r["lag"],r=r["r"],n=r["n"],tier=t,dependence=r["dependence_status"],
                         z_mean_all_rows=float(np.mean([np.arctanh(x["r"]) for x in cr])),
                         r_min_tier=min(x["r"] for x in cand), r_max_tier=max(x["r"] for x in cand))
                if r["dependence_status"]=="repeated_observations_pooled":
                    rec["excluded_reason"]="only n available is pooled observations (n=%d); not an independent athlete-level n"%r["n"]; excl.append(rec)
                else: sel.append(rec)
            d["selected"]=sel; d["excluded_from_pooling"]=excl
            d["k_pooled"]=len(sel)
            if len(sel)>=3:
                y=np.array([np.arctanh(s["r"]) for s in sel]); v=np.array([1/(s["n"]-3) for s in sel])
                m=meta(y,v); d["pooled"]=m
                d["loo"]=[dict(left_out=sel[i]["cohort"],**meta(np.delete(y,i),np.delete(v,i))) for i in range(len(sel))]
                ya=np.array([s["z_mean_all_rows"] for s in sel]); d["alt_mean_of_rows"]=meta(ya,v)
                lo=np.array([np.arctanh(s["r_min_tier"]) for s in sel]); hi=np.array([np.arctanh(s["r_max_tier"]) for s in sel])
                d["alt_min_pick"]=meta(lo,v); d["alt_max_pick"]=meta(hi,v)
            else: d["pooled"]=None; d["note"]="gate met on cohort count but <3 poolable cohorts after exclusion"
        out[sid]=d
    return out

def jsonable(o):
    if isinstance(o,dict): return {k:jsonable(v) for k,v in o.items()}
    if isinstance(o,(list,tuple)): return [jsonable(v) for v in o]
    if isinstance(o,(np.floating,np.integer)): return o.item()
    return o

def forest(sid,d,path_base):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    sel=d["selected"]; m=d["pooled"]
    fig,ax=plt.subplots(figsize=(170/25.4,2.6+0.3*len(sel)),dpi=300)
    ys=list(range(len(sel)+1,1,-1))
    for yv,s in zip(ys,sel):
        z=math.atanh(s["r"]); se=1/math.sqrt(s["n"]-3)
        lo,hi=math.tanh(z-1.96*se),math.tanh(z+1.96*se)
        ax.plot([lo,hi],[yv,yv],color="#333",lw=1); ax.plot(s["r"],yv,"s",color="#1f4e79",ms=4+s["n"]**.5)
        ax.text(-2.08,yv,f'{s["record_id"]} (n={int(s["n"])}, {s["construct"].replace("_"," ")}, {s["lag_hours"]:.0f} h)',ha="left",va="center",fontsize=7)
    c=m["r"];lo,hi=m["ci_r"]
    ax.fill([lo,c,hi,c],[1,1.25,1,0.75],color="#c0504d",zorder=3)
    if m["pi_r"]: ax.plot(m["pi_r"],[0.45,0.45],color="#c0504d",lw=1.5)
    ax.text(-2.08,1,"EXPLORATORY pooled (REML, HK)",ha="left",va="center",fontsize=7,color="#c0504d",fontweight="bold")
    ax.text(-2.08,0.45,"95% prediction interval",ha="left",va="center",fontsize=7,color="#c0504d")
    ax.axvline(0,color="grey",lw=.8,ls="--"); ax.set_xlim(-2.1,1.02); ax.set_xticks([-1,-.5,0,.5,1]); ax.set_ylim(0,len(sel)+2.3)
    ax.set_yticks([]); ax.set_xlabel("Correlation r (load -> later CMJ height); bars: 95% CI of each cohort",fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title(f"{sid}: POST HOC exploratory sensitivity analysis, not a primary result (k={len(sel)})",fontsize=8,fontweight="bold")
    for s in ("top","right","left"): ax.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("png","pdf"): fig.savefig(f"{path_base}.{ext}")
    plt.close(fig)

def main():
    prim=run("primary"); sens=run("sensitivity"); chk=crosscheck()
    res={"status":"POST HOC, not prespecified","crosscheck":chk,"primary_coding":prim,"sensitivity_coding_counts":{k:{x:v[x] for x in("n_rows","n_reports","n_cohorts","cohorts","gate_met")} for k,v in sens.items()}}
    json.dump(jsonable(res),open(os.path.join(OUTD,"sa_results.json"),"w"),indent=1)
    with open(os.path.join(OUTD,"sa_results.csv"),"w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["stratum","definition","defensibility","n_rows","n_reports","n_cohorts_primary","n_cohorts_sensitivity","gate_met_primary","gate_met_sensitivity","k_pooled","pooled_r","ci_lo","ci_hi","pi_lo","pi_hi","tau2_z","I2_pct","loo_r_range","cohorts_primary","note"])
        for sid,d in prim.items():
            p=d.get("pooled"); s=sens[sid]
            loo=d.get("loo"); lr=("%.2f to %.2f"%(min(l["r"] for l in loo),max(l["r"] for l in loo))) if loo else ""
            f=lambda x:"" if x is None else round(x,3)
            w.writerow([sid,d["label"],d["defensibility"],d["n_rows"],d["n_reports"],d["n_cohorts"],s["n_cohorts"],d["gate_met"],s["gate_met"],
              d.get("k_pooled",""),f(p and p["r"]),f(p and p["ci_r"][0]),f(p and p["ci_r"][1]),f(p and p["pi_r"] and p["pi_r"][0]),f(p and p["pi_r"] and p["pi_r"][1]),f(p and p["tau2"]),f(p and p["I2"]),lr,";".join(d["cohorts"]),d.get("note","")])
    for sid,d in prim.items():
        pass
    for sid,d in prim.items():
        print(sid,d["n_rows"],d["n_reports"],"cohorts prim",d["n_cohorts"],"sens",sens[sid]["n_cohorts"],"gate",d["gate_met"],"k_pooled",d.get("k_pooled"))
        if d.get("pooled"):
            p=d["pooled"]; print("   r=%.3f CI(%.3f,%.3f) PI%s tau2=%.4f I2=%.1f"%(p["r"],*p["ci_r"],None if not p["pi_r"] else tuple(round(x,3) for x in p["pi_r"]),p["tau2"],p["I2"]))
    for c in chk: print(c["case"],c["ok"],c["tau2_mine"],c["tau2_direct_REML"],c["se_hk_mine"],c["se_hk_statsmodels_WLS"])
main()
