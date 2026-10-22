#!/usr/bin/env python3
"""Reproducible effect-, report- and cohort-level attrition cascade.
Run from the repository root:  python3 scripts/build_effect_cascade.py   (writes into results/cascade/)
"""
import json, re, itertools, sys
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
D = HERE.parent/"data"
OUT = HERE.parent/"results"/"cascade"
OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------ load
eff = pd.read_csv(D/"03_effect_level_data.csv", dtype=str)           # v14: rows of included reports only (removed rows are in 03b)
lst = pd.read_csv(D/"02_report_decisions.csv")
link = pd.read_csv(D/"05_cohort_linkage.csv")
strata = pd.read_csv(D/"08_synthesis_strata.csv")
cls = pd.DataFrame(json.load(open(D/"08_classified_correlation_rows.json")))
cls["rid"] = cls.row_id.astype(int)
gate = json.load(open(D/"09_meta_analysis_gate.json"))

# ------------------------------------------------------------------ headline facts
fact = {}
fact["reports_included"] = int(lst.disposition_category.str.startswith("included").sum())
fact["reports_by_stream"] = lst[lst.disposition_category.str.startswith("included")].disposition_category.value_counts().to_dict()
inc_ids = set(lst[lst.disposition_category.str.startswith("included")].record_id)
stream = lst.set_index("record_id").disposition_category.str.replace("included_", "", regex=False)
eff["row"] = eff.effect_row_id.astype(int)
eff["included_report"] = eff.record_id.isin(inc_ids)
fact["effect_rows_archived_file"] = len(eff)
fact["effect_rows_from_excluded_records"] = int((~eff.included_report).sum())
fact["effect_rows_included"] = int(eff.included_report.sum())
fact["reports_with_rows"] = int(eff[eff.included_report].record_id.nunique())
fact["contextual_without_rows"] = fact["reports_included"] - fact["reports_with_rows"]
N_REP = fact["reports_included"]; N_ROWS = fact["effect_rows_included"]; N_REP_ROWS = fact["reports_with_rows"]
assert fact["effect_rows_from_excluded_records"] == 0 and eff.effect_row_id.is_unique
assert sum(fact["reports_by_stream"].values()) == N_REP
fact["reports_with_rows_by_stream"] = eff[eff.included_report].drop_duplicates("record_id").record_id.map(stream).value_counts().to_dict()
fact["contextual_with_rows"] = int(fact["reports_with_rows_by_stream"].get("contextual", 0))
fact["contextual_total"] = int(fact["reports_by_stream"].get("included_contextual", 0))

# ------------------------------------------------------------------ join classified correlation rows (positional row_id)
assert cls.rid.is_unique
cm = cls.set_index("rid")
j = eff.merge(cls.drop(columns=["record_id"]), left_on="row", right_on="rid", how="left")
chk = eff.set_index("row").loc[cls.rid, "record_id"].values == cls.record_id.values
assert chk.all(), "effect_row_id alignment with record_id failed"
j["in_classified"] = j.rid.notna()
inc = j[j.included_report].copy()
assert len(inc) == N_ROWS
fact["classified_rows_total"] = int(j.in_classified.sum())
fact["classified_rows_in_included_reports"] = int(inc.in_classified.sum())
fact["classified_rows_in_excluded_records"] = int((j.in_classified & ~j.included_report).sum())
fact["poolable_yes"] = int((inc.poolable_candidate == "yes").sum())
N_PRIM = fact["poolable_yes"]
assert N_PRIM == int(gate["n_effects_in_primary_set"])

# ------------------------------------------------------------------ cohort key
pm2rec = {}
for r in lst.itertuples():
    if pd.notna(r.pmid): pm2rec[int(r.pmid)] = r.record_id
def link_rec(row):
    if row.report_id in set(lst.record_id): return row.report_id
    return pm2rec.get(int(row.pmid))
link["record_id"] = link.apply(link_rec, axis=1)
assert link.record_id.notna().all()
rec2cluster = dict(zip(link.record_id, link.cohort_id))
inc["cohort_key"] = inc.record_id.map(lambda r: rec2cluster.get(r, "ind_" + r))
inc["stream"] = inc.record_id.map(stream)

# ------------------------------------------------------------------ estimand family (rows NOT in classified file) + diagnostics for classified rows
def estimand_family(s):
    t = (s or "").lower()
    regress_words = r"regress|beta|β|slope|coefficient|\bgee\b|\br2\b|r²|\bb\)"
    if re.search(r"pearson|spearman|correlation|\brho\b|kendall", t) and not re.search(regress_words, t):
        return "correlation_type_unclassified_in_source"
    if re.search(regress_words, t):
        return "regression_or_model_coefficient_or_R2"
    if re.search(r"anova|ancova|interaction|t-test|main effect|differen|change|smd|cohen|hedges|mean|effect size|\bes\b|%|vs |median|percent|compar|decrement|delta|magnitude|descriptive|responder|pre/post|trend|chi-square|pillai|f-test|mbi|likelihood", t):
        return "group_or_time_contrast"
    return "other_or_unspecified_model_or_test"
inc["estimand_family"] = inc.estimand.map(estimand_family)
nonj = inc[~inc.in_classified]
est_map = (nonj.groupby(["estimand", "estimand_family"]).size().reset_index(name="n_rows")
           .sort_values(["estimand_family", "n_rows"], ascending=[True, False]))
est_map.to_csv(OUT/"estimand_mapping.csv", index=False)
assert est_map.n_rows.sum() == len(nonj)
fact["rows_not_in_classified_file"] = int(len(nonj))

# derived outcome (only used for the 9 unclassified correlation rows and the non-sequential matrix)
def derived_outcome_is_jh(s):
    t = (s or "").lower()
    if re.search(r"flight|power|force|velocity|impulse|rsi|time|ratio|depth|strategy|concentric|eccentric|braking|propulsive|stiff|contraction|rfd|phase", t): return False
    return bool(re.search(r"height|jump height|cmj", t) and re.search(r"height", t))
inc["derived_outcome_jh"] = inc.outcome_variable.map(derived_outcome_is_jh)
cc = inc[inc.in_classified]
agree = (cc.derived_outcome_jh == cc.outcome_class.str.startswith("jump_height"))
fact["derived_outcome_agreement_with_source_on_classified_rows"] = round(float(agree.mean()), 4)
fact["derived_outcome_agreement_n"] = [int(agree.sum()), int(len(agree))]

# ------------------------------------------------------------------ flag rules for classified rows (reason_if_no + structured fields)
RULES = {
 "A": ("estimand not a usable exposure-response association (no single usable coefficient / not a correlation / pre-post or group contrast)",
       [r"not (a )?correlation", r"r2 only", r"beta/r2 not r", r"pre-post", r"\bno (single |exact )?r\b", r"\bno r reported",
        r"\br (values )?not reported", r"rho not reported", r"coefficient not reported", r"r and n not reported", r"range of r",
        r"only as range", r"\(range only\)", r"no load-cmj correlation", r"many cells", r"multiple r$", r"multiple r\b"]),
 "C": ("exposure not a preceding quantified load / timing not eligible (concurrent, change score, biomarker, lab protocol, lag unclear/ambiguous, outside window)",
       [r"concurrent", r"not (a )?(preceding )?load", r"not training/match load", r"is not load", r"change score", r"biomarker",
        r"nutritional", r"lab ", r"wellness change", r"baseline measure preceding", r"predictor", r"temporal (relation|order)",
        r"lag (direction )?(unclear|ambiguous)", r"precedence unclear", r"multi-week", r"protocol-window", r"soreness",
        r"potentiation", r"priming", r"strength is not", r"exposure composite", r"exposure is", r"exposure not", r"exposure measured"]),
 "D": ("unit of analysis / dependence (pooled repeated observations, aggregated team means, effective n not determinable, within- vs between-athlete r)",
       [r"pooled", r"dependence", r"aggregated team means", r"microcycles", r"observations", r"rmcorr", r"rm-r", r"repeated",
        r"mean of individual correlations", r"\bn (obs )?unclear", r"\bn not reported", r"n subset unclear", r"recheck n", r"tiny n",
        r"athlete-sessions", r"between-athlete", r"n obs"]),
 "G": ("other source-recorded data-integrity problem (sign ambiguity, internal inconsistency, unverified source)",
       [r"sign ambiguity", r"internally inconsistent", r"full text not verified", r"decrement sign", r"typo"]),
}
def reason_flags(reason):
    t = reason.lower(); out = {}
    for k, (_, pats) in RULES.items():
        m = [p for p in pats if re.search(p, t)]
        out[k] = m
    return out
no = inc[inc.in_classified & (inc.poolable_candidate == "no")]
rm_rows = []
for rs, g in no.groupby("reason_if_no"):
    f = reason_flags(rs)
    rm_rows.append(dict(reason_if_no=rs, n_rows=len(g), flag_A_estimand="A" if f["A"] else "", flag_C_exposure_timing="C" if f["C"] else "",
        flag_D_unit_dependence="D" if f["D"] else "", flag_G_other_integrity="G" if f["G"] else "",
        matched_patterns=" | ".join(f"{k}:{','.join(v)}" for k, v in f.items() if v)))
reason_map = pd.DataFrame(rm_rows).sort_values("n_rows", ascending=False)

# row-level flags
def row_flags(r):
    fl = {k: False for k in "ABCDEFG"}; basis = {}
    if not r.in_classified:
        fam = r.estimand_family
        if fam != "correlation_type_unclassified_in_source":
            fl["A"] = True; basis["A"] = "estimand_family:" + fam
        else:
            basis["A"] = "not_assessed(correlation-type absent from classified file)"
            if not r.derived_outcome_jh: fl["B"] = True; basis["B"] = "derived_from_outcome_variable"
        fl["B"] = fl["B"] or (not r.derived_outcome_jh)   # derived outcome (matrix only)
        if fam != "correlation_type_unclassified_in_source": basis["B"] = "derived_from_outcome_variable"
        return fl, basis
    if r.poolable_candidate == "yes": return fl, {"all": "poolable_candidate=yes in source"}
    f = reason_flags(r.reason_if_no)
    estfam = estimand_family(r.estimand)
    if f["A"]: fl["A"] = True; basis["A"] = "reason_if_no"
    if re.search(r"mann-whitney|wilcoxon", (r.estimand or "").lower()) or re.search(r"r-squared|r²", (r.estimand or "").lower()) and re.fullmatch(r".*r-squared$", (r.estimand or "").lower()):
        if not fl["A"]: fl["A"] = True; basis["A"] = "estimand_string(rank-test effect size r / R2 only)"
    if not r.outcome_class.startswith("jump_height"): fl["B"] = True; basis["B"] = "outcome_class=" + r.outcome_class
    if r.precedes_outcome != "yes": fl["C"] = True; basis["C"] = "precedes_outcome=" + r.precedes_outcome
    elif f["C"]: fl["C"] = True; basis["C"] = "reason_if_no"
    if f["D"]: fl["D"] = True; basis["D"] = "reason_if_no"
    if f["G"]: fl["G"] = True; basis["G"] = "reason_if_no"
    return fl, basis
fl_rows = inc.apply(row_flags, axis=1)
for k in "ABCDEFG": inc["fail_" + k] = [x[0][k] for x in fl_rows]
inc["flag_basis"] = [json.dumps(x[1]) for x in fl_rows]

# reason map: add counts of structured-field flags per reason (transparency)
sub = inc[inc.in_classified & (inc.poolable_candidate == "no")]
extra = sub.groupby("reason_if_no").agg(n_outcome_not_JH_by_class=("fail_B", "sum"), n_precedes_not_yes=("fail_C", "sum"),
                                       n_estimand_string_A=("fail_A", "sum")).reset_index()
reason_map = reason_map.merge(extra, on="reason_if_no")
reason_map["mapping_note"] = ("Flags A/C/D/G from regex on reason_if_no (patterns in build_cascade.py RULES); B from outcome_class; "
                              "C additionally from precedes_outcome != yes. Sequential assignment uses first failed criterion in order A,B,C,D,G,E,F.")
reason_map.to_csv(OUT/"reason_mapping.csv", index=False)

# pooled-flag consistency
assert inc.loc[inc.poolable_candidate == "yes", [f"fail_{k}" for k in "ABCDEFG"]].sum().sum() == 0
# unclassifiable
anyfail = inc[[f"fail_{k}" for k in "ABCDEFG"]].any(axis=1)
retained = (inc.poolable_candidate == "yes")
unclass = inc[~anyfail & ~retained]
fact["rows_unclassifiable"] = int(len(unclass))

# ------------------------------------------------------------------ sequential cascade
ORDER = ["A", "B", "C", "D", "G", "E", "F"]
NAMES = {
 "A": "Estimand not a usable exposure-response association",
 "B": "Outcome not CMJ jump height",
 "C": "Exposure not a preceding quantified load / timing not eligible",
 "D": "Unit of analysis or dependence not usable",
 "G": "Other source-recorded data-integrity problem",
 "E": "Not selected by effect-selection hierarchy (within cohort x construct x lag family)",
 "F": "Duplicate / overlapping cohort (linked reports)",
}
def first_failed(r, order=ORDER):
    for k in order:
        if r["fail_" + k]: return k
    return "RETAINED" if r.poolable_candidate == "yes" else "UNCLASSIFIED"
inc["first_failed"] = inc.apply(first_failed, axis=1)
counts = inc.first_failed.value_counts().to_dict()
assert sum(counts.values()) == N_ROWS
assert counts.get("RETAINED", 0) == N_PRIM and counts.get("UNCLASSIFIED", 0) == 0, counts

# E/F diagnostics among the 59 (retained but clustered / linked)
r59 = inc[inc.first_failed == "RETAINED"].copy()
r59["fam"] = list(zip(r59.cohort_key, r59.exposure_domain, r59.exposure_construct, r59.lag_window, r59.outcome_class))
famsize = r59.groupby("fam").size()
fact["within_cohort_construct_lag_outcome_families"] = int(len(famsize))
fact["rows_sharing_family_with_another_retained_row"] = int(r59.fam.map(famsize).gt(1).sum())
fact["rows_beyond_first_in_family"] = int((famsize - 1).sum())
linked_members = set(link.record_id)
fact["linked_member_reports_with_retained_rows"] = sorted(set(r59.record_id) & linked_members)
fact["linked_reports_with_both_members_retained"] = 0
for cid, g in link.groupby("cohort_id"):
    if sum(r in set(r59.record_id) for r in g.record_id) > 1: fact["linked_reports_with_both_members_retained"] += 1
# E/F row-level removals: none (documented). Linked-member rows in the 2200:
fact["rows_in_linked_reports"] = int(inc.record_id.isin(linked_members).sum())

# cascade table
steps = [("0", "All effect rows of the 422 included reports", None)]
rows = []
remaining = inc.copy()
def snap(df):  return len(df), df.record_id.nunique(), df.cohort_key.nunique()
n0, r0, c0 = snap(remaining)
rows.append(dict(step=0, criterion=f"Effect rows of the {N_REP} included reports ({N_REP_ROWS} reports with rows)", rows_removed=0, rows_remaining=n0, reports_remaining=r0, cohorts_remaining=c0))
for i, k in enumerate(ORDER, 1):
    rem = remaining[remaining.first_failed == k]
    remaining = remaining[remaining.first_failed != k]
    n, rp, co = snap(remaining)
    rows.append(dict(step=i, criterion=NAMES[k], rows_removed=len(rem), rows_remaining=n, reports_remaining=rp, cohorts_remaining=co))
assert len(remaining) == N_PRIM
cas = pd.DataFrame(rows)
assert (cas.rows_removed.sum() == N_ROWS - N_PRIM) and cas.rows_remaining.iloc[-1] == N_PRIM
assert all(cas.rows_remaining.iloc[i] == cas.rows_remaining.iloc[i-1] - cas.rows_removed.iloc[i] for i in range(1, len(cas)))
cas.to_csv(OUT/"cascade_table.csv", index=False)

# step A sub-split
A_rows = inc[inc.first_failed == "A"]
a_split = {
 "non_classified_group_or_time_contrast": int(((~A_rows.in_classified) & (A_rows.estimand_family == "group_or_time_contrast")).sum()),
 "non_classified_regression_or_model_coefficient_or_R2": int(((~A_rows.in_classified) & (A_rows.estimand_family == "regression_or_model_coefficient_or_R2")).sum()),
 "non_classified_other_or_unspecified_model_or_test": int(((~A_rows.in_classified) & (A_rows.estimand_family == "other_or_unspecified_model_or_test")).sum()),
 "classified_correlation_type_but_no_usable_coefficient_or_not_association": int(A_rows.in_classified.sum()),
}
assert sum(a_split.values()) == len(A_rows)
# B sub-split by outcome
B_rows = inc[inc.first_failed == "B"]
b_split = B_rows.groupby(B_rows.outcome_class.fillna("derived_non_JH_unclassified_in_source")).size().to_dict()

# ------------------------------------------------------------------ non-sequential matrix
mat = []
for k in "ABCDGEF":
    assessed = {"A": N_ROWS, "B": N_ROWS, "C": int(inc.in_classified.sum()), "D": int(inc.in_classified.sum()), "G": int(inc.in_classified.sum()), "E": N_PRIM, "F": N_ROWS}[k]
    n_fail = int(inc["fail_" + k].sum())
    only = int((inc["fail_" + k] & (inc[[f"fail_{x}" for x in "ABCDGEF" if x != k]].sum(axis=1) == 0)).sum())
    mat.append(dict(criterion_code=k, criterion=NAMES[k], rows_assessable=assessed, rows_failing_regardless_of_order=n_fail,
                    rows_failing_ONLY_this_criterion=only,
                    first_failed_in_main_order=int(counts.get(k, 0)),
                    note={"A":f"all {N_ROWS:,} assessed (estimand string / reason_if_no)", "B":f"{int((~inc.in_classified).sum()):,} non-classified rows use outcome_variable-derived JH flag (validated; see json)",
                          "C":f"assessable only for the classified correlation-type rows (n={int(inc.in_classified.sum())})", "D":"as C", "G":"as C",
                          "E":f"no row-level removal recorded in data; {N_PRIM} retained rows retained-but-clustered (see json)",
                          "F":"no row-level removal recorded; linked-report members contribute to the retained set as listed in the json"}[k]))
matrix = pd.DataFrame(mat)
# combination patterns
inc["fail_pattern"] = inc.apply(lambda r: "+".join(k for k in "ABCDG" if r["fail_" + k]) or ("RETAINED" if r.first_failed == "RETAINED" else "NONE"), axis=1)
pat = inc.fail_pattern.value_counts().reset_index(); pat.columns = ["failed_criteria_set", "n_rows"]
matrix.to_csv(OUT/"nonsequential_matrix.csv", index=False)
pat.to_csv(OUT/"failure_combination_patterns.csv", index=False)

# order sensitivity: all 120 permutations of A,B,C,D,G (E,F last; no rows)
perm_res = {k: [] for k in "ABCDG"}
for p in itertools.permutations("ABCDG"):
    c = inc.apply(lambda r: first_failed(r, list(p) + ["E", "F"]), axis=1).value_counts()
    assert c.get("RETAINED", 0) == N_PRIM and sum(c.values) == N_ROWS
    for k in "ABCDG": perm_res[k].append(int(c.get(k, 0)))
sens = pd.DataFrame([dict(criterion_code=k, criterion=NAMES[k], removed_in_main_order=int(counts.get(k, 0)),
                          min_over_120_orders=min(v), max_over_120_orders=max(v)) for k, v in perm_res.items()])
sens.to_csv(OUT/"order_sensitivity.csv", index=False)

# ------------------------------------------------------------------ row-level assignment
keep = ["row", "record_id", "stream", "cohort_key", "estimand", "estimand_family", "outcome_variable", "outcome_class", "precedes_outcome",
        "lag_window", "dependence_status", "poolable_candidate", "reason_if_no", "review_status", "in_classified"] + [f"fail_{k}" for k in "ABCDGEF"] + ["first_failed", "fail_pattern", "flag_basis"]
inc[keep].rename(columns={"in_classified": "in_classified_correlation_file"}).to_csv(OUT/"row_level_assignment.csv", index=False)

# ------------------------------------------------------------------ REPORT LEVEL
rep = pd.DataFrame({"record_id": sorted(inc_ids)})
rep["stream"] = rep.record_id.map(stream)
g = inc.groupby("record_id")
rep["n_rows"] = rep.record_id.map(g.size()).fillna(0).astype(int)
rep["n_retained"] = rep.record_id.map(inc[inc.first_failed == "RETAINED"].groupby("record_id").size()).fillna(0).astype(int)
rep["n_classified_JH_rows"] = rep.record_id.map(inc[inc.in_classified & inc.outcome_class.fillna("").str.startswith("jump_height")].groupby("record_id").size()).fillna(0).astype(int)
rep["n_classified_rows"] = rep.record_id.map(inc[inc.in_classified].groupby("record_id").size()).fillna(0).astype(int)
rep["n_nonclassified_derivedJH_rows"] = rep.record_id.map(inc[~inc.in_classified & inc.derived_outcome_jh].groupby("record_id").size()).fillna(0).astype(int)
def terminal(r):
    if r.n_retained > 0: return "T1 contributes >=1 of the retained primary jump-height effects"
    if r.n_classified_JH_rows > 0: return "T2 has correlation-type jump-height rows, none retained"
    if r.n_nonclassified_derivedJH_rows > 0 and r.n_classified_rows == 0: return "T3 jump-height data only as other estimand (contrast/regression), none usable"
    if r.n_classified_rows > 0 and r.n_nonclassified_derivedJH_rows == 0: return "T4 correlation-type rows for non-jump-height outcomes only"
    if r.n_nonclassified_derivedJH_rows > 0: return "T3 jump-height data only as other estimand (contrast/regression), none usable"
    if r.n_rows > 0: return "T5 other-estimand data for non-jump-height outcomes only"
    return "T6 no effect row (characterised only)"
rep["terminal_category"] = rep.apply(terminal, axis=1)
# T3 refinement: reports with classified non-JH rows AND non-classified JH rows are placed in T3 (closest to the primary stream); documented
rep.to_csv(OUT/"report_level_terminal_category.csv", index=False)
rtab = pd.crosstab(rep.terminal_category, rep.stream, margins=True)
rtab.to_csv(OUT/"report_level_terminal_summary.csv")

# T2 dominant reason (first failed among JH classified rows) for T2 reports
t2 = rep[rep.terminal_category.str.startswith("T2")].record_id
t2_first = inc[inc.record_id.isin(t2) & inc.outcome_class.fillna("").str.startswith("jump_height")].groupby("record_id").first_failed.agg(lambda s: "+".join(sorted(set(s))))
t2_summary = t2_first.value_counts().to_dict()

# ------------------------------------------------------------------ reports -> cohorts
rep["cohort_key"] = rep.record_id.map(lambda r: rec2cluster.get(r, "ind_" + r))
rep.to_csv(OUT/"report_level_terminal_category.csv", index=False)
ret_reports = rep[rep.n_retained > 0]
tok = r59.groupby("record_id").cohort.agg(lambda s: sorted(set(s)))
tok_map = {r: t[0] for r, t in tok.items()}
assert all(len(t) == 1 for t in tok.values)
n_tokens = len(set(tok_map.values()))
fact["retained_reports"] = int(len(ret_reports)); fact["retained_reports_by_stream"] = ret_reports.stream.value_counts().to_dict()
fact["retained_cohort_tokens_source"] = n_tokens
fact["retained_cohort_keys_derived"] = int(r59.cohort_key.nunique())
assert n_tokens == fact["retained_cohort_keys_derived"] == gate["effective_independent_cohort_count_corpus_wide"]
# linkage within all 422
fact["linkage_rows"] = int(len(link)); fact["linkage_clusters"] = int(link.cohort_id.nunique())
fact["linked_reports_included"] = int(link.record_id.isin(inc_ids).sum())
fact["cohorts_among_all_included_reports"] = int(len(inc_ids) - fact["linked_reports_included"] + link.cohort_id.nunique())
fact["linkage_status_clusters"] = link.groupby("cohort_id").linkage_status.first().value_counts().to_dict()
fact["retained_reports_that_are_linked_members"] = {r: rec2cluster[r] for r in ret_reports.record_id if r in rec2cluster}
fact["partner_reports_without_retained_rows"] = [r for r in link.record_id if r not in set(ret_reports.record_id) and rec2cluster[r] in {rec2cluster[x] for x in ret_reports.record_id if x in rec2cluster}]
# 280 primary -> cohorts
prim = rep[rep.stream == "primary"]
fact["primary_reports"] = int(len(prim))
fact["primary_reports_linked_members"] = int(prim.record_id.isin(link.record_id).sum())
fact["primary_report_cohorts"] = int(prim.cohort_key.nunique())
fact["primary_reports_retained"] = int((prim.n_retained > 0).sum())

# strata
r59["stratum_domain"] = r59.exposure_domain.map(lambda d: d.split("_")[0] + "_load")
s_chk = r59.groupby(["stratum_domain", "exposure_construct", "lag_window", "outcome_class"]).agg(n_rows=("row", "size"), n_coh=("cohort_key", "nunique"))
assert len(s_chk) == len(strata) and s_chk.n_rows.sum() == N_PRIM
assert strata.n_rows.sum() == N_PRIM and strata.n_independent_cohorts.max() == s_chk.n_coh.max()
fact["strata"] = int(len(strata)); fact["max_cohorts_per_stratum"] = int(s_chk.n_coh.max())
fact["strata_by_cohort_count"] = s_chk.n_coh.value_counts().sort_index().to_dict()
fact["strata_rows_per_stratum"] = {"min": int(s_chk.n_rows.min()), "max": int(s_chk.n_rows.max()), "median": float(s_chk.n_rows.median())}
fact["strata_with_single_row"] = int((s_chk.n_rows == 1).sum())
# recon table 2
fact["table2_domain"] = strata.groupby("exposure_domain").agg(strata=("n_rows", "size"), rows=("n_rows", "sum"), maxc=("n_independent_cohorts", "max")).reset_index().to_dict("records")
fact["table2_timing"] = strata.groupby("lag_window").agg(strata=("n_rows", "size"), rows=("n_rows", "sum"), maxc=("n_independent_cohorts", "max")).reset_index().to_dict("records")
fact["outcome_class_59"] = r59.outcome_class.value_counts().to_dict()
fact["dependence_status_59"] = r59.dependence_status.value_counts().to_dict()
fact["lag_window_59"] = r59.lag_window.value_counts().to_dict()
fact["retained_by_stream_rows"] = r59.stream.value_counts().to_dict()
fact["rows_per_cohort_59"] = r59.groupby("cohort").size().sort_values(ascending=False).to_dict()
fact["cohorts_of_retained_rows_with_one_report_each"] = bool(r59.groupby("cohort").record_id.nunique().max() == 1)
fact["review_status_all_rows"] = inc.review_status.value_counts().to_dict()
fact["review_status_by_first_failed"] = pd.crosstab(inc.first_failed, inc.review_status).to_dict("index")
fact["A_split"] = a_split; fact["B_split"] = b_split
fact["T2_reports_by_first_failed_of_JH_rows"] = t2_summary
fact["unique_estimand_strings_nonclassified"] = int(est_map.estimand.nunique())
fact["nonclassified_estimand_families"] = est_map.groupby("estimand_family").n_rows.sum().to_dict()

# reasons check: reason-text/outcome agreement diagnostics
sub = inc[inc.in_classified & (inc.poolable_candidate == "no")]
fact["classified_no_rows"] = int(len(sub))
fact["reason_flag_rows_C_by_reason_only_with_precedes_yes"] = int(((sub.precedes_outcome == "yes") & sub.fail_C).sum())


# ------------------------------------------------------------------ chain table (reports -> cohorts -> strata)
_rs = fact["reports_by_stream"]; _ws = fact["reports_with_rows_by_stream"]
_pc = r59.groupby("stream").record_id.nunique().to_dict()
_rows_ch = int((r59.outcome_class == "jump_height_change").sum()); _rows_ab = int((r59.outcome_class == "jump_height_absolute").sum())
_famsz = r59.groupby("fam").size()
chain = pd.DataFrame([
 dict(level="reports", stage="Included reports", n=N_REP, detail=f"{_rs.get('included_primary',0)} primary / {_rs.get('included_secondary',0)} secondary / {_rs.get('included_contextual',0)} contextual"),
 dict(level="reports", stage="Reports with >=1 effect row", n=N_REP_ROWS, detail=f"{_ws.get('primary',0)} primary / {_ws.get('secondary',0)} secondary / {_ws.get('contextual',0)} contextual; {N_REP-N_REP_ROWS} without row"),
 dict(level="reports", stage="Reports with >=1 correlation-type row in classified file", n=int(inc[inc.in_classified].record_id.nunique()), detail=f"classified file = {int(inc.in_classified.sum())} rows of included reports"),
 dict(level="reports", stage="Reports with >=1 correlation-type jump-height row", n=int(rep[(rep.n_classified_JH_rows>0)].shape[0]), detail="by stream: " + json.dumps(rep[rep.n_classified_JH_rows>0].stream.value_counts().to_dict())),
 dict(level="reports", stage=f"Reports contributing >=1 of the {N_PRIM} effects", n=int(r59.record_id.nunique()), detail=" / ".join(f"{v} {k}" for k, v in _pc.items())),
 dict(level="cohorts", stage=f"Cohort labels of those {int(r59.record_id.nunique())} reports", n=int(r59.cohort_key.nunique()), detail=f"{int(r59.record_id.nunique())} reports map to {int(r59.cohort_key.nunique())} cohort labels; linked-pair partners without retained effects: {fact['partner_reports_without_retained_rows']}"),
 dict(level="strata", stage="Effects", n=N_PRIM, detail=f"{_rows_ch} change / {_rows_ab} absolute jump height"),
 dict(level="strata", stage="Compatibility strata", n=int(len(strata)), detail=f"{int((s_chk.n_rows==1).sum())} strata hold one effect; {int(s_chk.n_rows.min())}-{int(s_chk.n_rows.max())} effects per stratum"),
 dict(level="strata", stage="Maximum independent cohorts in any stratum", n=int(s_chk.n_coh.max()), detail=", ".join(f"{v} strata with {k} cohort(s)" for k, v in sorted(fact["strata_by_cohort_count"].items())) + f"; strata with >=3 cohorts: {int((s_chk.n_coh>=3).sum())}"),
])
chain.to_csv(OUT/"chain_table.csv", index=False)
fact["reports_with_classified_JH_rows"] = int(chain.loc[3, "n"]); fact["reports_with_classified_rows"] = int(chain.loc[2, "n"])
fact["secondary_report_with_retained_rows"] = sorted(r59[r59.stream=="secondary"].record_id.unique().tolist())
fact["retained_internal_constructs"] = r59[r59.exposure_domain=="internal_load"].exposure_construct.value_counts().to_dict()
fact["retained_rows_lag_unclear_or_rolling"] = int(r59.lag_window.isin(["unclear","rolling_or_cumulative"]).sum())
fact["retained_rows_dependence_not_independent_athletes"] = int((r59.dependence_status!="independent_athletes").sum())

# ------------------------------------------------------------------ v14 supplementary: first-failed criterion by row flag / report stream
if "v14_row_action" in inc.columns:
    pd.crosstab(inc.first_failed, inc.v14_row_action).to_csv(OUT/"firstfailed_by_v14_row_action.csv")
    pd.crosstab(inc.first_failed, inc.stream).to_csv(OUT/"firstfailed_by_report_stream.csv")
    fact["rows_flagged_remove"] = int((inc.v14_row_action == "flagged_remove_treatment_or_non_CMJ_row").sum())
    fact["rows_needs_reextraction"] = int((inc.needs_reextraction == "yes").sum())
    fact["rows_eligible_primary_secondary_flag_yes"] = int((inc.v14_eligible_primary_secondary_synthesis == "yes").sum())
    fact["rows_net_of_flagged_remove"] = int(len(inc) - fact["rows_flagged_remove"])
# ------------------------------------------------------------------ JSON
cas_out = cas.to_dict("records")
res = dict(facts=fact, cascade_order=[dict(code=k, name=NAMES[k]) for k in ORDER], cascade_table=cas_out,
           first_failed_counts={k: int(v) for k, v in counts.items()},
           nonsequential_matrix=matrix.to_dict("records"), order_sensitivity=sens.to_dict("records"),
           report_level={"all_included": rtab.to_dict(), "primary_stream": rep[rep.stream == "primary"].terminal_category.value_counts().to_dict()})
def conv(o):
    import numpy as np
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, (set, tuple)): return list(o)
    return str(o)
json.dump(res, open(OUT/"cascade_numbers.json", "w"), indent=2, default=conv)
print(cas.to_string(index=False)); print(matrix.drop(columns="note").to_string(index=False)); print(sens.to_string(index=False))
print(rtab.to_string()); print(json.dumps(fact, indent=1, default=conv)[:6000])
print(reason_map[["n_rows","flag_A_estimand","flag_C_exposure_timing","flag_D_unit_dependence","flag_G_other_integrity"]].describe())
print("UNCLASSIFIED rows:", len(unclass)); print(pat.head(25).to_string())
