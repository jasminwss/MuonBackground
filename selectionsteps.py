from collections import defaultdict
from tabulate import tabulate
# ---------- module-level config variables (set via init_config) ----------#
channel = None
partial_IP_cut = None
SBTVeto = None
ana_part = None
x0 = None
x1 = None
y0 = None
y1 = None
options_tag = None
veto_thresholds = []
cut_eff_counts_per_veto = {}
counts_per_veto = {}
event_log_lines = []
xyz_groups = {}
pdg_name_map = {}
selection_steps = []
region_labels = []

def init_config(cfg_dict):
    import selectionsteps as _self
    for k, v in cfg_dict.items():
        setattr(_self, k, v)

region_labels = ['helium','LiSc','inner_wall','outer_wall','rib','UBT']

#set PID after has a reco. but this is for comparison reasons with anupama
selection_steps = [
    'DIS in region',
    'has reco candidate',
    'incl. final state PID',
    'nDoF>25 & chi2/ndf<5 & p>1GeV',
    'exactly 1 reco candidate',
    'fiducial (walls>5cm, entrance>20cm)',
    'DOCA<1cm',
    'impact parameter',
    'SBT Veto + Skipping SBT m if applied'
]

 # add combined cut as last step, with different name if not applying cuts
sbt_region_names = {'LiSc','VetoInnerWall','VetoOuterWall','VetoLongitRib'}

# ---------- setting histograms and dictionaries ----------#
counts = {reg: [0.0 for _ in selection_steps] for reg in region_labels} #counts[region][step_index] creates empty table
counts_raw = {reg: [0 for _ in selection_steps] for reg in region_labels} #same but unweighted counts
h = {} # will store the histogram objects
xy_weight_sums = defaultdict(float)  # track total weights for xy_at_UBT histograms
xy_weight_inside = defaultdict(float)  # track inside-UBT weight totals for xy_at_UBT histograms
xy_weight_outside = defaultdict(float)  # track outside-UBT weight totals for xy_at_UBT histograms
xy_weight_counts = defaultdict(int)  # track total candidate counts for xy_at_UBT histograms
xy_weight_inside_counts = defaultdict(int)  # track inside-UBT candidate counts for xy_at_UBT histograms
xy_weight_outside_counts = defaultdict(int)  # track outside-UBT candidate counts for xy_at_UBT histograms
pid_eff_rows = ['all candidates', 'ee', 'mu mu', 'e mu', 'eX', 'mu X', 'll', 'lx']
pid_eff_counts = {row: {'He': 0.0, 'SBT': 0.0} for row in pid_eff_rows}
pid_eff_event_counts = {'He': 0.0, 'SBT': 0.0}
mumu_origin_counts = defaultdict(float) #Tracks the origin of mu mu events 
mu_origin_counts = [] #Tracks the origin of remaining events after cuts and veto are applied
vtx_eff_counts = {'reconstructible': 0, 'found': 0, 'charge_misid': 0} # global counters for overall efficiency number
cut_eff_rows = ['has a reco candidate', 'good daughter', 'has reco cand + good daughters', '1reco cand', 'fiducial', 'DOCA', 'IP<10', 'IP<250', 'IP<IP(z)', 'SBT veto', 'UBT Veto', 'mass > 0.15 GeV', 'basic cuts + IP<10', 'basic cuts + IP<250', 'basic cuts + IP<IP(z)', 'basic cuts + mass > 0.15 GeV', 'basic cuts + IP<10 f', 'basic cuts + IP<250 f', 'basic cuts + IP<IP(z) f', 'basic cuts + mass > 0.15 GeV f', 'basic cuts + IP<10 + SBT veto', 'basic cuts + IP<250 + SBT veto', 'basic cuts + IP<IP(z) + SBT veto', 'basic cuts + mass > 0.15 GeV + SBT veto', 'basic cuts + IP<10 + SBT veto f', 'basic cuts + IP<250 + SBT veto f', 'basic cuts + IP<IP(z) + SBT veto f', 'basic cuts + mass > 0.15 GeV + SBT veto f']
cand_type_labels = ["ee", "mumu", "emu", "ex", "mux", "ll", "lx", "all"]
cut_eff_counts = {
    row: {
        region: {cand: 0.0 for cand in cand_type_labels}
        for region in ['He', 'SBT']
    }
    for row in cut_eff_rows
}
cut_eff_event_counts = {'He': 0.0, 'SBT': 0.0}


def record_xy_weight(key, x_val, y_val, weight):
    """Accumulate weight sums and inside/outside counts for xy_at_UBT plots."""
    xy_weight_sums[key] += weight
    xy_weight_counts[key] += 1
    inside = (x_val >= x0) and (x_val <= x1) and (y_val >= y0) and (y_val <= y1)
    if inside:
        xy_weight_inside[key] += weight
        xy_weight_inside_counts[key] += 1
    else:
        xy_weight_outside[key] += weight
        xy_weight_outside_counts[key] += 1

def persist_xy_weight_sums(outfile_base):
    """
    Collect sum of weights for xy_at_UBT histograms, print them, and save to a text file.
    The text file shares the analysis output prefix (same as the ROOT file) with suffix 'xy_weight_sums.txt'.
    """
    if not xy_weight_sums:
        return
    lines = ["Sum of weights for xy_at_UBT histograms (total / inside UBT / outside UBT):\n"]
    for key in sorted(xy_weight_sums.keys()):
        total = xy_weight_sums[key]
        inside = xy_weight_inside.get(key, 0.0)
        outside = xy_weight_outside.get(key, 0.0)
        total_count = xy_weight_counts.get(key, 0)
        inside_count = xy_weight_inside_counts.get(key, 0)
        outside_count = xy_weight_outside_counts.get(key, 0)
        lines.append(
            f"{key}: total={total} ({total_count}), inside={inside} ({inside_count}), "
            f"outside={outside} ({outside_count})\n"
        )
    txt_path = outfile_base + 'xy.txt'
    try:
        with open(txt_path, 'w') as txt_out:
            txt_out.writelines(lines)
        print(f"Saved xy_at_UBT weight sums to {txt_path}")
    except Exception as e:
        print(f"Warning: could not save xy_at_UBT weight sums to {txt_path} ({e})")

def persist_selection_table(outfile_base):
    """Write weighted selection cutflow per region to txt and print as table.
    
    If a cell in the helium or SBT(sum) row hits 0, all remaining cells in that
    row are backfilled using:
        backfilled = last_nonzero * (cut_eff_count / has_a_reco_candidate_count)
    The backfilled number is suffixed with 'f' so it is visually distinct.
    Steps without an EffCuts mapping (e.g. 'incl. final state PID') are left as 0.
    This could of course be extended but is too complicated for now and it has not been the case yet that those steps are the ones with missing counts.
    """
    if not counts:
        return

    # --- build aggregated SBT row ---
    sbt_regions = ['LiSc', 'inner_wall', 'outer_wall', 'rib']
    sbt_vals = [
        sum(counts.get(r, [0.0] * len(selection_steps))[i] for r in sbt_regions)
        for i in range(len(selection_steps))
    ]

    # --- choose IP row name based on active channel ---
    if channel == 'fullreco':
        ip_effcuts_row = 'IP<10'
    elif channel == 'partialreco':  # partialreco
        ip_effcuts_row = 'IP<IP(z)'
    if channel == 'partialreco' and partial_IP_cut is not None:
        ip_effcuts_row = f'IP<{partial_IP_cut}'
    

    # --- choose SBT veto row name based on active threshold ---
    sbt_veto_effcuts_row = 'SBT veto' if SBTVeto is not None else None

    # Map each selection_step index to the EffCuts row name (None = no mapping)
    # selection_steps indices:
    #  0 DIS in region
    #  1 has reco candidate          <- denominator anchor
    #  2 incl. final state PID       <- no EffCuts row
    #  3 nDoF>25 & chi2/ndf<5 & p>1GeV  -> good daughter
    #  4 exactly 1 reco candidate    -> 1reco cand
    #  5 fiducial                    -> fiducial
    #  6 DOCA<1cm                    -> DOCA
    #  7 impact parameter            -> IP<10 or IP<IP(z)
    #  8 SBT Veto + ...              -> SBT veto X MeV
    # (9 mass cut if present         -> no direct row)
    step_to_effcuts = {
        3: 'good daughter',
        4: '1reco cand',
        5: 'fiducial',
        6: 'DOCA',
        7: ip_effcuts_row,
    }
    if sbt_veto_effcuts_row is not None:
        step_to_effcuts[8] = sbt_veto_effcuts_row

    def backfill_row(raw_vals, region_key):
        #denom_default = cut_eff_counts['has a reco candidate'][region_key]['all']
        denom_default = cut_eff_counts['has reco cand + good daughters'][region_key]['all']
        denom_sbt     = cut_eff_counts['has reco cand + good daughters'][region_key]['all']

        step_denom = {
            3: denom_default,
            4: denom_default,
            5: denom_default,
            6: denom_default,
            7: denom_default,
            8: denom_sbt,      # SBT veto
        }

        result = []
        last_nonzero = None
        backfilling = False

        for i, v in enumerate(raw_vals):
            denom = step_denom.get(i, denom_default)  # <-- look up per step here
            if not backfilling:
                if v != 0.0:
                    last_nonzero = v
                    result.append(f"{v:.6g}")
                else:
                    effcuts_row = step_to_effcuts.get(i)
                    if effcuts_row and denom and last_nonzero is not None:
                        eff = cut_eff_counts[effcuts_row][region_key]['all'] / denom
                        filled = last_nonzero * eff
                        last_nonzero = filled
                        result.append(f"{filled:.6g}f")
                        backfilling = True
                    else:
                        result.append("0")
                        backfilling = True
            else:
                effcuts_row = step_to_effcuts.get(i)
                if effcuts_row and denom and last_nonzero is not None and last_nonzero != 0.0:
                    eff = cut_eff_counts[effcuts_row][region_key]['all'] / denom
                    filled = last_nonzero * eff
                    last_nonzero = filled
                    result.append(f"{filled:.6g}f")
                else:
                    result.append("0")
        return result

    headers = ['region'] + selection_steps
    rows = []
    any_nonzero = False

    # --- individual region rows (no backfill, just display as-is) ---
    for reg in region_labels:
        vals = counts.get(reg, [0.0] * len(selection_steps))
        any_nonzero = any_nonzero or any(v != 0 for v in vals)
        rows.append([reg] + [f"{v:.6g}" for v in vals])

    # --- SBT(sum) raw row ---
    any_nonzero = any_nonzero or any(v != 0 for v in sbt_vals)

    # --- backfilled rows for helium and SBT(sum) ---
    # Only add backfilled rows if any zero appears after a nonzero
    he_vals = counts.get('helium', [0.0] * len(selection_steps))
    he_has_zero_after_nonzero = any(
        he_vals[i] == 0.0 and any(he_vals[j] != 0.0 for j in range(i))
        for i in range(len(he_vals))
    )
    sbt_has_zero_after_nonzero = any(
        sbt_vals[i] == 0.0 and any(sbt_vals[j] != 0.0 for j in range(i))
        for i in range(len(sbt_vals))
    )

    rows.append(['SBT(sum)'] + [f"{v:.6g}" for v in sbt_vals])

    if he_has_zero_after_nonzero:
        he_backfilled = backfill_row(he_vals, 'He')
        rows.append(['helium_bf'] + he_backfilled)

    if sbt_has_zero_after_nonzero:
        sbt_backfilled = backfill_row(sbt_vals, 'SBT')
        rows.append(['SBT(sum)_bf'] + sbt_backfilled)

    if not any_nonzero:
        return

    table_str = tabulate(rows, headers=headers, tablefmt='grid')
    txt_path = outfile_base + 'Sel.txt'
    try:
        with open(txt_path, 'w') as txt_out:
            txt_out.write(table_str + '\n')
            if he_has_zero_after_nonzero or sbt_has_zero_after_nonzero:
                txt_out.write(
                    "\n[_bf rows: zeros replaced by last_nonzero * (EffCuts_count / has_reco_candidate_count); "
                    f"IP eff used: {ip_effcuts_row}, SBT veto eff used: {sbt_veto_effcuts_row}]\n"
                )
        print(f"Selection counts table to {txt_path}")
    except Exception as e:
        print(f"Warning: could not save selection counts to {txt_path} ({e})")

def persist_selection_rawtable(outfile_base):
    """Write UNweighted selection cutflow per region to txt and print as table.
    """
    if not counts_raw:
        return

    # --- build aggregated SBT row ---
    sbt_regions = ['LiSc', 'inner_wall', 'outer_wall', 'rib']
    sbt_vals = [
        sum(counts_raw.get(r, [0.0] * len(selection_steps))[i] for r in sbt_regions)
        for i in range(len(selection_steps))
    ]

    # --- choose IP row name based on active channel ---
    if channel == 'fullreco':
        ip_effcuts_row = 'IP<10'
    elif channel == 'partialreco':  # partialreco
        ip_effcuts_row = 'IP<IP(z)'
    if channel == 'partialreco' and partial_IP_cut is not None:
        ip_effcuts_row = f'IP<{partial_IP_cut}'
    

    # --- choose SBT veto row name based on active threshold ---
    sbt_veto_effcuts_row = 'SBT veto' if SBTVeto is not None else None

    # Map each selection_step index to the EffCuts row name (None = no mapping)
    # selection_steps indices:
    #  0 DIS in region
    #  1 has reco candidate          <- denominator anchor
    #  2 incl. final state PID       <- no EffCuts row
    #  3 nDoF>25 & chi2/ndf<5 & p>1GeV  -> good daughter
    #  4 exactly 1 reco candidate    -> 1reco cand
    #  5 fiducial                    -> fiducial
    #  6 DOCA<1cm                    -> DOCA
    #  7 impact parameter            -> IP<10 or IP<IP(z)
    #  8 SBT Veto + ...              -> SBT veto X MeV
    # (9 mass cut if present         -> no direct row)
    step_to_effcuts = {
        3: 'good daughter',
        4: '1reco cand',
        5: 'fiducial',
        6: 'DOCA',
        7: ip_effcuts_row,
    }
    if sbt_veto_effcuts_row is not None:
        step_to_effcuts[8] = sbt_veto_effcuts_row


    headers = ['region'] + selection_steps
    rows = []
    any_nonzero = False

    # --- individual region rows (no backfill, just display as-is) ---
    for reg in region_labels:
        vals = counts_raw.get(reg, [0.0] * len(selection_steps))
        any_nonzero = any_nonzero or any(v != 0 for v in vals)
        rows.append([reg] + [f"{v:.6g}" for v in vals])

    # --- SBT(sum) raw row ---
    any_nonzero = any_nonzero or any(v != 0 for v in sbt_vals)

    rows.append(['SBT(sum)'] + [f"{v:.6g}" for v in sbt_vals])

    table_str = tabulate(rows, headers=headers, tablefmt='grid')
    txt_path = outfile_base + 'SelRaw.txt'
    try:        
        with open(txt_path, 'w') as txt_out:
            txt_out.write(table_str + '\n')
        print(f"Selection rawcounts table to {txt_path}")
    except Exception as e:
        print(f"Warning: could not save selection rawcounts to {txt_path} ({e})")


def persist_pid_efficiencies(outfile_base):
    """Write PID efficiencies table to txt."""
    he_total_all = pid_eff_counts['all candidates']['He']
    sbt_total_all = pid_eff_counts['all candidates']['SBT']
    he_total_events = pid_eff_event_counts['He']
    sbt_total_events = pid_eff_event_counts['SBT']
    rows = []
    any_nonzero = False
    eventwise_rows = {'ee', 'mu mu', 'e mu', 'eX', 'mu X', 'll', 'lx'}
    for row in pid_eff_rows:
        he_total = pid_eff_counts[row]['He']
        sbt_total = pid_eff_counts[row]['SBT']
        any_nonzero = any_nonzero or (he_total != 0 or sbt_total != 0)
        if row in eventwise_rows:
            he_denom = he_total_events
            sbt_denom = sbt_total_events
        else:
            he_denom = he_total_all
            sbt_denom = sbt_total_all
        he_pct = (he_total / he_denom * 100.0) if he_denom else 0.0
        sbt_pct = (sbt_total / sbt_denom * 100.0) if sbt_denom else 0.0
        rows.append([row, f"{he_total:.2g}", f"{sbt_total:.2g}", f"{he_pct:.2f}", f"{sbt_pct:.2f}"])
    if not any_nonzero:
        return
    headers = ['PID final state', 'mu IS in He (total)', 'mu IS in SBT (total)', 'He %', 'SBT %']
    table_str = tabulate(rows, headers=headers, tablefmt='grid')
    txt_path = outfile_base + 'EffPID.txt'
    try:
        with open(txt_path, 'w') as txt_out:
            txt_out.write(table_str + '\n')
        print(f"PID efficiencies table to {txt_path}")
    except Exception as e:
        print(f"Warning: could not save PID efficiencies to {txt_path} ({e})")

def persist_cut_efficiencies(outfile_base):
    """Write cut efficiencies tables for all cand types to a single txt."""
    all_tables = []

    for cand_label in cand_type_labels:  # ["ee", "mumu", "ex", "mux", "ll", "lx", "all"]
        he_total_events = cut_eff_counts['has reco cand + good daughters']['He'][cand_label]
        sbt_total_events = cut_eff_counts['has reco cand + good daughters']['SBT'][cand_label]
        rows = []

        any_nonzero = False
        for row in cut_eff_rows:
            if not (row.endswith('f') or row.endswith('SBT veto')):
                he_total = cut_eff_counts[row]['He'][cand_label]
                sbt_total = cut_eff_counts[row]['SBT'][cand_label]
                any_nonzero = any_nonzero or (he_total != 0 or sbt_total != 0)
                he_denom = he_total_events
                sbt_denom = sbt_total_events
                he_pct = (he_total / he_denom * 100.0) if he_denom else 0.0
                sbt_pct = (sbt_total / sbt_denom * 100.0) if sbt_denom else 0.0
                rows.append([row, f"{he_total:.5e}", f"{sbt_total:.5e}", f"{he_pct:.2f}", f"{sbt_pct:.2f}"])
            if row == 'SBT veto':
                he_total = cut_eff_counts[row]['He'][cand_label]
                sbt_total = cut_eff_counts[row]['SBT'][cand_label]
                any_nonzero = any_nonzero or (he_total != 0 or sbt_total != 0)
                he_denom = cut_eff_counts['has reco cand + good daughters']['He'][cand_label]
                sbt_denom = cut_eff_counts['has reco cand + good daughters']['SBT'][cand_label]
                he_pct = (he_total / he_denom * 100.0) if he_denom else 0.0
                sbt_pct = (sbt_total / sbt_denom * 100.0) if sbt_denom else 0.0
                rows.append([row, f"{he_total:.5e}", f"{sbt_total:.5e}", f"{he_pct:.2f}", f"{sbt_pct:.2f}"])

        he_denom = he_total_events
        sbt_denom = sbt_total_events

        He_basic = (cut_eff_counts['good daughter']['He'][cand_label] *
                    cut_eff_counts['DOCA']['He'][cand_label] *
                    cut_eff_counts['fiducial']['He'][cand_label] *
                    cut_eff_counts['1reco cand']['He'][cand_label])
        HE_IP10      = He_basic * cut_eff_counts['IP<10']['He'][cand_label]           / (he_denom**5) if he_denom else 0.0
        HE_IP250     = He_basic * cut_eff_counts['IP<250']['He'][cand_label]          / (he_denom**5) if he_denom else 0.0
        HE_IPz       = He_basic * cut_eff_counts['IP<IP(z)']['He'][cand_label]        / (he_denom**5) if he_denom else 0.0
        HE_mass      = He_basic * cut_eff_counts['mass > 0.15 GeV']['He'][cand_label] / (he_denom**5) if he_denom else 0.0
        HE_IP10_SBT  = HE_IP10  * cut_eff_counts['SBT veto']['He'][cand_label]  / he_denom  if he_denom else 0.0
        HE_IP250_SBT = HE_IP250 * cut_eff_counts['SBT veto']['He'][cand_label]  / he_denom  if he_denom else 0.0
        HE_IPz_SBT   = HE_IPz   * cut_eff_counts['SBT veto']['He'][cand_label]  / he_denom  if he_denom else 0.0
        HE_mass_SBT  = HE_mass  * cut_eff_counts['SBT veto']['He'][cand_label]  / he_denom  if he_denom else 0.0

        SBT_basic = (cut_eff_counts['good daughter']['SBT'][cand_label] *
                     cut_eff_counts['DOCA']['SBT'][cand_label] *
                     cut_eff_counts['fiducial']['SBT'][cand_label] *
                     cut_eff_counts['1reco cand']['SBT'][cand_label])
        SBT_IP10      = SBT_basic * cut_eff_counts['IP<10']['SBT'][cand_label]           / (sbt_denom**5) if sbt_denom else 0.0
        SBT_IP250     = SBT_basic * cut_eff_counts['IP<250']['SBT'][cand_label]          / (sbt_denom**5) if sbt_denom else 0.0
        SBT_IPz       = SBT_basic * cut_eff_counts['IP<IP(z)']['SBT'][cand_label]        / (sbt_denom**5) if sbt_denom else 0.0
        SBT_mass      = SBT_basic * cut_eff_counts['mass > 0.15 GeV']['SBT'][cand_label] / (sbt_denom**5) if sbt_denom else 0.0
        SBT_IP10_SBT  = SBT_IP10  * cut_eff_counts['SBT veto']['SBT'][cand_label] / sbt_denom if sbt_denom else 0.0
        SBT_IP250_SBT = SBT_IP250 * cut_eff_counts['SBT veto']['SBT'][cand_label] / sbt_denom if sbt_denom else 0.0
        SBT_IPz_SBT   = SBT_IPz   * cut_eff_counts['SBT veto']['SBT'][cand_label] / sbt_denom if sbt_denom else 0.0
        SBT_mass_SBT  = SBT_mass  * cut_eff_counts['SBT veto']['SBT'][cand_label] / sbt_denom if sbt_denom else 0.0

        for row in ['basic cuts + IP<10 f', 'basic cuts + IP<250 f', 'basic cuts + IP<IP(z) f',
                    'basic cuts + mass > 0.15 GeV f', 'basic cuts + IP<10 + SBT veto f',
                    'basic cuts + IP<250 + SBT veto f', 'basic cuts + IP<IP(z) + SBT veto f',
                    'basic cuts + mass > 0.15 GeV + SBT veto f']:
            if row.endswith('IP<10 + SBT veto f'):              he_eff, sbt_eff = HE_IP10_SBT,  SBT_IP10_SBT
            elif row.endswith('IP<250 + SBT veto f'):           he_eff, sbt_eff = HE_IP250_SBT, SBT_IP250_SBT
            elif row.endswith('IP<IP(z) + SBT veto f'):         he_eff, sbt_eff = HE_IPz_SBT,   SBT_IPz_SBT
            elif row.endswith('mass > 0.15 GeV + SBT veto f'):  he_eff, sbt_eff = HE_mass_SBT,  SBT_mass_SBT
            elif row.endswith('IP<10 f'):                        he_eff, sbt_eff = HE_IP10,  SBT_IP10
            elif row.endswith('IP<250 f'):                       he_eff, sbt_eff = HE_IP250, SBT_IP250
            elif row.endswith('IP<IP(z) f'):                     he_eff, sbt_eff = HE_IPz,   SBT_IPz
            elif row.endswith('mass > 0.15 GeV f'):              he_eff, sbt_eff = HE_mass,  SBT_mass
            else: continue
            any_nonzero = any_nonzero or (he_eff != 0 or sbt_eff != 0)
            rows.append([row, "-", "-", f"{he_eff*100:.5e}", f"{sbt_eff*100:.5e}"])

        if any_nonzero:
            headers = ['Cut', 'surviving mu, IS in He (total)', 'surviving mu, IS in SBT (total)', 'He %', 'SBT %']
            table_str = tabulate(rows, headers=headers, tablefmt='grid')
            all_tables.append(f"=== {cand_label} ===\n{table_str}")

    if not all_tables:
        return

    txt_path = outfile_base + 'EffCuts.txt'
    try:
        with open(txt_path, 'w') as txt_out:
            txt_out.write('\n\n'.join(all_tables) + '\n')
        print(f"Cut efficiencies table written to {txt_path}")
    except Exception as e:
        print(f"Warning: could not save cut efficiencies to {txt_path} ({e})")

def persist_mumu_origin_counts(outfile_base):
    """Write mu-mu origin counts (mother PDG pairs) to txt."""
    if not mumu_origin_counts:
        return
    rows = []
    for key, total in sorted(mumu_origin_counts.items(), key=lambda item: item[1], reverse=True):
        pdg1, pdg2 = key
        name1 = pdg_name_map.get(pdg1, str(pdg1))
        name2 = pdg_name_map.get(pdg2, str(pdg2))
        rows.append([f"{pdg1}, {pdg2}", f"{name1}, {name2}", f"{total:.6g}"])
    headers = ['origin of muons (mother PDG, sorted)', 'origin of muons (names, sorted)', 'weighted count']
    table_str = tabulate(rows, headers=headers, tablefmt='grid')
    txt_path = outfile_base + 'MuMuOrigin.txt'
    try:
        with open(txt_path, 'w') as txt_out:
            txt_out.write(table_str + '\n')
        print(f"Mu-mu origin counts table to {txt_path}")
    except Exception as e:
        print(f"Warning: could not save mu-mu origin counts to {txt_path} ({e})")

def persist_mu_origin_counts(outfile_base):
    """Write muon origin list (DIS region + mother PDG pairs + mass) to txt."""
    if not mu_origin_counts:
        return
    rows = []
    for region_label, pdg1, pdg2, mass in mu_origin_counts:
        name1 = pdg_name_map.get(pdg1, str(pdg1))
        name2 = pdg_name_map.get(pdg2, str(pdg2))
        rows.append([region_label, f"{pdg1}, {pdg2}", f"{name1}, {name2}", mass])
    headers = ['DIS region', 'origin of muons (mother PDG, sorted)', 'origin of muons (names, sorted)', 'invariant mass (GeV)']
    table_str = tabulate(rows, headers=headers, tablefmt='grid')
    txt_path = outfile_base + 'MuOrigin.txt'
    try:
        with open(txt_path, 'w') as txt_out:
            txt_out.write(table_str + '\n')
        print(f"Mu origin counts table to {txt_path}")
    except Exception as e:
        print(f"Warning: could not save mu origin counts to {txt_path} ({e})")

def persist_SBT_stats(outfile_base):
    """Write infos on remaimning events about SBT hits and vertex positions to a text file."""
    if not event_log_lines:
        return
    txt_path = outfile_base + 'SBTstats.txt'
    try:
        with open(txt_path, 'w') as txt_out:
            txt_out.writelines(event_log_lines)
        print(f"Saved SBT stats to {txt_path}")
    except Exception as e:
        print(f"Warning: could not save SBT stats to {txt_path} ({e})")

def region_label_from_basename(base):
    region_map = {
        'DecayVacuum': 'helium',
        'LiSc': 'LiSc',
        'VetoInnerWall': 'inner_wall',
        'VetoOuterWall': 'outer_wall',
        'VetoLongitRib': 'rib',
        'VetoVerticalRib': 'rib',
        'glass': 'UBT'
    }
    return region_map.get(base)

def update_selection_counts(region_label, step_index,weight):
    """Accumulate weighted counts for the given region and selection step."""
    if region_label not in counts:
        return
    if step_index < 0 or step_index >= len(selection_steps):
        return
    counts[region_label][step_index] += weight

def update_selection_rawcounts(region_label, step_index):
    """Accumulate UNWEIGHTED counts for the given region and selection step."""
    if region_label not in counts:
        return
    if step_index < 0 or step_index >= len(selection_steps):
        return
    counts_raw[region_label][step_index] += 1



def _xyz_ensure(key):
    """Ensure key exists in xyz_groups."""
    if key not in xyz_groups:
        xyz_groups[key] = {
            'pts': {
                'reco': {'zx': [[], []], 'zy': [[], []], 'xy': [[], []]},
                'ip':   {'zx': [[], []], 'zy': [[], []], 'xy': [[], []]},
            },
            'has_ip': True,
            'files': set(),
        }
    return xyz_groups[key]

def _factorized_step_value(raw_vals, target_step, region_key):
    #denom_default = cut_eff_counts['has a reco candidate'][region_key]['all']
    denom_default = cut_eff_counts['has reco cand + good daughters'][region_key]['all']
    denom_sbt     = cut_eff_counts['has reco cand + good daughters'][region_key]['all']

    # which denominator to use for each step
    step_denom = {
        3: denom_default,  # good daughter
        4: denom_default,  # 1reco cand
        5: denom_default,  # fiducial
        6: denom_default,  # DOCA
        7: denom_default,  # IP
        8: denom_sbt,      # SBT veto  <-- corrected
    }

    last_nonzero = None
    backfilling  = False
    result       = 0.0

    for i in range(target_step + 1):
        v = raw_vals[i]
        denom = step_denom.get(i, denom_default)

        if not backfilling:
            if v != 0.0:
                last_nonzero = v
                result = v
            else:
                effcuts_row = step_to_effcuts.get(i)
                if effcuts_row and denom and last_nonzero is not None:
                    eff    = cut_eff_counts[effcuts_row][region_key]['all'] / denom
                    filled = last_nonzero * eff
                    last_nonzero = filled
                    result = filled
                else:
                    result = 0.0
                backfilling = True
        else:
            effcuts_row = step_to_effcuts.get(i)
            if effcuts_row and denom and last_nonzero is not None and last_nonzero != 0.0:
                eff    = cut_eff_counts[effcuts_row][region_key]['all'] / denom
                filled = last_nonzero * eff
                last_nonzero = filled
                result = filled
            else:
                result = 0.0

    return result

def make_TestVetos_plots(outfile_base):
    cand_label = "all"
    sbt_regions = ['LiSc', 'inner_wall', 'outer_wall', 'rib']
    n = len(veto_thresholds)

    x       = array('d', [float(t) for t in veto_thresholds])
    he_rem  = array('d', [0.0] * n)
    sbt_rem = array('d', [0.0] * n)
    he_eff  = array('d', [0.0] * n)
    sbt_eff = array('d', [0.0] * n)
    he_is_factorized  = []   # bool per threshold
    sbt_is_factorized = []
    he_pre_veto_vals  = []
    sbt_pre_veto_vals = []
    he_raw_pre_zero   = []   # True if real step-7 was 0
    sbt_raw_pre_zero  = []

    sbt_raw = [
        sum(counts.get(r, [0.0] * len(selection_steps))[s] for r in sbt_regions)
        for s in range(len(selection_steps))
    ]

    for i, thr in enumerate(veto_thresholds):
        he_total_events  = cut_eff_counts['has reco candidate + good daughters']['He'][cand_label]
        sbt_total_events = cut_eff_counts['has reco candidate + good daughters']['SBT'][cand_label]

        # ---- Helium --------------------------------------------------------
        he_raw       = counts.get('helium', [0.0] * len(selection_steps))
        he_dis       = _factorized_step_value(he_raw, 0, 'He')
        he_pre_veto  = _factorized_step_value(he_raw, 7, 'He')
        real_he_pre  = he_raw[7]

        he_veto_eff  = (cut_eff_counts_per_veto[thr]['SBT veto']['He'][cand_label] / he_total_events) if he_total_events else 0.0

        real_he_after = counts.get('helium', [0.0] * len(selection_steps))[8]
        he_after      = real_he_after if real_he_after != 0.0 else he_pre_veto * he_veto_eff
        factorized_he = (real_he_after == 0.0)

        he_rem[i]  = he_after
        he_eff[i]  = he_veto_eff * 100.0
        he_is_factorized.append(factorized_he)
        print(f"  thr={thr} MeV | He : pre_veto={he_pre_veto:.4g}"
              f"{'(f)' if real_he_pre == 0.0 else ''}, "
              f"veto_eff={he_veto_eff*100:.4g}%, "
              f"after={he_after:.4g}{'(f)' if factorized_he else ''}")

        # ---- SBT(sum) ------------------------------------------------------
        sbt_dis      = _factorized_step_value(sbt_raw, 0, 'SBT')
        sbt_pre_veto = _factorized_step_value(sbt_raw, 7, 'SBT')
        real_sbt_pre = sbt_raw[7]

        sbt_veto_eff  = (cut_eff_counts_per_veto[thr]['SBT veto']['SBT'][cand_label] / sbt_total_events) if sbt_total_events else 0.0

        real_sbt_after = sum(counts.get(r, [0.0] * len(selection_steps))[8] for r in sbt_regions)
        sbt_after      = real_sbt_after if real_sbt_after != 0.0 else sbt_pre_veto * sbt_veto_eff
        factorized_sbt = (real_sbt_after == 0.0)

        sbt_rem[i]  = sbt_after
        sbt_eff[i]  = sbt_veto_eff * 100.0
        sbt_is_factorized.append(factorized_sbt)
        print(f"  thr={thr} MeV | SBT: pre_veto={sbt_pre_veto:.4g}"
              f"{'(f)' if real_sbt_pre == 0.0 else ''}, "
              f"veto_eff={sbt_veto_eff*100:.4g}%, "
              f"after={sbt_after:.4g}{'(f)' if factorized_sbt else ''}")
        he_pre_veto_vals.append(he_pre_veto)
        he_raw_pre_zero.append(real_he_pre == 0.0)
        sbt_pre_veto_vals.append(sbt_pre_veto)
        sbt_raw_pre_zero.append(real_sbt_pre == 0.0)

    # ---- split each series into real vs factorized arrays ------------------
    def split_real_fact(x_all, y_all, is_fact_flags):
        """Return (x_real, y_real, x_fact, y_fact) as array('d') pairs."""
        xr, yr, xf, yf = [], [], [], []
        for xi, yi, fact in zip(x_all, y_all, is_fact_flags):
            if fact:
                xf.append(xi); yf.append(yi)
            else:
                xr.append(xi); yr.append(yi)
        to_arr = lambda lst: array('d', lst)
        return to_arr(xr), to_arr(yr), to_arr(xf), to_arr(yf)

    he_xr,  he_yr,  he_xf,  he_yf  = split_real_fact(x, he_rem,  he_is_factorized)
    sbt_xr, sbt_yr, sbt_xf, sbt_yf = split_real_fact(x, sbt_rem, sbt_is_factorized)
    he_xr2, he_yr2, he_xf2, he_yf2   = split_real_fact(x, he_eff,  he_is_factorized)
    sbt_xr2,sbt_yr2,sbt_xf2,sbt_yf2  = split_real_fact(x, sbt_eff, sbt_is_factorized)

    STAR   = 29   # ROOT filled star (5-point)
    CIRCLE = 20
    SQUARE = 21

    def make_graph(xarr, yarr, color, marker, marker_size=1.5):
        if len(xarr) == 0:
            return None
        g = ROOT.TGraph(len(xarr), xarr, yarr)
        g.SetLineColor(color);    g.SetMarkerColor(color)
        g.SetMarkerStyle(marker); g.SetMarkerSize(marker_size)
        g.SetLineWidth(2)
        return g

    # ---- Plot 1: Remaining events ------------------------------------------
    c1 = ROOT.TCanvas('c1_veto_rem', 'SBT Veto Remaining Events', 800, 600)
    c1.SetGrid()
    mg1 = ROOT.TMultiGraph()
    g_he_rem_r  = make_graph(he_xr,  he_yr,  ROOT.kBlue, CIRCLE)
    g_he_rem_f  = make_graph(he_xf,  he_yf,  ROOT.kBlue, STAR, marker_size=2.5)
    g_sbt_rem_r = make_graph(sbt_xr, sbt_yr, ROOT.kRed,  SQUARE)
    g_sbt_rem_f = make_graph(sbt_xf, sbt_yf, ROOT.kRed,  STAR, marker_size=2.5)

    for g in (g_he_rem_r, g_he_rem_f, g_sbt_rem_r, g_sbt_rem_f):
        if g: mg1.Add(g, 'LP')

    mg1.Draw('A')
    mg1.SetTitle('Remaining events after SBT veto (* = factorized);'
                 'SBT Veto Threshold [MeV];Remaining events (weighted)')

    leg1 = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
    if g_he_rem_r:  leg1.AddEntry(g_he_rem_r,  'Helium',            'LP')
    if g_he_rem_f:  leg1.AddEntry(g_he_rem_f,  'Helium (fact.)',     'LP')
    if g_sbt_rem_r: leg1.AddEntry(g_sbt_rem_r, 'SBT(sum)',           'LP')
    if g_sbt_rem_f: leg1.AddEntry(g_sbt_rem_f, 'SBT(sum) (fact.)',   'LP')
    leg1.Draw()
    #c1.SaveAs(outfile_base + 'TestVetos_rem.png')

    # ---- Plot 2: Absolute efficiency ---------------------------------------
    c2 = ROOT.TCanvas('c2_veto_eff', 'SBT Veto Efficiency', 800, 600)
    c2.SetGrid()
    mg2 = ROOT.TMultiGraph()

    g_he_eff_r  = make_graph(he_xr2,  he_yr2,  ROOT.kBlue, CIRCLE)
    g_he_eff_f  = make_graph(he_xf2,  he_yf2,  ROOT.kBlue, STAR, marker_size=2.5)
    g_sbt_eff_r = make_graph(sbt_xr2, sbt_yr2, ROOT.kRed,  SQUARE)
    g_sbt_eff_f = make_graph(sbt_xf2, sbt_yf2, ROOT.kRed,  STAR, marker_size=2.5)

    for g in (g_he_eff_r, g_he_eff_f, g_sbt_eff_r, g_sbt_eff_f):
        if g: mg2.Add(g, 'LP')

    mg2.Draw('A')
    mg2.SetTitle('SBT veto efficiency (* = factorized);'
                 'SBT Veto Threshold [MeV]; efficiency [%]')

    leg2 = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
    if g_he_eff_r:  leg2.AddEntry(g_he_eff_r,  'Helium',            'LP')
    if g_he_eff_f:  leg2.AddEntry(g_he_eff_f,  'Helium (fact.)',     'LP')
    if g_sbt_eff_r: leg2.AddEntry(g_sbt_eff_r, 'SBT(sum)',           'LP')
    if g_sbt_eff_f: leg2.AddEntry(g_sbt_eff_f, 'SBT(sum) (fact.)',   'LP')
    leg2.Draw()
    #c2.SaveAs(outfile_base + 'TestVetos_eff.png')

    # ---- save to ROOT file -------------------------------------------------
    out = ROOT.TFile(outfile_base + 'TestVetos.root', 'RECREATE')
    c1.Write(); c2.Write()
    for g, name in [
        (g_he_rem_r,  'g_he_rem_real'),  (g_he_rem_f,  'g_he_rem_fact'),
        (g_sbt_rem_r, 'g_sbt_rem_real'), (g_sbt_rem_f, 'g_sbt_rem_fact'),
        (g_he_eff_r,  'g_he_eff_real'),  (g_he_eff_f,  'g_he_eff_fact'),
        (g_sbt_eff_r, 'g_sbt_eff_real'), (g_sbt_eff_f, 'g_sbt_eff_fact'),
    ]:
        if g: g.Write(name)

    # --- save per-threshold info to txt ---
    txt_path = outfile_base + 'TestVetos.txt'
    with open(txt_path, 'w') as f:
        f.write(f"options_tag: {options_tag}\n")
        f.write(f"channel: {channel}\n")
        f.write(f"SBTVeto: {SBTVeto} MeV\n")
        f.write(f"he_reco_candidate:  {cut_eff_counts['has a reco candidate']['He'][cand_label]:.5e}\n")
        f.write(f"sbt_reco_candidate: {cut_eff_counts['has a reco candidate']['SBT'][cand_label]:.5e}\n")
        f.write('\n')
        headers = ['thr[MeV]', 'he_pre(f?)', 'he_after(f?)', 'he_veff%', 'sbt_pre(f?)', 'sbt_after(f?)', 'sbt_veff%']
        rows = []
        for i, thr in enumerate(veto_thresholds):
            rows.append([
                thr,
                f"{he_pre_veto_vals[i]:.4g}{'(f)' if he_raw_pre_zero[i] else ''}",
                f"{he_rem[i]:.4g}{'(f)' if he_is_factorized[i] else ''}",
                f"{he_eff[i]:.4g}",
                f"{sbt_pre_veto_vals[i]:.4g}{'(f)' if sbt_raw_pre_zero[i] else ''}",
                f"{sbt_rem[i]:.4g}{'(f)' if sbt_is_factorized[i] else ''}",
                f"{sbt_eff[i]:.4g}",
            ])
        f.write(tabulate(rows, headers=headers, tablefmt='grid'))
        f.write('\n')
    print(f"TestVetos info saved to {txt_path}")

    out.Close()
    print(f"TestVetos plots saved to {outfile_base}TestVetos.root")