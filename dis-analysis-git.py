import rootUtils as ut
from rootpyPickler import Unpickler
import ROOT, os
from argparse import ArgumentParser
from collections import defaultdict
from tabulate import tabulate
import yaml
from ShipGeoConfig import AttrDict

# ---------- setting up argument parser ----------#
parser = ArgumentParser()
parser.add_argument('--test', dest='testing_code', help='Run Test', required=False, action='store_true',default=False) # test option for the code
parser.add_argument('--path', dest='path', help='path to the DIS files', required=False, action='store')  # Made not required since we now can have multiple paths 
parser.add_argument('--paths', dest='paths', help='Comma-separated list of paths to the DIS files', required=False, action='store') # argument for multiple paths, you can write --paths path1,path2,... when running the sript
parser.add_argument('--trackinfo', dest='trackinfo', help='Should be activated when a print-output should be produced about all particles passing selection cuts', action='store_true', required=False) # more 
parser.add_argument('--cutSBTDIS', dest='cutSBTDIS', help='Argument to apply a cut only considering particles with DIS in SBT in a distance (cm) to decay vessel entrance', default=None, action='store', required=False)
parser.add_argument('--onlySBTDIS', dest='onlySBTDIS', help='Should be activated when only particles are considered with DIS in SBT', action='store_true', required=False)
parser.add_argument('--tag', dest='tag', action='store', default='', help='add a tag to produced file name', required=False) # argument for giving additional tag to file name
parser.add_argument('--ip_cut', dest='ip_cut', help='argument to apply a cut on the impact parameter, ip in cm', default=None, action='store', required=False)
## next arguments are not compatible, so hardwire that here
sin_cut = parser.add_mutually_exclusive_group()
sin_cut.add_argument('--single', dest='single', help='flag for single track analysis', default=False, action='store_true', required=False) # argument declaring single track analysis
sin_cut.add_argument('--cut', dest='cut', help='flag if signal selection cuts (apart from impact parameter) should be applied, needs candidate events, so not compatible with single tracks', default=False, action='store_true', required=False) # argument to apply all cuts considered for candidate events, except impact parameter (there is a plot of ip distribution)
sin2_cut = parser.add_mutually_exclusive_group()
sin2_cut.add_argument('--no_rectangle', dest='no_rectangle', help='Should be activated to NOT draw a rectangle for UBT boundaries in xy at z of UBT plots', default=False, action='store_true', required=False)   
sin2_cut.add_argument('--UBTdimension', dest='UBTdimension', help='Comma-separated x,y dimension of UBT', action='store', required=False)   

options = parser.parse_args() # with that you can then access the arguments by options.argument[dest]

ana_part = 'muon'

# ---------- define regions  ----------#

region_labels = ['helium','LiSc','inner_wall','outer_wall','rib','UBT']
selection_steps = [
    'DIS in region',
    'has reco candidate',
    'nDoF>25 & chi2/ndf<5 & p>1GeV',
    'exactly 1 reco candidate',
    'fiducial (walls>5cm, entrance>20cm)',
    'DOCA<1cm',
    'impact parameter'
]

# ---------- setting histograms and dictionaries ----------#
# counts[region][step_index] creates empty table
counts = {reg: [0.0 for _ in selection_steps] for reg in region_labels}
h = {} # will store the histogram objects
xy_weight_sums = defaultdict(float)  # track total weights for xy_at_UBT histograms
xy_weight_inside = defaultdict(float)  # track inside-UBT weight totals for xy_at_UBT histograms
xy_weight_outside = defaultdict(float)  # track outside-UBT weight totals for xy_at_UBT histograms

# ---------- setting up test option ----------#
if options.testing_code:
     directory = './test_'
     print('test option')
else:
     directory = './'
     print('no test')

# ---------- extracting information from arguments ----------#
paths = options.paths # getting multiple paths if given 
ip_cut = float(options.ip_cut) if options.ip_cut else None # ip_cut is used later but not defined 
tag = options.tag # tag for file name
cut = options.cut
trackinfo = options.trackinfo
no_rectangle = options.no_rectangle
UBTdimension = options.UBTdimension
cutSBTDIS = float(options.cutSBTDIS) if options.cutSBTDIS else None 
print("cutSBTDIS in cm:", cutSBTDIS)
onlySBTDIS = options.onlySBTDIS
print("onlySBTDIS flag:", onlySBTDIS)

sf = ''
if options.single: sf = ' (single tracks)'

options_tag = ''
if cut:
    options_tag += 'cut_'
if options.single:
    options_tag += 'single_tracks_'
if options.ip_cut:
    options_tag += f'ip_{options.ip_cut}cm_'
if paths:
    options_tag += 'multiple_paths_'
if UBTdimension:
    options_tag += f'UBTdim{UBTdimension}_'
if no_rectangle:
    options_tag += f'noRectangle_'
if cutSBTDIS:
    options_tag += f'cutSBTDIS_{cutSBTDIS}cm_'
if onlySBTDIS:
    options_tag += 'onlySBTDIS_'

# ---------- load geometry ----------#
fairship = ROOT.gSystem.Getenv("FAIRSHIP")
with open(fairship + "/geometry/veto_config_helium.yaml", "r") as file:
    config = yaml.safe_load(file)
    veto_geo = AttrDict(config)

# ---------- setting up plots ----------#
bin_coord = 100 # #bins for all coordinates
low_x = -150 # lower limit x
high_x = 150 # upper limit x (and so on)
low_y = -150
high_y = 150

if UBTdimension:
    dims = options.UBTdimension.split(',')
    if len(dims) == 2:
        x_half = float(dims[0])/2
        y_half = float(dims[1])/2
        x0, x1 = -x_half, x_half
        y0, y1 = -y_half, y_half
        print(f'Using UBT dimensions: x: {x0} to {x1} cm, y: {y0} to {y1} cm')
    else:
        print('Error: UBTdimension argument must contain two comma-separated values for x and y dimensions.')
else:
    x0, x1 = -60, 60
    y0, y1 = -148.5, 148.5

keylist_hist = []
keylist_colz = []

item_list = ['all','helium','LiSc','inner_wall','outer_wall','rib','UBT']

for item in item_list:
    ut.bookHist(h, f'Massdistr_at_UBT_{item}', f'Massdistribution of Muons passing UBT{sf} - DIS in {item}; Mass [GeV]',bin_coord,0,5)
    keylist_hist.append(f'Massdistr_at_UBT_{item}')

    ut.bookHist(h, f'xy_at_UBT_{item}', f'Muon crossing point at z plane of UBT{sf} - DIS in {item}; x [cm]; y[cm]',bin_coord,low_x,high_x,bin_coord,low_y,high_y)
    keylist_colz.append(f'xy_at_UBT_{item}')
    xy_weight_sums[f'xy_at_UBT_{item}'] = 0.0
    if not no_rectangle:
        hist = h[f'xy_at_UBT_{item}']
        box = ROOT.TBox(x0, y0, x1, y1)
        box.SetLineColor(ROOT.kPink+3)
        box.SetLineWidth(3)
        box.SetFillStyle(0)
        hist.GetListOfFunctions().Add(box)


# ---------- functions used for analysis ----------#

def helium_rhoL(event, track_index=0,
                rho_he=1.78e-4, tol=6e-5, eps=1e-4, max_hops=256): #Anupama copy paste 
    
    """
    Return helium-only rho*L (g/cm^2) for a DIS vertex.
    Uses TGeoNavigator to follow the track in both directions from the vertex.
    """

    nav = ROOT.gGeoManager.GetCurrentNavigator()
    if not nav:
        raise RuntimeError("No TGeoNavigator; load geometry first.")

    # --- get vertex & direction ---
    v = ROOT.TVector3()
    event.MCTrack[track_index].GetStartVertex(v)
    xv, yv, zv = v.X(), v.Y(), v.Z()

    px = event.MCTrack[track_index].GetPx()
    py = event.MCTrack[track_index].GetPy()
    pz = event.MCTrack[track_index].GetPz()
    pm = ROOT.TMath.Sqrt(px*px + py*py + pz*pz) or 1.0
    dx, dy, dz = px/pm, py/pm, pz/pm

    # --- helper: is this point in helium? ---
    def in_helium(x, y, z):
        nav.SetCurrentPoint(x, y, z)
        nav.FindNode()
        node = nav.GetCurrentNode()
        if not node: return False
        try:
            rho = node.GetMedium().GetMaterial().GetDensity()
        except Exception:
            return False
        return abs(rho - rho_he) <= tol
    
    # quick check that vertex is in helium (nudge if needed)
    if not (in_helium(xv, yv, zv) or
            in_helium(xv + dx*eps, yv + dy*eps, zv + dz*eps) or
            in_helium(xv - dx*eps, yv - dy*eps, zv - dz*eps)):
        return 0.0

    # --- helper: integrate helium length from point along direction ---
    def helium_len_from(x0, y0, z0, dx, dy, dz):
        nav.SetCurrentPoint(x0 + dx*eps, y0 + dy*eps, z0 + dz*eps)
        nav.SetCurrentDirection(dx, dy, dz)
        nav.FindNode()
        total, seen_he = 0.0, False
        for _ in range(max_hops):
            node = nav.GetCurrentNode()
            if not node: break
            try:
                rho = node.GetMedium().GetMaterial().GetDensity()
            except Exception:
                rho = -1
            in_he = abs(rho - rho_he) <= tol
            nav.FindNextBoundaryAndStep()
            step = nav.GetStep()
            if in_he:
                total += step
                seen_he = True
            elif seen_he:
                break
            if nav.IsOutside(): break
            cp = nav.GetCurrentPoint()
            nav.SetCurrentPoint(cp[0] + dx*eps, cp[1] + dy*eps, cp[2] + dz*eps)
        return total

    # --- integrate forward + backward from vertex ---
    L_fwd = helium_len_from(xv, yv, zv,  dx,  dy,  dz)
    L_bwd = helium_len_from(xv, yv, zv, -dx, -dy, -dz)
    L_he  = L_fwd + L_bwd

    return rho_he * L_he  # g/cm^2


def define_muon_weight(event,SHiP_running=15,w_DIS=None):
    """Calculate event weight in 15 years."""    
    
    w_mu=event.MCTrack[0].GetWeight()  #weight of the incoming muon*DIS multiplicity normalised to a full spill   sum(w_mu) = nMuons_perspill = number of muons in a spill. w_mu is not the same as N_muperspill/N_gen, where N_gen = nEvents*DISmultiplicity ( events enhanced in Pythia to increase statistics) .

    cross=event.CrossSection
    
    if w_DIS==None:
        rho_l=event.MCTrack[2].GetWeight() # What is the difference to
    else:
        rho_l=w_DIS

    N_a=6.022e+23 

    sigma_DIS=cross*1e-27*N_a #cross section cm^2 per mole
    
    nPOTinteraction     =(2.e+20)*(SHiP_running/5) #in years
    nPOTinteraction_perspill =5.e+13
    
    n_Spill  = nPOTinteraction/nPOTinteraction_perspill  #Number of Spills in SHiP running( default=5) years  
        
    weight_i = rho_l*sigma_DIS*w_mu*n_Spill 

    return weight_i   


def dist2InnerWall(part_vtx,sgeo):
    #  dist = 0#X = ROOT.TVector3()
    #X.SetXYZ(vtx_x, vtx_y, vtx_z)
    nsteps = 8
    dalpha = 2*ROOT.TMath.Pi()/nsteps
    minDistance = float("inf")#100 *u.m
    node = sgeo.FindNode(part_vtx.X(), part_vtx.Y(), part_vtx.Z())
    if not node:
        return 0
    for n in range(nsteps):
        alpha = n * dalpha
        sdir  = (ROOT.TMath.Sin(alpha),ROOT.TMath.Cos(alpha),0.)
        node = sgeo.InitTrack(part_vtx.X(), part_vtx.Y(), part_vtx.Z(), sdir[0], sdir[1], sdir[2])
        nxt = sgeo.FindNextBoundary()
        if not nxt:
            continue
        distance = sgeo.GetStep()
        minDistance = min(minDistance, distance)
    return minDistance if minDistance < 100000 else 0

def dist2Entrance(part_vtx):
    return part_vtx.Z() - veto_geo.z0

def impact_parameter(part_vtx,part_mom,ShipGeo):
  P = ROOT.TVector3(part_mom.X(),part_mom.Y(),part_mom.Z())
  Pmag=P.Mag()
  target_point = ROOT.TVector3(0, 0, ShipGeo.target.z0) 
  X = part_vtx
  t = 0
  for i in range(3):   t += P(i)/Pmag*(target_point(i)-X(i))
  dist = 0
  for i in range(3):   dist += (target_point(i)-X(i)-t*P(i)/Pmag)**2
  dist = ROOT.TMath.Sqrt(dist)
  return dist

def find_zUBT_xy(event,sgeo,z_UBT=-2497): ### remove hardcoding!!!!!
    x_DIS = event.MCTrack[0].GetStartX()
    y_DIS = event.MCTrack[0].GetStartY()
    z_DIS = event.MCTrack[0].GetStartZ()
    p_x = event.MCTrack[0].GetPx()
    p_y = event.MCTrack[0].GetPy()
    p_z = event.MCTrack[0].GetPz()
    x_UBT = x_DIS + (z_UBT-z_DIS)*(p_x/p_z)
    y_UBT = y_DIS + (z_UBT-z_DIS)*(p_y/p_z)
    return x_UBT, y_UBT

def fill_SBT_plots(event, sgeo, ShipGeo, part_vtx=None, part_mom=None, weight=1, mass=None):
    # Always fill 'all' histograms
    h['Massdistr_at_UBT_all'].Fill(mass, weight)
    x_UBT, y_UBT = find_zUBT_xy(event, sgeo)
    h['xy_at_UBT_all'].Fill(x_UBT, y_UBT, weight)
    record_xy_weight('xy_at_UBT_all', x_UBT, y_UBT, weight)

    if trackinfo:
        try:
            tid = getattr(event.MCTrack[0], "GetTrackID", lambda: None)()
            startX = event.MCTrack[0].GetStartX()
            startY = event.MCTrack[0].GetStartY()
            startZ = event.MCTrack[0].GetStartZ()
            px = event.MCTrack[0].GetPx()
            py = event.MCTrack[0].GetPy()
            pz = event.MCTrack[0].GetPz()
            cross = getattr(event, "CrossSection", None)
            w_mu = getattr(event.MCTrack[0], "GetWeight", lambda: None)()
            print(
                f"XY_FILL ALL: TrackID={tid} start=({startX:.2f},{startY:.2f},{startZ:.2f}) "
                f"p=({px:.3f},{py:.3f},{pz:.3f}) x_UBT={x_UBT:.2f} y_UBT={y_UBT:.2f} "
                f"weight={weight} w_mu={w_mu} cross={cross}"
            )
            DISx = event.MCTrack[1].GetStartX()
            DISy = event.MCTrack[1].GetStartY()
            DISz = event.MCTrack[1].GetStartZ()
            print(f"DIS interaction point: {DISx:.4f}, {DISy:.4f}, {DISz:.4f}")
            if part_vtx is not None and part_mom is not None:
                print("impact parameter:", impact_parameter(part_vtx, part_mom, ShipGeo))
        except Exception:
            print("XY_FILL ALL: (failed to print debug info)")

    # Fill region-specific histograms based on DIS interaction point
    origin_node = sgeo.FindNode(
        event.MCTrack[0].GetStartX(), event.MCTrack[0].GetStartY(), event.MCTrack[0].GetStartZ()
    ).GetName()
    baseName = origin_node.split("_")[0]
    if baseName[:4] == 'LiSc':
        baseName = 'LiSc'
    if baseName == 'VetoVerticalRib':
        baseName = 'VetoLongitRib'

    region_to_item = {
        'DecayVacuum': 'helium',
        'LiSc': 'LiSc',
        'VetoInnerWall': 'inner_wall',
        'VetoOuterWall': 'outer_wall',
        'VetoLongitRib': 'rib',
        'glass': 'UBT',
    }
    item = region_to_item.get(baseName)
    if not item:
        return
    h[f'Massdistr_at_UBT_{item}'].Fill(mass, weight)
    h[f'xy_at_UBT_{item}'].Fill(x_UBT, y_UBT, weight)
    record_xy_weight(f'xy_at_UBT_{item}', x_UBT, y_UBT, weight)

def dis_region_basename(event, sgeo):
    node = sgeo.FindNode(event.MCTrack[0].GetStartX(),
                            event.MCTrack[0].GetStartY(),
                            event.MCTrack[0].GetStartZ())
    if not node:
            return None
    base = node.GetName().split("_")[0]
    if base[:4] == 'LiSc':
            base = 'LiSc'
    if base == 'VetoVerticalRib':
            base = 'VetoLongitRib'
    return base

def record_xy_weight(key, x_val, y_val, weight):
    """Accumulate weight sums and inside/outside counts for xy_at_UBT plots."""
    xy_weight_sums[key] += weight
    inside = (x_val >= x0) and (x_val <= x1) and (y_val >= y0) and (y_val <= y1)
    if inside:
        xy_weight_inside[key] += weight
    else:
        xy_weight_outside[key] += weight

def persist_xy_weight_sums(outfile_base):
    """
    Collect sum of weights for xy_at_UBT histograms, print them, and save to a text file.
    The text file shares the analysis output prefix (same as the ROOT file) with suffix 'xy_weight_sums.txt'.
    """
    if ana_part != 'muon' or not xy_weight_sums:
        return
    print("\nSum of weights for xy_at_UBT histograms:")
    lines = ["Sum of weights for xy_at_UBT histograms (total / inside UBT / outside UBT):\n"]
    for key in sorted(xy_weight_sums.keys()):
        total = xy_weight_sums[key]
        inside = xy_weight_inside.get(key, 0.0)
        outside = xy_weight_outside.get(key, 0.0)
        print(f"  {key}: total={total}, inside={inside}, outside={outside}")
        lines.append(f"{key}: total={total}, inside={inside}, outside={outside}\n")
    txt_path = outfile_base + 'xy_weight_sums.txt'
    try:
        with open(txt_path, 'w') as txt_out:
            txt_out.writelines(lines)
        print(f"Saved xy_at_UBT weight sums to {txt_path}")
    except Exception as e:
        print(f"Warning: could not save xy_at_UBT weight sums to {txt_path} ({e})")

def persist_selection_table(outfile_base):
    """Write weighted selection cutflow per region to txt and print as table."""
    if not counts:
        return
    headers = ['region'] + selection_steps
    rows = []
    any_nonzero = False
    for reg in region_labels:
        vals = counts.get(reg, [0.0]*len(selection_steps))
        any_nonzero = any_nonzero or any(v != 0 for v in vals)
        rows.append([reg] + [f"{v:.6g}" for v in vals])
    # add aggregated SBT row (sum of individual SBT regions)
    sbt_regions = ['LiSc','inner_wall','outer_wall','rib']
    sbt_vals = [sum(counts.get(r, [0.0]*len(selection_steps))[i] for r in sbt_regions) for i in range(len(selection_steps))]
    any_nonzero = any_nonzero or any(v != 0 for v in sbt_vals)
    rows.append(['SBT(sum)'] + [f"{v:.6g}" for v in sbt_vals])
    if not any_nonzero:
        return
    table_str = tabulate(rows, headers=headers, tablefmt='grid')
    print("\nWeighted selection counts per region (cutflow):\n" + table_str)

    txt_path = outfile_base + 'selection_counts.txt'
    try:
        with open(txt_path, 'w') as txt_out:
            txt_out.write(table_str + '\\n')
        print(f"Selection counts table to {txt_path}")
    except Exception as e:
        print(f"Warning: could not save selection counts to {txt_path} ({e})")

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

def update_selection_counts(region_label, step_index, weight):
    """Accumulate weighted counts for the given region and selection step."""
    if region_label not in counts:
        return
    if step_index < 0 or step_index >= len(selection_steps):
        return
    counts[region_label][step_index] += weight

def is_in_fiducial(candidate, event, sgeo, ShipGeo):
    """Check if the candidate is within the Fiducial Volume and has hits in all four tracking stations"""

    def tracks_in_fiducial(t1, t2):
        """
        Return True if BOTH daughter tracks (t1, t2) have hits in all
        four straw‐tube stations (1, 2, 3, 4).  Return False otherwise.
        """
        required_stations = {1, 2, 3, 4}

        for track_index in (t1, t2):
            mc_id = event.fitTrack2MC[track_index]
            seen_stations = set()
            for hit in event.strawtubesPoint:
                if hit.GetTrackID() == mc_id:
                    det_id_str = str(hit.GetDetectorID())
                    station = int(det_id_str[0])
                    seen_stations.add(station)
            
            if not required_stations.issubset(seen_stations):
                return False
        return True  # both tracks are fine

    candidate_pos = ROOT.TVector3()
    candidate.GetVertex(candidate_pos)

    if candidate_pos.Z() > ShipGeo.TrackStation1.z:
        return False
    if candidate_pos.Z() < veto_geo.z0:
        return False

    vertex_node = ROOT.gGeoManager.FindNode(
        candidate_pos.X(), candidate_pos.Y(), candidate_pos.Z()
    )
    vertex_elem = vertex_node.GetVolume().GetName()
    if not vertex_elem.startswith("DecayVacuum_"):
        return False

    t1, t2 = candidate.GetDaughter(0), candidate.GetDaughter(1)
    if not tracks_in_fiducial(t1, t2):
        return False
    return True


# ---------- analysis ----------#
def main_analysis(event, sgeo, ShipGeo, rescale_fn=None, eventNr=None, counts=None):
    cat = dis_region_basename(event,sgeo)
    if cat=='DecayVacuum':
        corrected_rhoL= helium_rhoL(event)
        rhoL=corrected_rhoL
        
    else:
        rhoL=event.MCTrack[2].GetWeight()
    weight = define_muon_weight(event,SHiP_running=15, w_DIS=rhoL)

    #### reject non-SBT DIS early
    origin_node = sgeo.FindNode(event.MCTrack[0].GetStartX(),event.MCTrack[0].GetStartY(),event.MCTrack[0].GetStartZ()).GetName()
    baseName = origin_node.split("_")[0]
    if baseName[:4] == 'LiSc': baseName = 'LiSc'
    if baseName == 'VetoVerticalRib': baseName = 'VetoLongitRib'

    region_label = region_label_from_basename(baseName)
    if region_label:
        update_selection_counts(region_label, 0, weight)  # DIS in region
        if len(event.Particles) > 0:
            update_selection_counts(region_label, 1, weight)  # has reco candidate

        # check quality and fiducial cuts on candidates; use first candidate satisfying chain
        selected_candidate = None
        selected_vtx = None
        selected_mom = None
        for part in event.Particles:
            part_vtx_tmp = ROOT.TVector3()
            part.GetVertex(part_vtx_tmp)
            status1 = event.FitTracks[part.GetDaughter(0)].getFitStatus()
            status2 = event.FitTracks[part.GetDaughter(1)].getFitStatus()
            rounded_status1 = int(round(status1.getNdf()))
            rounded_status2 = int(round(status2.getNdf()))
            if rounded_status1 <= 25 or rounded_status2 <= 25:
                continue
            if status1.getChi2()/status1.getNdf() >= 5 or status2.getChi2()/status2.getNdf() >= 5:
                continue
            if event.FitTracks[part.GetDaughter(0)].getFittedState().getMom().Mag() <= 1 or event.FitTracks[part.GetDaughter(1)].getFittedState().getMom().Mag() <= 1:
                continue
            selected_candidate = part
            selected_vtx = part_vtx_tmp
            selected_mom = ROOT.TLorentzVector()
            part.Momentum(selected_mom)
            break

        if selected_candidate:
            update_selection_counts(region_label, 2, weight)  # quality cuts
            if len(event.Particles) == 1:
                update_selection_counts(region_label, 3, weight)  # exactly 1 reco candidate
                if dist2InnerWall(selected_vtx,sgeo) > 5 and dist2Entrance(selected_vtx) > 20 and is_in_fiducial(selected_candidate, event, sgeo, ShipGeo):
                    update_selection_counts(region_label, 4, weight)  # fiducial
                    if selected_candidate.GetDoca() < 1:
                        update_selection_counts(region_label, 5, weight)  # DOCA
                        if options.ip_cut:
                            if impact_parameter(selected_vtx,selected_mom,ShipGeo) <= ip_cut:
                                update_selection_counts(region_label, 6, weight)  # IP
        
    sbt_region_names = {'LiSc','VetoInnerWall','VetoOuterWall','VetoLongitRib'}
    if onlySBTDIS or cutSBTDIS:
        if baseName not in sbt_region_names:
        # skip event entirely if DIS did not occur in SBT
            return
    
    for part in event.Particles: ## This means all particle candidates that are found in reconstruction -> everything with two tracks

        if cutSBTDIS:
            # skip candidate if z-vertex is deeper than 5 m downstream
            # geometry is in **cm**, so 5 m = 500 cm
            #decayvessel_start = -2478.0  # cm
            #z_downstream_limit = decayvessel_start + cutSBTDIS  # 5 m downstream of decay vessel start      
            #if part_vtx.Z() < z_downstream_limit: #couldnt it be if dist2entrance(part_vtx) < cutSBTDIS
            if dist2Entrance(part_vtx) < cutSBTDIS:
                continue
        
        ## define a vector to get vertex positon
        part_vtx = ROOT.TVector3()
        part.GetVertex(part_vtx)

        if cut:
            # pass
            if part.GetDoca() >= 1: continue # doca
            #if part.GetMass() < 0.15: continue # invariant mass (in GeV)
            status1 = event.FitTracks[part.GetDaughter(0)].getFitStatus() # some fit status of one of the reconstructed tracks I guess
            status2 = event.FitTracks[part.GetDaughter(1)].getFitStatus() # and the other
            if status1.getNdf() <= 25 or status2.getNdf() <= 25: continue # something with degrees of freedom
            if status1.getChi2()/status1.getNdf() >= 5 or status2.getChi2()/status2.getNdf() >= 5: continue # some statistic stuff ### need to learn again about that
            if event.FitTracks[part.GetDaughter(0)].getFittedState().getMom().Mag() <= 1 or event.FitTracks[part.GetDaughter(1)].getFittedState().getMom().Mag() <= 1: continue # particle momentum
            if dist2InnerWall(part_vtx,sgeo) <= 5: continue
            if dist2Entrance(part_vtx) < 20: continue #### 100 (vacuum) or 20 (helium)!
            if is_in_fiducial(part, event, sgeo, ShipGeo) == False: continue
            if len(event.Particles) != 1: continue

            ## at the moment all cuts are applied except for the impact parameter

        ## define a vector to get momentum
        part_mom = ROOT.TLorentzVector()
        part.Momentum(part_mom)
        mass = part_mom.M()

        if options.ip_cut:
            if impact_parameter(part_vtx,part_mom,ShipGeo) > ip_cut: continue

        fill_SBT_plots(event, sgeo, ShipGeo, part_vtx, part_mom, weight, mass)



def Main_function():
    global h
    
    files = 0
    f = None
    fgeo = None
    sgeo = None
    exception_issues = {}
    
    # Create list of paths to process
    paths_to_process = []
    if options.paths:
        # allow comma or semicolon separated lists, strip whitespace, drop empty entries
        raw_paths = [p.strip() for p in options.paths.replace(';', ',').split(',')]
        paths_to_process.extend([p for p in raw_paths if p])
    if options.path:
        if options.path.strip():
            paths_to_process.append(options.path.strip())

    if not paths_to_process:
        print("No input path(s) provided. Use --path or --paths.")
        return

    # iterate over validated paths, skip ones that don't exist
    for current_path in paths_to_process:
        if not os.path.isdir(current_path):
            print(f"Warning: path does not exist or is not a directory: '{current_path}'. Skipping.")
            continue
        print(f"Processing path: {current_path}")
        
        # Process each job directory in current path
        for jobDir in os.listdir(current_path):
            try:
                inputFile = f'{current_path}/{jobDir}/ship.conical.muonDIS-TGeant4_rec.root'
                
                f = ROOT.TFile.Open(inputFile)
                tree = f.cbmsim
                
                if not sgeo:
                    geoFile = f'{current_path}/{jobDir}/geofile_full.conical.muonDIS-TGeant4.root'
                    fgeo = ROOT.TFile(geoFile)
                    upkl = Unpickler(fgeo)
                    ShipGeo = upkl.load('ShipGeo')
                    sgeo = fgeo.FAIRGeom
                
                if options.testing_code and files > 5:
                    break
                
                print(files, jobDir)
                files += 1

                for eventNr, event in enumerate(tree):
                    try:
                        main_analysis(event, sgeo, ShipGeo)
                    except Exception as e:
                        print(f'Except called for (Reason :{e})')
                        exception_issues[jobDir] = e
                        continue
                
                f.Close()
                fgeo.Close()
                
            except Exception as e:
                if f:
                    f.Close()
                if fgeo:
                    fgeo.Close()
                print(f'Except called for (Reason :{e})')
                exception_issues[jobDir] = e
                continue

    # ... rest of existing code for saving histograms ...
    output_base = directory+tag+ana_part+'_Mass_ana_'+options_tag
    persist_xy_weight_sums(output_base)
    persist_selection_table(output_base)
    for key in keylist_hist:
        h[key].SetOption('HIST')
    for key in keylist_colz:
        h[key].SetOption('COLZ')
    ut.writeHists(h, directory+tag+ana_part+'_ana_'+options_tag+'plots.root')
    print('done')

# ---------- run analysis ----------#

Main_function()
