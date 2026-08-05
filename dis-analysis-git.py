import rootUtils as ut
from rootpyPickler import Unpickler
import ROOT, os
import geomGeant4
from argparse import ArgumentParser
from collections import defaultdict
from tabulate import tabulate
import yaml
from ShipGeoConfig import AttrDict
import dis_surviving_xyzplots as xyzplots
from array import array
from vertexeff import compute_vertexing_efficiency, persist_vertexing_efficiency, check_charge_misid, is_reconstructible_mc
import selectionsteps
from selectionsteps import (region_labels, selection_steps, sbt_region_names, h,
    counts, counts_raw, xy_weight_sums, xy_weight_inside, xy_weight_outside,
    xy_weight_counts, xy_weight_inside_counts, xy_weight_outside_counts,
    pid_eff_rows, pid_eff_counts, pid_eff_event_counts,
    mumu_origin_counts, mu_origin_counts, vtx_eff_counts,
    cut_eff_rows, cand_type_labels, cut_eff_counts, cut_eff_event_counts,
    record_xy_weight, persist_xy_weight_sums, persist_selection_table,
    persist_selection_rawtable, persist_pid_efficiencies, persist_cut_efficiencies,
    persist_mumu_origin_counts, persist_mu_origin_counts, persist_SBT_stats,
    region_label_from_basename, update_selection_counts, update_selection_rawcounts,
    _xyz_ensure, _factorized_step_value, make_TestVetos_plots
)

# ---------- setting up argument parser ----------#
parser = ArgumentParser()
parser.add_argument('--test', dest='testing_code', help='Run Test', required=False, action='store_true',default=False) # test option for the code
parser.add_argument('--path', dest='path', help='path to the DIS files', required=False, action='store')  # Made not required since we now can have multiple paths 
parser.add_argument('--paths', dest='paths', help='Comma-separated list of paths to the DIS files', required=False, action='store', default='/eos/experiment/ship/simulation/bkg/MuonDIS_2024helium/8070735/Tr,/eos/experiment/ship/simulation/bkg/MuonDIS_2024helium/8070735/SBT') # argument for multiple paths, you can write --paths path1,path2,... when running the sript
parser.add_argument('--trackinfo', dest='trackinfo', help='Should be activated when a print-output should be produced about all particles passing selection cuts', action='store_true', required=False) # more 
parser.add_argument('--cutSBTDIS', dest='cutSBTDIS', help='Argument to apply a cut only considering particles where DIS in SBT is only considered after a distance (cm) to decay vessel entrance (also these cells cannot apply veto).', default=None, action='store', required=False) #and Helium DIS can appear everywhere
parser.add_argument('--cutSBTDISiny', dest='cutSBTDISiny', help='Removes sclices of the SBT using the slice_sides_SBT function. There, no IS can appear in the SBT but also hits in there cannot Veto. ', default=None, action='store_true', required=False) #and Helium DIS can appear everywhere
parser.add_argument('--onlySBTDIS', dest='onlySBTDIS', help='Should be activated when only particles are considered with DIS in SBT', action='store_true', required=False)
parser.add_argument('--tag', dest='tag', action='store', default='', help='add a tag to produced file name', required=False) # argument for giving additional tag to file name
#parser.add_argument('--ip_cut', dest='ip_cut', help='argument to apply a cut on the impact parameter, ip in cm', default=None, action='store', required=False) ####new will be overwritten if a channel is chosen!
# does not work anymore because ip_cut is set by choosing the channel
parser.add_argument('--channel', dest='channel', help='argument to choose between fullreco and partialreco (applies ip cut automatically)', default=None, action='store', required=True)
parser.add_argument('--PID', dest='PID', help='Should be activated when PID selection should be applied', default=False, action='store_true', required=False) # argument to apply PID selection 
parser.add_argument('--SBTVeto', dest='SBTVeto', help='argument to choose between threshhold 0, 45, 90 MeV', default=None, action='store', required=False)
parser.add_argument('--dist2iWall', dest='dist2iWall', help='argument to choose from which distance to decay vessel inner wall (cm) vertices will be considerd (standard: 5cm)', default=None, action='store', required=False)
parser.add_argument('--mass_cut', dest='mass_cut', help='Should be activated when a mass cut at 0.15 GeV should be applied', default=False, action='store_true', required=False) # argument to apply a mass cut > 0.15 GeV
parser.add_argument('--partial_IP_cut', dest='partial_IP_cut', help='Argument to apply a cut on the impact parameter for the partial reco. channel, ip in cm', default=None, action='store', required=False) ####can be used instead of the dileptonic_ip_tresh for the partial reco. channel
parser.add_argument('--TestVetos', dest='TestVetos', help='Test SBT veto thresholds. Format: min,max,step (MeV). E.g. --TestVetos 40,90,5', default=None, action='store', required=False)
parser.add_argument('--PID_NO_CONFUSION', dest='PID_NO_CONFUSION', help='Should be activated to apply only PID efficiency without confusion matrix', default=False, action='store_true', required=False) # argument to apply only PID efficiency without confusion matrix
sin_cut = parser.add_mutually_exclusive_group()
sin_cut.add_argument('--single', dest='single', help='flag for single track analysis', default=False, action='store_true', required=False) # argument declaring single track analysis
sin_cut.add_argument('--cut', dest='cut', help='flag if signal selection cuts (apart from impact parameter) should be applied, needs candidate events, so not compatible with single tracks', default=False, action='store_true', required=False) # argument to apply all cuts considered for candidate events, except impact parameter (there is a plot of ip distribution)
sin2_cut = parser.add_mutually_exclusive_group()
sin2_cut.add_argument('--no_rectangle', dest='no_rectangle', help='Should be activated to NOT draw a rectangle for UBT boundaries in xy at z of UBT plots', default=False, action='store_true', required=False)   
sin2_cut.add_argument('--UBTdimension', dest='UBTdimension', help='Comma-separated x,y dimension of UBT', action='store', required=False)   


options = parser.parse_args() # with that you can then access the arguments by options.argument[dest]

ana_part = 'muon'

# ---------- define regions  ----------#

if options.mass_cut:
    selection_steps.append('mass > 0.15 GeV + SBT veto')

if options.partial_IP_cut:
    _dynamic_ip_row = f'IP<{int(options.partial_IP_cut)}'
    if _dynamic_ip_row not in cut_eff_rows:
        cut_eff_rows.append(_dynamic_ip_row)
        cut_eff_counts[_dynamic_ip_row] = {
            region: {cand: 0.0 for cand in cand_type_labels}
            for region in ['He', 'SBT']
        }


event_log_lines = []
pdg_name_map = { #maps PDG codes to names for common particles
    0: 'unknown',
    11: 'e-',
    -11: 'e+',
    13: 'mu-',
    -13: 'mu+',
    22: 'photon',
    111: 'pi0',
    211: 'pi+',
    -211: 'pi-',
    321: 'K+',
    -321: 'K-',
    2212: 'p',
    -2212: 'pbar',
    2112: 'n',
    -2112: 'nbar',
} 
xyz_groups = {}
_genfit_field_ready = False

# ---------- setting up test option ----------#
if options.testing_code:
     directory = '/afs/cern.ch/work/j/jaweiss/private/test_'
     print('test option')
else:
     directory = '/eos/user/j/jaweiss/results_muonbackground/'
     print('no test')

# ---------- extracting information from arguments ----------#
paths = options.paths # getting multiple paths if given 
#ip_cut = float(options.ip_cut) if options.ip_cut else None
tag = options.tag # tag for file name
cut = options.cut
trackinfo = options.trackinfo
no_rectangle = options.no_rectangle
UBTdimension = options.UBTdimension
cutSBTDIS = float(options.cutSBTDIS) if options.cutSBTDIS else None 
print("cutSBTDIS in cm:", cutSBTDIS)
onlySBTDIS = options.onlySBTDIS
print("onlySBTDIS flag:", onlySBTDIS)
PID = options.PID
print("Using PID selection?:", PID)
SBTVeto = int(options.SBTVeto) if options.SBTVeto else None
print("SBTVeto threshhold (MeV):", SBTVeto)
dist2iWall = float(options.dist2iWall) if options.dist2iWall else 5.0
print("Distance to inner wall cut (cm):", dist2iWall)
cutSBTDISiny = options.cutSBTDISiny
if cutSBTDISiny:
    print(r"cutSBTDISiny along $|y|<(z-z_{entrance}) \cdot 25 cm/(z_{exit}-z_{entrance})+15 cm$")
sf = ''
if options.single: sf = ' (single tracks)'
PID_NO_CONFUSION = options.PID_NO_CONFUSION
if PID_NO_CONFUSION:
    print("No confusion matrix will be applied for PID, only efficiency.")
partial_IP_cut = int(options.partial_IP_cut) if options.partial_IP_cut else None

if options.TestVetos:
    parts = options.TestVetos.split(',')
    veto_thresholds = list(range(int(parts[0]), int(parts[1]) + int(parts[2]), int(parts[2])))
    print(f"TestVetos: testing thresholds {veto_thresholds} MeV")
    cut_eff_counts_per_veto = {
        thr: {
            'SBT veto': {
                region: {cand: 0.0 for cand in cand_type_labels}
                for region in ['He', 'SBT']
            }
        }
        for thr in veto_thresholds
    }
    counts_per_veto = {thr: {reg: [0.0 for _ in selection_steps] for reg in region_labels} for thr in veto_thresholds}
else:
    veto_thresholds = []

options_tag = ''
options_tag += 'cut_'
if options.single:
    options_tag += 'single_tracks_'
if UBTdimension:
    options_tag += f'UBTdim{UBTdimension}_'
if no_rectangle:
    options_tag += f'noRectangle_'
if cutSBTDIS:
    options_tag += f'cutSBT{cutSBTDIS}cm_'
if onlySBTDIS:
    options_tag += 'onlySBT_'
if PID:
    options_tag += 'PID_'
if SBTVeto:
    options_tag += f'SBTVeto{SBTVeto}MeV_'
if dist2iWall:
    options_tag += f'dist2iWall{dist2iWall}cm_'
if options.mass_cut:
    options_tag += 'mass0.15_'
if cutSBTDISiny is not None:
    options_tag += f'cutSBTDISiny_'
if PID_NO_CONFUSION:
    options_tag += 'PIDnoConfusion_'
if partial_IP_cut:
    options_tag += f'IP{partial_IP_cut}cm_'

# ---------- setting up channel dependent parameters----------#
channel = options.channel
if channel == "partialreco":
    print(f"Partial Reco. (l l ν) channel Analysis starts now ")
    ip_cut=250
    # usually not used because the dileptonic_ip_tresh is more efficient
    # can by applied using --partial_IP_cut argument, which applies a cut on the impact parameter for the partial reco. channel, ip in cm
    # still needs to be activated because in cutflow i ask for "if ip_cut"
    # should be changed 
    finalstate='dileptonic'
    options_tag += 'dilep_'
elif channel == "fullreco":
    print(f"Fully Reco. (l π) channel Analysis starts now ")
    ip_cut=10
    finalstate='semileptonic'
    options_tag += 'semilep_'
else: 
    raise RuntimeError("Unknown channel!")

# --- unchique necesseties due to naming differences ---
if channel == 'fullreco':
    ip_effcuts_row = 'IP<10'
else:
    ip_effcuts_row = 'IP<IP(z)'

sbt_veto_effcuts_row = 'SBT veto' if SBTVeto is not None else None

step_to_effcuts = {
    3: 'good daughter',
    4: '1reco cand',
    5: 'fiducial',
    6: 'DOCA',
    7: ip_effcuts_row,
}
if sbt_veto_effcuts_row is not None:
    step_to_effcuts[8] = sbt_veto_effcuts_row

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

selectionsteps.init_config({
    'channel': channel,
    'partial_IP_cut': partial_IP_cut,
    'SBTVeto': SBTVeto,
    'ana_part': ana_part,
    'x0': x0, 'x1': x1, 'y0': y0, 'y1': y1,
    'options_tag': options_tag,
    'veto_thresholds': veto_thresholds,
    'event_log_lines': event_log_lines,
    'xyz_groups': xyz_groups,
    'pdg_name_map': pdg_name_map,
})

for item in item_list:
    ut.bookHist(h, f'Massdistr_at_UBT_{item}', f'Massdistribution of Muons passing UBT{sf} - DIS in {item}, {finalstate}; Mass [GeV]',bin_coord,0,5)
    keylist_hist.append(f'Massdistr_at_UBT_{item}')
    ut.bookHist(h, f'xy_at_UBT_{item}', f'Muon crossing point at z plane of UBT{sf} - DIS in {item}, {finalstate}; x [cm]; y[cm]',bin_coord,low_x,high_x,bin_coord,low_y,high_y)
    keylist_colz.append(f'xy_at_UBT_{item}')
    xy_weight_sums[f'xy_at_UBT_{item}'] = 0.0
    if not no_rectangle:
        hist = h[f'xy_at_UBT_{item}']
        box = ROOT.TBox(x0, y0, x1, y1)
        box.SetLineColor(ROOT.kPink+3)
        box.SetLineWidth(3)
        box.SetFillStyle(0)
        hist.GetListOfFunctions().Add(box)

ut.bookHist(h, f'Dist2WallvsVtx_z', f'Distance to Wall vs z of IS vertex; z_vtx [cm]; d [cm]',50,-2500,2500,40,0,200)

ut.bookHist(h, f'Energy_distribution', f'Energy distribution of Muons before any cuts; abs(p) [GeV]', 100,0,350)

ut.bookHist(
    h,
    'Massdistr_at_UBT_helium_inside',
    f'Massdistribution of Muons passing UBT{sf} - DIS in helium, inside UBT, {finalstate}; Mass [GeV]',
    bin_coord, 0, 5)
ut.bookHist(h, f'Energy_distribution2', f'Energy distribution of Muons before PID cut; abs(p) [GeV]', 100,0,350)

ut.bookHist(
    h,
    'Massdistr_at_UBT_helium_inside',
    f'Massdistribution of Muons passing UBT{sf} - DIS in helium, inside UBT, {finalstate}; Mass [GeV]',
    bin_coord, 0, 5)

ut.bookHist(
    h,
    'Massdistr_at_UBT_helium_outside',
    f'Massdistribution of Muons passing UBT{sf} - DIS in helium, outside UBT, {finalstate}; Mass [GeV]',
    bin_coord, 0, 5
)

ut.bookHist(h, f'x-y-IS-SBT', f'x-y-distribution of vertices in SBT before cuts; x [cm]; y[cm]',bin_coord,low_x,high_x,bin_coord,low_y,high_y)
ut.bookHist(h, f'z-IS-SBT', f'z-distribution of vertices in SBT before cuts; z [cm]; Events',50,-2500,2500)
ut.bookHist(h, f'z-IS-SBT-after-cuts', f'z-distribution of vertices in SBT after cuts; z [cm]; Events',50,-2500,2500)

# ---------- histogram setup (in the booking section at the top) ----------#
ut.bookHist(h, 'vtx_eff_z_reconstructible', 
            'Reconstructible events vs z; z_{vtx,MC} [cm]; Events', 
            50, -2500, 2500)
ut.bookHist(h, 'vtx_eff_z_found',           
            'Vertex found events vs z; z_{vtx,MC} [cm]; Events',     
            50, -2500, 2500)
ut.bookHist(h, 'vtx_eff_z_ratio',           
            'Vertexing efficiency vs z; z_{vtx,MC} [cm]; Efficiency',
            50, -2500, 2500)
ut.bookHist(h, 'charge_misid_z',            
            'Charge mis-ID events vs z; z_{vtx,MC} [cm]; Events',    
            50, -2500, 2500)



keylist_hist.append('Massdistr_at_UBT_helium_inside')
keylist_hist.append('Massdistr_at_UBT_helium_outside')
keylist_hist.append('Energy_distribution')
keylist_hist.append('Energy_distribution2')
keylist_colz.append('Dist2WallvsVtx_z')
keylist_colz.append('x-y-IS-SBT')
keylist_hist.append('z-IS-SBT')
keylist_hist.append('z-IS-SBT-after-cuts')


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

def fill_SBT_plots(event, sgeo, ShipGeo, part_vtx=None, part_mom=None, weight=1, mass=None, baseName=None):
    # Always fill 'all' histograms
    h['Massdistr_at_UBT_all'].Fill(mass, weight)
    x_UBT, y_UBT = find_zUBT_xy(event, sgeo)
    h['xy_at_UBT_all'].Fill(x_UBT, y_UBT, weight)
    record_xy_weight('xy_at_UBT_all', x_UBT, y_UBT, weight, mass=mass)



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
    #origin_node = sgeo.FindNode(
     #   event.MCTrack[0].GetStartX(), event.MCTrack[0].GetStartY(), event.MCTrack[0].GetStartZ()
    #).GetName()
    #baseName = origin_node.split("_")[0]
    #if baseName[:4] == 'LiSc':
     #   baseName = 'LiSc'
    #if baseName == 'VetoVerticalRib':
     #   baseName = 'VetoLongitRib'

    #region_to_item = {
     #   'DecayVacuum': 'helium',
      #  'LiSc': 'LiSc',
       # 'VetoInnerWall': 'inner_wall',
        #'VetoOuterWall': 'outer_wall',
        #'VetoLongitRib': 'rib',
        #'glass': 'UBT',
    #}
    item = region_label_from_basename(baseName)
    if not item:
        return

    h[f'Massdistr_at_UBT_{item}'].Fill(mass, weight)
    h[f'xy_at_UBT_{item}'].Fill(x_UBT, y_UBT, weight)
    record_xy_weight(f'xy_at_UBT_{item}', x_UBT, y_UBT, weight, mass=mass)
    if item == 'helium':
        inside = (x_UBT >= x0) and (x_UBT <= x1) and (y_UBT >= y0) and (y_UBT <= y1)
        if inside:
            h['Massdistr_at_UBT_helium_inside'].Fill(mass, weight)
        else:
            h['Massdistr_at_UBT_helium_outside'].Fill(mass, weight)

    # --- xyz accumulation (for surviving_xyzplots) ---
    if part_vtx is not None:
        key = f"{item}_pos"   # e.g. 'LiSc_pos', 'helium_pos', etc.
        g = _xyz_ensure(key)
        rx, ry, rz = part_vtx.X(), part_vtx.Y(), part_vtx.Z()
        g['pts']['reco']['zx'][0].append(rz); g['pts']['reco']['zx'][1].append(rx)
        g['pts']['reco']['zy'][0].append(rz); g['pts']['reco']['zy'][1].append(ry)
        g['pts']['reco']['xy'][0].append(rx); g['pts']['reco']['xy'][1].append(ry)

        try:
            IPx = event.MCTrack[1].GetStartX() #IP = interaction point of IS (should be the same as event.MCTrack[0] but IDK
            IPy = event.MCTrack[1].GetStartY()
            IPz = event.MCTrack[1].GetStartZ()
            g['pts']['ip']['zx'][0].append(IPz); g['pts']['ip']['zx'][1].append(IPx)
            g['pts']['ip']['zy'][0].append(IPz); g['pts']['ip']['zy'][1].append(IPy)
            g['pts']['ip']['xy'][0].append(IPx); g['pts']['ip']['xy'][1].append(IPy)
        except Exception:
            g['has_ip'] = False


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
# Confusion matrix: PID_CONFUSION[true_type][reco_type] = probability
# true_type and reco_type: 'e', 'mu', 'hadron'
# "x" / "->"" is predicted, "y" / "|^"" is true
# Data from Matei presented at the 35th collaboration meeting. 10 000 e, mu, nu were simulated at detector 10-150 GeV

PID_CONFUSION = {
    'e':      {'e': 0.9978, 'mu': 0.0, 'hadron': 0.0022},
    'mu':     {'e': 0.0001, 'mu': 0.9966, 'hadron': 0.0033},
    'hadron': {'e': 0.0061, 'mu': 0.0062, 'hadron': 0.9877},
}


def pid_apply_confusion(event, track_index, true_type):
    """
    Given the true particle type, return the reconstructed type
    according to the confusion matrix. Deterministic per (event, track).
    true_type: 'e', 'mu', or 'hadron'
    returns:   'e', 'mu', or 'hadron'
    """

    if PID_NO_CONFUSION:
        return true_type

    evt_time_ns = int(event.ShipEventHeader.GetEventTime())
    seed = (evt_time_ns * 131071 + track_index + 999983) & 0x7FFFFFFF
    ROOT.gRandom.SetSeed(seed)
    rnd = ROOT.gRandom.Rndm()  # uniform in [0, 1)

    probs = PID_CONFUSION[true_type]
    cumulative = 0.0
    for reco_type, prob in probs.items():
        cumulative += prob
        if rnd < cumulative:
            return reco_type
    return 'hadron'  # fallback, should never be reached if rows sum to 1

def pid_decision(event,candidate):  
    """
    Interim solution for PID check:Uses track truth info and aoplies PID_Confusion.

    pid_code:   0     = hadronic,
                1     = dileptonic (any leptons),
                1.1   = dileptonic ee,
                1.2   = dileptonic μμ,
                1.3   = dileptonic eμ,         
                2     = semileptonic(any lepton),
                2.1   = semileptonic containing an e,
                2.2   = semileptonic containing a μ,
                3     = at least one track has unknown PID, #never used since truth but added for historical reasons
                4     = fewer than two PID tracks available #two track candidates

    """

    if(len(event.Pid)<2):
        print("Pid is less than 2 particles!") #sanity check
        return 4

    d1_mc=event.MCTrack[event.fitTrack2MC[candidate.GetDaughter(0)]]
    d1_pdg=d1_mc.GetPdgCode()
    d2_mc=event.MCTrack[event.fitTrack2MC[candidate.GetDaughter(1)]]
    d2_pdg=d2_mc.GetPdgCode()

    LEPTON_PDGS = {11, 13}      # 11 = electron, 13 = muon

    d1_is_lepton = abs(d1_pdg) in LEPTON_PDGS
    d2_is_lepton = abs(d2_pdg) in LEPTON_PDGS

    d1_is_mu= (abs(d1_pdg)==13)
    d2_is_mu= (abs(d2_pdg)==13)

    d1_is_e= (abs(d1_pdg)==11)
    d2_is_e= (abs(d2_pdg)==11)    

    # Track 1
    d1_reco      = pid_apply_confusion(event, candidate.GetDaughter(0), true_type(d1_pdg))
    #d1_reco = true_type(d1_pdg) #for now, no confusion matrix applied since we want to see the effect of the PID efficiency alone. Will be added in the future
    d1_is_lepton = d1_reco in ('e', 'mu')
    d1_is_mu     = (d1_reco == 'mu')
    d1_is_e      = (d1_reco == 'e')

    # Track 2
    d2_reco      = pid_apply_confusion(event, candidate.GetDaughter(1), true_type(d2_pdg))
    #d2_reco = true_type(d2_pdg) #for now, no confusion matrix applied since we want to see the effect of the PID efficiency alone. Will be added in the future
    d2_is_lepton = d2_reco in ('e', 'mu')
    d2_is_mu     = (d2_reco == 'mu')
    d2_is_e      = (d2_reco == 'e')


    if d1_is_lepton and d2_is_lepton:                       # dileptonic final state
        if d1_is_e and d2_is_e:                             # ee
            return 1.1                          
        if d1_is_mu and d2_is_mu:                           # mu mu 
            return 1.2
        if (d1_is_e and d2_is_mu) or (d1_is_mu and d2_is_e):    # mu e / e mu 
            return 1.3
        return 1

    if d1_is_lepton or d2_is_lepton:                        # semileptonic
        if (d1_is_e or d2_is_e):
            return 2.1
        if (d1_is_mu or d2_is_mu):
            return 2.2
        return 2

    return 0 


def mother_pdg_from_fittrack(event, fit_track_index):
    """Return PDG code of the mother of a fit track's MC match; 0 if unavailable."""
    try:
        mc_index = event.fitTrack2MC[fit_track_index]
    except Exception:
        return 0
    if mc_index < 0:
        return 0
    mc_track = event.MCTrack[mc_index]
    mid = mc_track.GetMotherId()
    if mid < 0:
        return 0
    return event.MCTrack[mid].GetPdgCode()

def UBT_decision(event): # not yet used 
    """Implementation of UBT veto. Simple MC check; no efficiency, no mom. check """
    nHits = 0
    for ahit in event.UpstreamTaggerPoint:
        nHits+=1

    if nHits:
        veto=True  
    else:
        veto=False 
    return veto, nHits

#### functions to call before vetos can be implemented 
# module globals
_genfit_field_ready = False
_field_maker = None
_genfit_bfield = None

def initial_AnalysisContext(ShipGeo):
    global _genfit_field_ready, _field_maker, _genfit_bfield
    if _genfit_field_ready:
        return
    ShipGeo.Bfield.fieldMap = os.path.join(fairship, "files/MainSpectrometerField.root")
    #ShipGeo.Bfield.fieldMap = "files/TRY_2025.root"
    # addVMCFields resolves relative file paths via $VMCWORKDIR; point it to FairShip so
    # 'files/MainSpectrometerField.root' resolves correctly instead of landing in FairRoot examples.
    os.environ['VMCWORKDIR'] = fairship
    # MuonDIS geofiles don't include EmuMagnet; addVMCFields accesses it without a hasattr guard.
    if not hasattr(ShipGeo, 'EmuMagnet'):
        ShipGeo.EmuMagnet = AttrDict({'MagneticField': False})
    _field_maker = geomGeant4.addVMCFields(ShipGeo, '', True, withVirtualMC=False)

    geoMat = ROOT.genfit.TGeoMaterialInterface()
    ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
    _genfit_bfield = ROOT.genfit.FairShipFields()
    _genfit_bfield.setField(_field_maker.getGlobalField())
    field_mgr = ROOT.genfit.FieldManager.getInstance()
    field_mgr.init(_genfit_bfield)

    _genfit_field_ready = True

def dileptonic_ip_tresh(part_vtx):
    z = part_vtx.Z()
    z_entrance = veto_geo.z0
    z_exit = - z_entrance
    #dileptonic_ip_tresh(z) = m * z + b = -240 / (z_entrance - z_exit) * z + b = -240 / (z_entrance - z_exit) * z + 10 + 240 / (z_entrance - z_exit) * z_entrance
    #dileptonic_ip_tresh(z_entrance) = 10 = m * z_entrance + b
    #dileptonic_ip_tresh(z_exit) = 250 = m * z_exit + b
    #dileptonic_ip_tresh(z_entrance) - dileptonic_ip_tresh(z_exit) = -240 = m * (z_entrance - z_exit) => -240 / (z_entrance - z_exit) = m
    #dileptonic_ip_tresh(z_entrance) =  -240 / (z_entrance - z_exit) * z_entrance + b = 10 => 10 + 240 / (z_entrance - z_exit) * z_entrance = b
    return -240 / (z_entrance - z_exit) * z + 10 + 240 / (z_entrance - z_exit) * z_entrance

def slice_function(z):
    '''Returns the y border of the slice at a given z. The slice is defined as |y|<(z-z_{entrance})•25 cm/(z_{exit}-z_{entrance})+15 cm.'''
    z_entrance =  veto_geo.z0
    z_exit = - z_entrance
    y_border =(z-z_entrance)*25/(z_exit-z_entrance)+15
    return y_border

def slice_sides_SBT(event):
    '''Returns true, if the truth IS vertex is in the slice |y|<(z-z_{entrance})•25 cm/(z_{exit}-z_{entrance})+15 cm. Than it is not included.'''
    z = event.MCTrack[0].GetStartZ()
    y = event.MCTrack[0].GetStartY() 
    y_border = slice_function(z)
    if abs(y) <= y_border:
        return True
    else:
        return False 

def true_type(pdg):
    if abs(pdg) == 11: return 'e'
    if abs(pdg) == 13: return 'mu'
    return 'hadron'

### veto functions 
SBT_EFFICIENCY = 0.99 
#UBTefficiency = 0.9   # Upstream background tagger
random = ROOT.TRandom(13)

def sbt_cell_fired(event, detID, efficiency=SBT_EFFICIENCY):
    """Deterministic live/dead state for SBT cells per event."""
    evt_time_ns = int(event.ShipEventHeader.GetEventTime())
    seed = (evt_time_ns * 131071 + detID) & 0x7FFFFFFF
    ROOT.gRandom.SetSeed(seed)
    return ROOT.gRandom.Rndm() < efficiency
   

def z_SBTcell(Zlayer): # currently not used
    """Return z position of the end of the SBT cell layer in regards to Vessel entrance."""
    thickness_cell1 = 80 # cm
    thickness_cells = 82 #cm
    z_pos =  thickness_cell1 + (Zlayer - 1) * thickness_cells 
    return z_pos

def extrapolateTrackToSBT(event, fitIndex, tol_cm=320.0, back_dist_m=60, n_steps=300,Digi_SBTHits=None):
        """
        
        Extrapolate a fitted GenFit track backwards onto SBT.
        Uniformly sample the trajectory in n_steps, stop when we first enter 
        any LiSc volume, and then match to the nearest digi hits within the tolerance (tol_cm).

        Returns:
          best_hits, xs, ys, zs
        
        """
        """
        if not self.fM:
            geoMat =  ROOT.genfit.TGeoMaterialInterface()
            ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
            bfield = ROOT.genfit.FairShipFields()
            bfield.setField(self.fieldMaker.getGlobalField())
            self.fM = ROOT.genfit.FieldManager.getInstance()
            self.fM.init(bfield)

        """
        track = event.FitTracks[fitIndex]
        fst = track.getFitStatus() 
        if not (fst.isFitConverged() and fst.getNdf() > 0): #exclude tracks which do not converge
            return [], [], [], [] 
        # get fitted state & build the RK rep
        fstate = track.getFittedState(0)
        pos0   = fstate.getPos()
        mom0   = fstate.getMom()
        rep    = ROOT.genfit.RKTrackRep(fstate.getPDG())
        state  = ROOT.genfit.StateOnPlane(rep)
        rep.setPosMom(state, pos0, mom0)

        nav = ROOT.gGeoManager.GetCurrentNavigator()

        dx, dy, dz = (-mom0.Unit()).X(), (-mom0.Unit()).Y(), (-mom0.Unit()).Z()
        nav.InitTrack(pos0.X(), pos0.Y(), pos0.Z(), dx, dy, dz)

        back_cm = back_dist_m * 100
        ds      = -back_cm / float(n_steps)
        xs = []; ys = []; zs = []
        predPos = None

        for i in range(n_steps+1):
            p = state.getPos()
            xs.append(p.X()); ys.append(p.Y()); zs.append(p.Z())

            node = nav.FindNode(p.X(), p.Y(), p.Z())
            if node and node.GetName().startswith(("LiSc", "VetoInnerWall", "VetoOuterWall","VetoVerticalRib","VetoLongitRib")):
                predPos = p
                predMom = state.getMom()
                break

            target = p + state.getMom().Unit()*ds
            try:
                #rep.extrapolateToPoint(state, target, 0)
                rep.extrapolateToPoint(state, target, False)
            except Exception as e:
                print(f"Exception at step {i}: {e}")
                break

        # if we never hit LiSc, still push to the first boundary
        
        if predPos is None:
            # reset & do one boundary‐stop propagate
            rep.setPosMom(state, pos0, mom0)
            full_target = pos0 + mom0.Unit()*ds*n_steps
            
            #rep.extrapolateToPoint(state, full_target, 1)
            rep.extrapolateToPoint(state, full_target, True)  
            predPos = state.getPos()
            predMom = state.getMom()
            return [], xs, ys, zs

        # match the nearest SBT hits
        
        hits_in_tol = []      # will hold tuples of (hit) within the 320 cm of the track 
        
        if Digi_SBTHits==None:
            Digi_SBTHits=event.Digi_SBTHits
        

        for hit in Digi_SBTHits:
            
            if not sbt_cell_fired(event,hit.GetDetectorID()):
                continue
            
            d = (hit.GetXYZ() - predPos).Mag()
            if d < tol_cm:

                hits_in_tol.append(hit)

        return hits_in_tol, xs, ys, zs

def Digi_Hit_beforeCutSBT(hit, cutSBTDIS): #maybe this could be changed by using hit.GetXYZ()
    '''Returns True if there is any SBT hit before cutSBTDIS cm from vessel entrance. Then the event cannot be vetoed.'''
    #detectorID = hit.GetDetectorID()
    #detIDstr = str(detectorID)
    #Zlayer = int(detIDstr[2:4])
    #if z_SBTcell(Zlayer) <= cutSBTDIS:
     #   print("z position of SBT cell:",z_SBTcell(Zlayer))
      #  return True 
    if (hit.GetXYZ().Z() - veto_geo.z0) <= cutSBTDIS:
        return True # cannot be vetoed if there is no SBT between [0,cutSBTDIS] cm
    else:
        return False



# ---------- analysis ----------#
def main_analysis(event, sgeo, ShipGeo, rescale_fn=None, eventNr=None, counts=None, finalstate=None):
    # ---------- weight calculation ----------#
    cat = dis_region_basename(event,sgeo)
    if cat=='DecayVacuum':
        corrected_rhoL= helium_rhoL(event)
        rhoL=corrected_rhoL
    else:
        rhoL=event.MCTrack[2].GetWeight()
    weight = define_muon_weight(event,SHiP_running=15, w_DIS=rhoL)
    
    # ---------- material of scattering point ----------#
    origin_node = sgeo.FindNode(event.MCTrack[0].GetStartX(),event.MCTrack[0].GetStartY(),event.MCTrack[0].GetStartZ()).GetName()
    baseName = origin_node.split("_")[0]
    if baseName[:4] == 'LiSc': baseName = 'LiSc'
    if baseName == 'VetoVerticalRib': baseName = 'VetoLongitRib'

    region_label = region_label_from_basename(baseName)

    # ---------- check vertex efficiency ----------#
    compute_vertexing_efficiency(event, sgeo, h, vtx_eff_counts, ShipGeo, veto_geo)

    # ---------- reject non-SBT IS early  ----------# 
    if onlySBTDIS and baseName not in sbt_region_names:
        # skip event entirely if IS did not occur in SBT
        return
    if cutSBTDIS and baseName in sbt_region_names and (event.MCTrack[0].GetStartZ()-veto_geo.z0) < cutSBTDIS:
        # skip event entirely if IS occurred in the SBT but before cutSBTDIS cm from vessel entrance
        return
    if cutSBTDISiny and baseName in sbt_region_names and slice_sides_SBT(event):
        # skip event entirely if IS occurred in the "slice" of the SBT
        return

    # ---------- calculate pure PID efficiency calculations - before any cuts (except for those events taking place we do not include when running with cutSBTDIS or so) ----------#  
    if len(event.Particles) > 0:
        in_he = (baseName == 'DecayVacuum')
        in_sbt = baseName in sbt_region_names
        if in_he or in_sbt:
            region_key = 'He' if in_he else 'SBT'
            #for PID eff
            pid_eff_event_counts[region_key] += weight
            pid_eff_counts['all candidates'][region_key] += weight
            ee_any,mumu_any,emu_any,ex_any,mux_any,ll_any,lx_any = False,False,False,False,False,False,False
            #for cut eff
            gd_any,reco_cand_any,fiducial_any,doca_any,ip10_any,ip250_any,ipz_any,sbt_veto_any,ubt_veto,mass_any,ip_partial_any = False,False,False,False,False,False,False,False,False,False,False
            any_recocand = True # asked for len(event.Particles)>0 already 
            sbt_veto_any_per_thr = {thr: False for thr in veto_thresholds} 
            


            for part in event.Particles: ## all particle candidates that are found in reconstruction -> everything with two tracks
                status1 = event.FitTracks[part.GetDaughter(0)].getFitStatus()
                status2 = event.FitTracks[part.GetDaughter(1)].getFitStatus()
                rounded_status1 = int(round(status1.getNdf()))
                rounded_status2 = int(round(status2.getNdf()))
                selected_vtx = ROOT.TVector3()
                part.GetVertex(selected_vtx)
                if rounded_status1 > 25 and rounded_status2 > 25 and status1.getChi2()/status1.getNdf() < 5 and status2.getChi2()/status2.getNdf() < 5 and event.FitTracks[part.GetDaughter(0)].getFittedState().getMom().Mag() > 1 and event.FitTracks[part.GetDaughter(1)].getFittedState().getMom().Mag() > 1: # has a reco and Good Daughters 
                    #for PID eff 
                    #fill a histogram with the energyspectrum of the events surviving the cuts 
                    selected_mom = ROOT.TLorentzVector()
                    part.Momentum(selected_mom)
                    Energy = abs(selected_mom.E())
                    h[f'Energy_distribution'].Fill(Energy, weight)
                    pid_code = pid_decision(event, candidate=part) # calculate the PID code based on truth info

                    if pid_code == 1.1 or int(pid_code) == 3:
                        ee_any = True
                    if pid_code == 1.2 or int(pid_code) == 3:
                        mumu_any = True
                        if pid_code == 1.2:
                            d1 = part.GetDaughter(0) #link to the first fitted track of the candidate
                            d2 = part.GetDaughter(1)
                            pdg1 = mother_pdg_from_fittrack(event, d1) # Return PDG code of the mother of a fit track's MC match
                            pdg2 = mother_pdg_from_fittrack(event, d2)
                            key = tuple(sorted((pdg1, pdg2))) # sorted such that (1,2)=(2,1)
                            mumu_origin_counts[key] += weight
                    if pid_code == 1.3 or int(pid_code) == 3:
                        emu_any = True

                    if pid_code == 2.1 or int(pid_code) == 3:
                        ex_any = True
                    if pid_code == 2.2 or int(pid_code) == 3:
                        mux_any = True
                    if int(pid_code) == 1 or int(pid_code) == 3:
                        ll_any = True
                    if int(pid_code) == 2 or int(pid_code) == 3:
                        lx_any = True
                    
                    # for cut eff
                    if rounded_status1 > 25 and rounded_status2 > 25 and status1.getChi2()/status1.getNdf() < 5 and status2.getChi2()/status2.getNdf() < 5 and event.FitTracks[part.GetDaughter(0)].getFittedState().getMom().Mag() > 1 and event.FitTracks[part.GetDaughter(1)].getFittedState().getMom().Mag() > 1: # Good Daughters 
                        gd_any = True
                        #SBT Extrapolation veto
                        xs, ys, zs, bestHits = [],[],[],[]
                        track_index_first,track_index_last = part.GetDaughter(0),part.GetDaughter(1)
                        for tr in [track_index_first,track_index_last]:
                            bestHit,xs_, ys_, zs_= extrapolateTrackToSBT(event,tr)
                            xs.append(xs_)
                            ys.append(ys_)
                            zs.append(zs_)
                            if len(bestHit):
                                bestHits.extend(bestHit)       
                        for hit in bestHits:
                            ELoss    = hit.GetEloss()
                            if SBTVeto is not None and (ELoss>=SBTVeto*0.001):
                                sbt_veto_any = True
                                for thr in veto_thresholds:
                                    if ELoss * 1000 >= thr:
                                        sbt_veto_any_per_thr[thr] = True
                    if (len(event.Particles) == 1):
                        reco_cand_any = True
                    if (dist2InnerWall(selected_vtx,sgeo) > dist2iWall and dist2Entrance(selected_vtx) > 20 and is_in_fiducial(part, event, sgeo, ShipGeo)):
                        fiducial_any = True
                    if (part.GetDoca() < 1):
                        doca_any = True
                    if (impact_parameter(selected_vtx, selected_mom, ShipGeo) < 10):
                        ip10_any = True
                    if (impact_parameter(selected_vtx, selected_mom, ShipGeo) < 250):
                        ip250_any = True
                    if (impact_parameter(selected_vtx, selected_mom, ShipGeo) < dileptonic_ip_tresh(selected_vtx)):
                        ipz_any = True
                    if partial_IP_cut and (impact_parameter(selected_vtx, selected_mom, ShipGeo) < partial_IP_cut):
                        ip_partial_any = True
                    if selected_mom.M() > 0.15: #
                        mass_any = True
                    
                # if (UBT_decision(event)):
                #    ubt_veto_any = True
                    nHits = 0
                    for ahit in event.UpstreamTaggerPoint:
                        nHits+=1
                    if nHits>1:
                        ubt_veto = True

            #eventwise counts for PID eff
            if ee_any:
                pid_eff_counts['ee'][region_key] += weight
            if mumu_any:
                pid_eff_counts['mu mu'][region_key] += weight
                #motherid = event.MCTrack[1].GetMotherId() 
                #print(motherid)
            if emu_any:
                pid_eff_counts['e mu'][region_key] += weight
            if ex_any:
                pid_eff_counts['eX'][region_key] += weight
            if mux_any:
                pid_eff_counts['mu X'][region_key] += weight
            if ll_any:
                pid_eff_counts['ll'][region_key] += weight
            if lx_any:
                pid_eff_counts['lx'][region_key] += weight

            #eventwise counts for cut eff

            cand_types = [
            (ee_any,   "ee"),
            (mumu_any, "mumu"),
            (emu_any,  "emu"),
            (ex_any,   "ex"),
            (mux_any,  "mux"),
            (ll_any,   "ll"),
            (lx_any,   "lx"),
            (True,     "all")]

            for cand_flag, cand_label in cand_types:
                if not cand_flag:
                    continue
        
                basic_cuts = False 
                if any_recocand:
                    cut_eff_counts['has a reco candidate'][region_key][cand_label] += weight
                if gd_any:
                    cut_eff_counts['good daughter'][region_key][cand_label] += weight
                if any_recocand and gd_any:
                    cut_eff_counts['has reco cand + good daughters'][region_key][cand_label] += weight
                if reco_cand_any:
                    cut_eff_counts['1reco cand'][region_key][cand_label] += weight
                if fiducial_any:
                    cut_eff_counts['fiducial'][region_key][cand_label] += weight
                if doca_any:
                    cut_eff_counts['DOCA'][region_key][cand_label] += weight
                if ip10_any:
                    cut_eff_counts['IP<10'][region_key][cand_label] += weight
                if ip250_any:
                    cut_eff_counts['IP<250'][region_key][cand_label] += weight
                if ip_partial_any:
                    _row = f'IP<{int(partial_IP_cut)}'
                    if _row in cut_eff_counts:
                        cut_eff_counts[_row][region_key][cand_label] += weight
                if ipz_any:
                    cut_eff_counts['IP<IP(z)'][region_key][cand_label] += weight
                if not sbt_veto_any and any_recocand and gd_any:
                    cut_eff_counts['SBT veto'][region_key][cand_label] += weight
                for thr in veto_thresholds:
                    if not sbt_veto_any_per_thr[thr]:
                        cut_eff_counts_per_veto[thr]['SBT veto'][region_key][cand_label] += weight
                if not ubt_veto:
                    cut_eff_counts['UBT Veto'][region_key][cand_label] += weight
                if mass_any: 
                    cut_eff_counts['mass > 0.15 GeV'][region_key][cand_label] += weight
                if any_recocand and gd_any and reco_cand_any and fiducial_any and doca_any:
                    basic_cuts=True 
                if basic_cuts and ip10_any:
                    cut_eff_counts['basic cuts + IP<10'][region_key][cand_label] += weight
                if basic_cuts and ip250_any:
                    cut_eff_counts['basic cuts + IP<250'][region_key][cand_label] += weight
                if basic_cuts and ipz_any:
                    cut_eff_counts['basic cuts + IP<IP(z)'][region_key][cand_label] += weight
                if basic_cuts and mass_any:
                    cut_eff_counts['basic cuts + mass > 0.15 GeV'][region_key][cand_label] += weight
                if basic_cuts and ip10_any and not sbt_veto_any:
                    cut_eff_counts['basic cuts + IP<10 + SBT veto'][region_key][cand_label] += weight
                if basic_cuts and ip250_any and not sbt_veto_any:
                    cut_eff_counts['basic cuts + IP<250 + SBT veto'][region_key][cand_label] += weight
                if basic_cuts and ipz_any and not sbt_veto_any:
                    cut_eff_counts['basic cuts + IP<IP(z) + SBT veto'][region_key][cand_label] += weight
                if basic_cuts and mass_any and not sbt_veto_any:
                    cut_eff_counts['basic cuts + mass > 0.15 GeV + SBT veto'][region_key][cand_label] += weight

    # ---------- where in X,Y does scatteirng in the SBT happen? ----------#    
    if baseName in sbt_region_names:
        #x_SBT =part_vtx.X()
        #y_SBT = part_vtx.Y() #reco values
        x_DIS = event.MCTrack[0].GetStartX()
        y_DIS = event.MCTrack[0].GetStartY()
        z_DIS = event.MCTrack[0].GetStartZ()
        h['x-y-IS-SBT'].Fill(x_DIS, y_DIS, weight)
        h['z-IS-SBT'].Fill(z_DIS, weight)

    # ---------- apply cuts ----------#  
    if region_label:
        update_selection_counts(region_label, 0, weight)  # DIS in region
        update_selection_rawcounts(region_label,0)
        selected_candidate = None
        selected_vtx = None
        selected_mom = None
        if len(event.Particles) > 0:
            update_selection_counts(region_label, 1, weight)  # has reco candidate
            update_selection_rawcounts(region_label,1)
            for part in event.Particles:
                # check quality cuts on candidates; use first candidate satisfying chain.
                # PID is applied later (after the IP cut) to match Anupama's cutflow ordering.
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
            #print("mass difference ", selected_mom.M()- part.GetMass()) these are essentially the same
            update_selection_counts(region_label, 2, weight)  # quality cuts (nDoF, chi2, p)
            update_selection_rawcounts(region_label,2)
            if len(event.Particles) == 1:
                update_selection_counts(region_label, 3, weight)  # exactly 1 reco candidate
                update_selection_rawcounts(region_label,3)
                #distance to inner wall vs z vertex
                z_vtx=event.MCTrack[0].GetStartZ()
                d2wall = dist2InnerWall(selected_vtx, sgeo)
                h['Dist2WallvsVtx_z'].Fill(z_vtx, d2wall)
                if dist2InnerWall(selected_vtx,sgeo) > dist2iWall and dist2Entrance(selected_vtx) > 20 and is_in_fiducial(selected_candidate, event, sgeo, ShipGeo):
                    update_selection_counts(region_label, 4, weight)  # fiducial
                    update_selection_rawcounts(region_label,4)
                    if selected_candidate.GetDoca() < 1:
                        update_selection_counts(region_label, 5, weight)  # DOCA
                        update_selection_rawcounts(region_label,5)
                        if ip_cut:
                        #  if impact_parameter(selected_vtx,selected_mom,ShipGeo) <= ip_cut:
                         #       update_selection_counts(region_label, 5, weight)  # IP
                            ip_satisfied = False
                            if channel == 'fullreco':
                                if impact_parameter(selected_vtx, selected_mom, ShipGeo) < ip_cut:
                                    update_selection_counts(region_label, 6, weight)  # IP
                                    update_selection_rawcounts(region_label,6)
                                    ip_satisfied = True
                            elif channel == 'partialreco':
                                if partial_IP_cut and impact_parameter(selected_vtx, selected_mom, ShipGeo) < partial_IP_cut:
                                    update_selection_counts(region_label, 6, weight)  # IP
                                    update_selection_rawcounts(region_label,6)
                                    ip_satisfied = True
                                elif impact_parameter(selected_vtx, selected_mom, ShipGeo) < dileptonic_ip_tresh(selected_vtx):
                                    ip_satisfied = True
                                    update_selection_counts(region_label, 6, weight)  # IP
                                    update_selection_rawcounts(region_label,6)

                            # PID cut, applied last (after IP) so the cutflow matches selection_steps order
                            pid_satisfied = False
                            if ip_satisfied:
                                pid_code = pid_decision(event, candidate=selected_candidate) #makes a decision with certain efficiency and confusion matrix
                                pid_leptonic = (int(pid_code) == 1 or int(pid_code) == 3)
                                pid_semileptonic = (int(pid_code) == 2 or int(pid_code) == 3)
                                if PID:
                                    if finalstate=='dileptonic':
                                        pid_satisfied = pid_leptonic
                                    elif finalstate=='semileptonic':
                                        pid_satisfied = pid_semileptonic
                                else:
                                    pid_satisfied = True  # PID not activated

                                if pid_satisfied:
                                    update_selection_counts(region_label, 7, weight)  # incl. final state PID
                                    update_selection_rawcounts(region_label,7)

                            #extrapolation veto
                            if pid_satisfied:
                                Energy = abs(selected_mom.E())
                                h[f'Energy_distribution2'].Fill(Energy, weight)
                                #SBT Extrapolation veto
                                xs, ys, zs, bestHits = [],[],[],[]
                                AdvSBT_TagVeto = False
                                track_index_first,track_index_last = selected_candidate.GetDaughter(0),selected_candidate.GetDaughter(1)
                                #print("Track indices of the two daughter tracks:", track_index_first, track_index_last)

                                for tr in [track_index_first,track_index_last]:
                                    #print('extrapolation print:',extrapolateTrackToSBT(event,tr))
                                    bestHit,xs_, ys_, zs_= extrapolateTrackToSBT(event,tr)
                                    xs.append(xs_)
                                    ys.append(ys_)
                                    zs.append(zs_)

                                    if len(bestHit):
                                        bestHits.extend(bestHit)

                                valid_hits = []
                                for hit in bestHits:
                                    #print ("Best SBT hit position:", hit.GetXYZ().X(), hit.GetXYZ().Y(), hit.GetXYZ().Z())
                                    #print("slice function return:", slice_function(hit.GetXYZ().Z()))
                                    if cutSBTDIS and Digi_Hit_beforeCutSBT(hit, cutSBTDIS): # for every hit we check whether its a hit before cutSBTDIS m -if so, then no veto can be applied
                                        continue 
                                    elif cutSBTDISiny and abs(hit.GetXYZ().Y()) <= slice_function(hit.GetXYZ().Z()): # the y of the digi hit is in the slice than veto cannot be applied 
                                        continue
                                    valid_hits.append(hit)
                                #main veto decision
                                for hit in valid_hits:
                                    ELoss    = hit.GetEloss()
                                    if SBTVeto is not None and ELoss>= SBTVeto*0.001:
                                        AdvSBT_TagVeto=True
                                #create plots for different extrapolation veto treshholds
                                if veto_thresholds and ip_satisfied:
                                    max_eloss_mev = max((hit.GetEloss() * 1000 for hit in valid_hits), default=0.0)
                                    for thr in veto_thresholds:
                                        veto_fired = (max_eloss_mev >= thr)
                                        if not veto_fired:
                                            counts_per_veto[thr][region_label][8] += weight

                                if SBTVeto is not None and AdvSBT_TagVeto:
                                    print(f"AdvSBT veto applied with threshold {SBTVeto} MeV")
                                    pass  # veto triggered, do not count
                                else:
                                    update_selection_counts(region_label, 8, weight)  # SBT Veto passed
                                    update_selection_rawcounts(region_label,8)
                                    if options.mass_cut and selected_mom.M()<= 0.15:
                                        pass  # veto triggered, do not count
                                    else:
                                        if options.mass_cut:
                                            update_selection_counts(region_label, 9, weight)  # potential mass cut passed
                                            update_selection_rawcounts(region_label,9)
                                        else:
                                            fill_SBT_plots(event, sgeo, ShipGeo, selected_vtx, selected_mom, weight, selected_mom.M(),baseName=baseName)
                                            d1 = selected_candidate.GetDaughter(0)
                                            d2 = selected_candidate.GetDaughter(1)
                                            pdg1 = mother_pdg_from_fittrack(event, d1)
                                            pdg2 = mother_pdg_from_fittrack(event, d2)
                                            mu_origin_counts.append((region_label, pdg1, pdg2, selected_mom.M()))
                                            mu_origin_counts.append((region_label, pdg1, pdg2, selected_mom.M()))
                                            if baseName in sbt_region_names:
                                                h['z-IS-SBT-after-cuts'].Fill(event.MCTrack[0].GetStartZ(), weight)
                                                # build a single multi‑line string instead of printing directly
                                                msg = []
                                                msg.append(f"DIS in SBT passing cuts - region label: {region_label}\n")
                                                msg.append(f"momentum: {selected_mom.Px()} {selected_mom.Py()} "
                                                        f"{selected_mom.Pz()} {selected_mom.E()}\n")
                                                msg.append(f"vertex position: {selected_vtx.X()} {selected_vtx.Y()} "
                                                        f"{selected_vtx.Z()}\n")
                                                msg.append(f"Mass: {selected_mom.M()}\n")
                                                #print("Best hits in SBT:", bestHits)   # optional
                                                for hit in bestHits:
                                                    msg.append("-- new SBT hit --\n")
                                                    msg.append(f"E-Loss (in SBT?): {hit.GetEloss()}\n")
                                                    msg.append(f"Where hit in SBT: {hit.GetXYZ().X()} "
                                                            f"{hit.GetXYZ().Y()} {hit.GetXYZ().Z()}\n")
                                                    msg.append("\n")
                                                # keep console output if you still want
                                                print(''.join(msg), end='')
                                                event_log_lines.append(''.join(msg))


                                    





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
    global _genfit_field_ready
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
                    if not _genfit_field_ready:
                        initial_AnalysisContext(ShipGeo)
                
                if options.testing_code and files > 1:
                    break
                
                print(files, jobDir)
                files += 1


                for eventNr, event in enumerate(tree):
                    try:
                        main_analysis(event, sgeo, ShipGeo, finalstate=finalstate)
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
    outdir = os.path.join(directory, tag + options_tag)
    os.makedirs(outdir, exist_ok=True)
    output_base = os.path.join(outdir, '')  # ergibt outdir + '/' als Präfix
    persist_xy_weight_sums(output_base)
    persist_selection_table(output_base)
    persist_pid_efficiencies(output_base)
    persist_cut_efficiencies(output_base)
    persist_mumu_origin_counts(output_base)
    persist_mu_origin_counts(output_base)
    persist_SBT_stats(output_base)
    persist_selection_rawtable(output_base)
    persist_vertexing_efficiency(output_base, h, vtx_eff_counts)
    keylist_colz.append('vtx_eff_z_ratio')
    keylist_colz.append('charge_misid_z')
    ut.writeHists(h, output_base + 'plots.root')
    print(f"Histograms saved to {output_base + 'plots.root'}")
    # Generate xyz scatter plots
    if xyz_groups:
        print("Generating xyz position plots...")
        xyzplots.generate_plots_from_data(
        xyz_groups,
        outdir_name=f"plots_{options_tag}",
        mass_cut=options.mass_cut,
        sbt_veto=SBTVeto,
        test=options.testing_code
    )
    # Generate Veto threshold plots 
    if veto_thresholds:
        for thr in veto_thresholds:
            for reg in region_labels:
                for step in range(8):
                    counts_per_veto[thr][reg][step] = counts[reg][step]
        make_TestVetos_plots(output_base)
    for key in keylist_hist:
        h[key].SetOption('HIST')
    for key in keylist_colz:
        h[key].SetOption('COLZ')
    ut.writeHists(h, output_base + 'plots.root')
    print('done')

# ---------- run analysis ----------#

Main_function()

