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

region_labels = ['helium','LiSc','inner_wall','outer_wall','rib','UBT']
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
if options.mass_cut:
    selection_steps.append('mass > 0.15 GeV + SBT veto')
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
     directory = '/afs/cern.ch/work/j/jaweiss/private/'
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


keylist_hist.append('Massdistr_at_UBT_helium_inside')
keylist_hist.append('Massdistr_at_UBT_helium_outside')
keylist_hist.append('Energy_distribution')
keylist_hist.append('Energy_distribution2')
keylist_colz.append('Dist2WallvsVtx_z')
keylist_colz.append('x-y-IS-SBT')


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
    record_xy_weight(f'xy_at_UBT_{item}', x_UBT, y_UBT, weight)
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
    if ana_part != 'muon' or not xy_weight_sums:
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
    ShipGeo.Bfield.fieldMap = "files/MainSpectrometerField.root"
    #ShipGeo.Bfield.fieldMap = "files/TRY_2025.root"
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
   

#def pid_detector_fired(event, track_index, efficiency=PID_EFFICIENCY):
 #   """Deterministic live/dead state for PID detector per track per event."""
  #  evt_time_ns = int(event.ShipEventHeader.GetEventTime())
   # seed = (evt_time_ns * 131071 + track_index) & 0x7FFFFFFF
    #ROOT.gRandom.SetSeed(seed)
    #return ROOT.gRandom.Rndm() < efficiency #returns true or false

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
                rep.extrapolateToPoint(state, target, False)
            except Exception as e:
                print(f"Exception at step {i}: {e}")
                break

        # if we never hit LiSc, still push to the first boundary
        
        if predPos is None:
            # reset & do one boundary‐stop propagate
            rep.setPosMom(state, pos0, mom0)
            full_target = pos0 + mom0.Unit()*ds*n_steps
            
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
        h['x-y-IS-SBT'].Fill(x_DIS, y_DIS, weight)

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
                pid_code = pid_decision(event, candidate=part) #makes a decision with ceratain efficiency and confusion matrix
                pid_leptonic = (int(pid_code) == 1 or int(pid_code) == 3)
                pid_semileptonic = (int(pid_code) == 2 or int(pid_code) == 3)
            
                # Check if PID requirements are satisfied (or PID is not activated)
                pid_satisfied = False
                if PID:
                    if finalstate=='dileptonic':
                        pid_satisfied = pid_leptonic
                    elif finalstate=='semileptonic':
                        pid_satisfied = pid_semileptonic
                else:
                    pid_satisfied = True  # PID not activated
            
                # Vetos (only apply if PID requirements are met or PID is not activated)
                if pid_satisfied:
                    #print("Candidate passed PID requirements")
                    update_selection_counts(region_label, 2, weight)
                    update_selection_rawcounts(region_label,2)
                    # check quality and fiducial cuts on candidates; use first candidate satisfying chain

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
            update_selection_counts(region_label, 3, weight)  # quality cuts
            update_selection_rawcounts(region_label,3)
            if len(event.Particles) == 1:
                update_selection_counts(region_label, 4, weight)  # exactly 1 reco candidate
                update_selection_rawcounts(region_label,4)
                #distance to inner wall vs z vertex
                z_vtx=event.MCTrack[0].GetStartZ()
                d2wall = dist2InnerWall(selected_vtx, sgeo)
                h['Dist2WallvsVtx_z'].Fill(z_vtx, d2wall)
                if dist2InnerWall(selected_vtx,sgeo) > dist2iWall and dist2Entrance(selected_vtx) > 20 and is_in_fiducial(selected_candidate, event, sgeo, ShipGeo):
                    update_selection_counts(region_label, 5, weight)  # fiducial
                    update_selection_rawcounts(region_label,5)
                    if selected_candidate.GetDoca() < 1:
                        update_selection_counts(region_label, 6, weight)  # DOCA
                        update_selection_rawcounts(region_label,6)
                        if ip_cut:
                        #  if impact_parameter(selected_vtx,selected_mom,ShipGeo) <= ip_cut:
                         #       update_selection_counts(region_label, 6, weight)  # IP
                            ip_satisfied = False
                            if channel == 'fullreco':
                                if impact_parameter(selected_vtx, selected_mom, ShipGeo) < ip_cut:
                                    update_selection_counts(region_label, 7, weight)  # IP
                                    update_selection_rawcounts(region_label,7)
                                    ip_satisfied = True
                            elif channel == 'partialreco': 
                                if partial_IP_cut and impact_parameter(selected_vtx, selected_mom, ShipGeo) < partial_IP_cut:
                                    update_selection_counts(region_label, 7, weight)  # IP
                                    update_selection_rawcounts(region_label,7)
                                    ip_satisfied = True
                                elif impact_parameter(selected_vtx, selected_mom, ShipGeo) < dileptonic_ip_tresh(selected_vtx):
                                    ip_satisfied = True
                                    update_selection_counts(region_label, 7, weight)  # IP
                                    update_selection_rawcounts(region_label,7)
                                    
                            #PID
                            if ip_satisfied:
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


                                       # if selected_candidate.GetMass() <= 0.15 and mass_cut:
                                        #    pass  # veto triggered, do not count
                                       # else:
                                        #    
                                           


                                    





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
                
                if options.testing_code and files > 5:
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
    output_base = directory+tag+options_tag #+ana_part
    persist_xy_weight_sums(output_base)
    persist_selection_table(output_base)
    persist_pid_efficiencies(output_base)
    persist_cut_efficiencies(output_base)
    persist_mumu_origin_counts(output_base)
    persist_mu_origin_counts(output_base)
    persist_SBT_stats(output_base)
    persist_selection_rawtable(output_base)
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
    ut.writeHists(h, directory+tag+options_tag+'plots.root')
    print('done')

# ---------- run analysis ----------#

Main_function()

