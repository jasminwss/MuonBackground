import rootUtils as ut
from rootpyPickler import Unpickler
import ROOT, os
import geomGeant4
from tabulate import tabulate
import yaml
from ShipGeoConfig import AttrDict
from array import array

Mp = 0.938272  # proton mass [GeV]
test = False


def main_analysis(event):
    inMuon = event.InMuon.At(0) # first / only muon of event 
    #print(type(inMuon))          # Python wrapper type
    #print(inMuon.ClassName())    # actual ROOT class, e.g. TVectorD
    #print(inMuon.GetNrows())     # if TVectorD, how many entries it holds
    #for i in range(inMuon.GetNrows()):
    #    print(i, inMuon[i])
    # this is from         mu = array("d",[pid,px,py,pz,E,x,y,z,w,isProton,xsec,time_muon,args.nDIS,nmuons,],
    pdg_in = int(round(inMuon[0])) # PDG code of this first muon (index 0 stores PDG)
    px, py, pz, E = inMuon[1], inMuon[2], inMuon[3], inMuon[4]
    p_in = ROOT.TVector3(px, py, pz) #3d momentum vector of incoming muon

    # scattered muon: outgoing particle with the same flavour as the incoming muon
    # (take the highest-energy match in case of extra mu pairs from trident-like processes)
    scattered = None
    for part in event.DISParticles: # for all particles in the ones from DIS
        pdg_out = int(round(part[0])) # muon or antimuon 
        if abs(pdg_out) == abs(pdg_in): # aka a muon
            #  take muon with higher energy
            if scattered is None or part[4] > scattered[4]: # keeps the candidate or takes the one with higher energy 
                scattered = part # update best candidate 

    if scattered is None:
        return  # no outgoing muon found -> abort mission

    px2, py2, pz2, E2 = scattered[1], scattered[2], scattered[3], scattered[4]
    p_out = ROOT.TVector3(px2, py2, pz2) # 3d momentum of best candidate outgoing muon

    costheta = p_in.Dot(p_out) / (p_in.Mag() * p_out.Mag())
    Q2 = 2 * E * E2 * (1 - costheta)
    nu = E - E2
    xBj = Q2 / (2 * Mp * nu) if nu > 0 else -1

    W2 = Mp**2 + 2 * Mp * nu - Q2
    W = W2**0.5 if W2 > 0 else -1

    return Q2, xBj, W


def Main_function():
    global h
    h = {}

    results = []  # (Q2, xBj, W) tuples collected before we know the axis maxima

    files = 0
    f = None
    fgeo = None
    sgeo = None
    exception_issues = {}
    
    # Create list of paths to process
    paths_to_process = ['/eos/experiment/ship/simulation/bkg/MuonDIS_2024helium/8070735/SBT', '/eos/experiment/ship/simulation/bkg/MuonDIS_2024helium/8070735/Tr']

    # iterate over validated paths, skip ones that don't exist
    global _genfit_field_ready
    for current_path in paths_to_process:
        if not os.path.isdir(current_path):
            print(f"Warning: path does not exist or is not a directory: '{current_path}'. Skipping.")
            continue
        print(f"Processing path: {current_path}")
        
        # Process each job directory in current path
        for jobDir in os.listdir(current_path):
            if test and files > 10 :
                break 
            try:
                inputFile = f'{current_path}/{jobDir}/muonDis.root'
                
                f = ROOT.TFile.Open(inputFile)
                tree = f.DIS
                
                print(files, jobDir)
                files += 1

                for eventNr, event in enumerate(tree):
                    try:
                        result = main_analysis(event)
                        if result is not None:
                            results.append(result)

                    except Exception as e:
                        print(f'Except called for (Reason :{e})')
                        exception_issues[jobDir] = e
                        continue

                
                f.Close()
                
            except Exception as e:
                if f:
                    f.Close()
                print(f'Except called for (Reason :{e})')
                exception_issues[jobDir] = e
                continue

    Q2_max = max(r[0] for r in results)
    xBj_max = max(r[1] for r in results)
    W_max = max(r[2] for r in results)

    ut.bookHist(h, 'Q2', 'Q^{2} of muon DIS interaction; Q^{2} [GeV^{2}]; Events', 100, 0, Q2_max)
    ut.bookHist(h, 'xBjorken', 'Bjorken x of muon DIS interaction; x_{Bj}; Events', 100, 0, xBj_max)
    ut.bookHist(h, 'Q2_vs_xBj', 'Q^{2} vs x_{Bj}; x_{Bj}; Q^{2} [GeV^{2}]', 100, 0, xBj_max, 100, 0, Q2_max)
    ut.bookHist(h, 'W', 'Invariant mass of hadronic system; W [GeV]; Events', 100, 0, W_max)

    for Q2, xBj, W in results:
        h['Q2'].Fill(Q2)
        h['xBjorken'].Fill(xBj)
        h['Q2_vs_xBj'].Fill(xBj, Q2)
        h['W'].Fill(W)

    ut.writeHists(h, '/afs/cern.ch/work/j/jaweiss/private/MuonBackground/generator_analysis/gen_max_plots.root')
    print('done')

# ---------- run analysis ----------#

Main_function()

