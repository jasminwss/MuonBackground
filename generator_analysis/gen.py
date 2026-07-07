import rootUtils as ut
from rootpyPickler import Unpickler
import ROOT, os
import geomGeant4
from tabulate import tabulate
import yaml
from ShipGeoConfig import AttrDict
from array import array

Mp = 0.938272  # proton mass [GeV]
test = True


def main_analysis(event):
    inMuon = event.InMuon.At(0)
    pdg_in = int(round(inMuon[0]))
    px, py, pz, E = inMuon[1], inMuon[2], inMuon[3], inMuon[4]
    p_in = ROOT.TVector3(px, py, pz)

    # scattered muon: outgoing particle with the same flavour as the incoming muon
    # (take the highest-energy match in case of extra mu pairs from trident-like processes)
    scattered = None
    for part in event.DISParticles:
        pdg_out = int(round(part[0]))
        if abs(pdg_out) == abs(pdg_in):
            if scattered is None or part[4] > scattered[4]:
                scattered = part

    if scattered is None:
        return  # no outgoing muon found

    px2, py2, pz2, E2 = scattered[1], scattered[2], scattered[3], scattered[4]
    p_out = ROOT.TVector3(px2, py2, pz2)

    costheta = p_in.Dot(p_out) / (p_in.Mag() * p_out.Mag())
    Q2 = 2 * E * E2 * (1 - costheta)
    nu = E - E2
    xBj = Q2 / (2 * Mp * nu) if nu > 0 else -1

    h['Q2'].Fill(Q2)
    h['xBjorken'].Fill(xBj)
    h['Q2_vs_xBj'].Fill(xBj, Q2)


def Main_function():
    global h
    h = {}
    ut.bookHist(h, 'Q2', 'Q^{2} of muon DIS interaction; Q^{2} [GeV^{2}]; Events', 100, 0, 5)
    ut.bookHist(h, 'xBjorken', 'Bjorken x of muon DIS interaction; x_{Bj}; Events', 100, 0, 0.2)
    ut.bookHist(h, 'Q2_vs_xBj', 'Q^{2} vs x_{Bj}; x_{Bj}; Q^{2} [GeV^{2}]', 100, 0, 0.2, 100, 0, 5)

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
                        main_analysis(event)

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
    ut.writeHists(h, '/afs/cern.ch/work/j/jaweiss/private/MuonBackground/generator_analysis/gen_plots.root')
    print('done')

# ---------- run analysis ----------#

Main_function()

