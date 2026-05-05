import ROOT

def check_charge_misid(event, candidate):
    """
    Check if either daughter track has a wrong reconstructed charge.
    Returns True if a mis-ID is found.
    """
    misid = False
    for daughter_idx in [candidate.GetDaughter(0), candidate.GetDaughter(1)]:
        mc_idx = event.fitTrack2MC[daughter_idx]
        if mc_idx < 0: # no MC match for reco track, cannot check charge, skip
            continue
        pdg = event.MCTrack[mc_idx].GetPdgCode()
        if abs(pdg) not in {11, 13}:  # only check leptons
            continue
        # for e- (pdg=11) and mu- (pdg=13): truth charge is -1
        # for e+ (pdg=-11) and mu+ (pdg=-13): truth charge is +1
        truth_charge = -1 if pdg > 0 else +1
        reco_charge  = event.FitTracks[daughter_idx].getFittedState().getCharge()
        if int(round(reco_charge)) != truth_charge:
            misid = True
    return misid


def is_reconstructible_mc(event, mc_idx, ShipGeo, veto_geo):
    """
    Check if a MC track is reconstructible:
    >= 25 straw hits AND hits in >= 3 stations and vertex in fiducial volume.
    """
    n_hits = 0
    stations_hit = set()
    for hit in event.strawtubesPoint:
        if hit.GetTrackID() == mc_idx:
            n_hits += 1
            det_id_str = str(hit.GetDetectorID())
            station = int(det_id_str[0])
            stations_hit.add(station)

    # vertex in fiducial volume?
    vtx = ROOT.TVector3(
        event.MCTrack[1].GetStartX(),
        event.MCTrack[1].GetStartY(),
        event.MCTrack[1].GetStartZ()
    )
    if vtx.Z() > ShipGeo.TrackStation1.z:
        return False
    if vtx.Z() < veto_geo.z0:
        return False
    vertex_node = ROOT.gGeoManager.FindNode(vtx.X(), vtx.Y(), vtx.Z())
    if not vertex_node:
        return False
    if not vertex_node.GetVolume().GetName().startswith("DecayVacuum_"):
        return False

    return n_hits >= 25 and len(stations_hit) >= 3


def compute_vertexing_efficiency(event, sgeo, h, vtx_eff_counts, ShipGeo, veto_geo):
    z_mc = event.MCTrack[0].GetStartZ()

    # --- finde alle Töchter des DIS-Vertex (MotherID == 1) ---
    reco_daughters = []
    for i, mc_track in enumerate(event.MCTrack):
        if mc_track.GetMotherId() == 1:
            reco_daughters.append(i)

    if len(reco_daughters) < 2:
        return  # nicht genug Töchter für einen Vertex

    # --- rekonstruierbar wenn mind. 2 Töchter die Kriterien erfüllen ---
    n_reconstructible = sum(
        is_reconstructible_mc(event, mc_idx, ShipGeo, veto_geo)
        for mc_idx in reco_daughters
    )
    if n_reconstructible < 2:
        return

    vtx_eff_counts['reconstructible'] += 1
    h['vtx_eff_z_reconstructible'].Fill(z_mc)

    # --- vertex found? ---
    vertex_found = False
    for part in event.Particles:
        s1 = event.FitTracks[part.GetDaughter(0)].getFitStatus()
        s2 = event.FitTracks[part.GetDaughter(1)].getFitStatus()
        if s1.isFitConverged() and s2.isFitConverged():
            vertex_found = True
            break

    if vertex_found:
        vtx_eff_counts['found'] += 1
        h['vtx_eff_z_found'].Fill(z_mc)

        for part in event.Particles:
            if check_charge_misid(event, part):
                vtx_eff_counts['charge_misid'] += 1
                h['charge_misid_z'].Fill(z_mc)
                break


def persist_vertexing_efficiency(outfile_base, h, vtx_eff_counts):
    """
    Compute and save efficiency histogram + overall number.
    Call this at the end in Main_function, before writeHists.
    """
    #global vtx_eff_counts
    n_reco = vtx_eff_counts['reconstructible']
    n_found = vtx_eff_counts['found']
    n_misid = vtx_eff_counts['charge_misid']

    if n_reco == 0:
        print("Vertexing efficiency: no reconstructible events found.")
        return

    overall_eff = n_found / n_reco * 100.0
    misid_rate  = n_misid / n_reco * 100.0

    print(f"\n{'='*50}")
    print(f"Vertexing efficiency summary:")
    print(f"  Reconstructible events : {n_reco}")
    print(f"  Vertex found           : {n_found}  ({overall_eff:.1f}%)")
    print(f"  Charge mis-ID          : {n_misid}  ({misid_rate:.1f}% of reconstructible)")
    print(f"{'='*50}\n")

    # --- compute efficiency ratio histogram ---
    eff_hist = h['vtx_eff_z_ratio']
    eff_hist.Divide(h['vtx_eff_z_found'], h['vtx_eff_z_reconstructible'], 1, 1, 'B')
    # 'B' uses binomial errors

    # --- save summary to txt ---
    txt_path = outfile_base + 'VtxEff.txt'
    try:
        with open(txt_path, 'w') as f:
            f.write(f"Vertexing efficiency summary\n")
            f.write(f"Reconstructible events : {n_reco}\n")
            f.write(f"Vertex found           : {n_found}  ({overall_eff:.1f}%)\n")
            f.write(f"Charge mis-ID          : {n_misid}  ({misid_rate:.1f}% of reconstructible)\n")
        print(f"Vertexing efficiency summary saved to {txt_path}")
    except Exception as e:
        print(f"Warning: could not save vertexing efficiency to {txt_path} ({e})")