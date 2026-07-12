#!/usr/bin/env python3
import argparse
import os
import sys
import math

try:
    import ROOT
except ImportError as e:
    print('ERROR: PyROOT is required to run this script.')
    print(e)
    sys.exit(1)

# Set ROOT to batch mode for faster processing without pop-ups
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

def analyze_single_file(file_path, output_dir, dataset_color=ROOT.kBlack, tree_name="mini"):
    """
    Analyzes a single ATLAS Open Data ROOT file, selects the Higgs 
    candidate (4 leptons), and computes the 11 kinematic variables 
    for the associated objects. Saves results in both PNG and ROOT formats.
    """
    root_file = ROOT.TFile.Open(file_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"Error: Cannot open file {file_path}")
        return
        
    tree = root_file.Get(tree_name)
    if not tree:
        print(f"Error: Cannot find tree '{tree_name}' in {file_path}")
        return
        
    base_name = os.path.basename(file_path).replace(".root", "")
    print(f"Processing {base_name}... Total events: {tree.GetEntries()}")

    Z_MASS_GEV = 91.18 
    
    # 1. Initialize the 11 Histograms as defined in 01_analyze_madgraph_sim.py
    h = {}
    
    # General observables
    h["Njets"] = ROOT.TH1F(f"Njets_{base_name}", f"Number of Jets; N_{{jets}}; Events", 10, -0.5, 9.5)
    h["pT_H"]  = ROOT.TH1F(f"pT_H_{base_name}", f"Higgs Boson p_{{T}}; p_{{T}}^{{H}} [GeV]; Events", 50, 0, 400)
    h["eta_H"] = ROOT.TH1F(f"eta_H_{base_name}", f"Higgs Boson #eta; #eta_{{H}}; Events", 50, -5, 5)
    h["HT"]    = ROOT.TH1F(f"HT_{base_name}", f"Scalar Sum of Jet p_{{T}}; H_{{T}} [GeV]; Events", 50, 0, 800)
    
    # VBF Specific Observables
    h["m_jj"]    = ROOT.TH1F(f"m_jj_{base_name}", f"Dijet Invariant Mass; m_{{jj}} [GeV]; Events", 50, 0, 2500)
    h["dEta_jj"] = ROOT.TH1F(f"dEta_jj_{base_name}", f"Pseudorapidity Separation; |#Delta#eta_{{jj}}|; Events", 50, 0, 8)
    h["dPhi_jj"] = ROOT.TH1F(f"dPhi_jj_{base_name}", f"Azimuthal Separation; #Delta#phi_{{jj}} [rad]; Events", 50, 0, ROOT.TMath.Pi())

    # VH Specific Observables (Associated Leptons and MET)
    h["N_l"]     = ROOT.TH1F(f"N_l_{base_name}", f"Number of Associated Leptons; N_{{extra l}}; Events", 6, -0.5, 5.5)
    h["Emiss_T"] = ROOT.TH1F(f"Emiss_T_{base_name}", f"Missing Transverse Energy; E_{{T}}^{{miss}} [GeV]; Events", 50, 0, 400)
    h["m_ll"]    = ROOT.TH1F(f"m_ll_{base_name}", f"Associated Dilepton Mass; m_{{ll}} [GeV]; Events", 50, 0, 200)
    h["m_WT"]    = ROOT.TH1F(f"m_WT_{base_name}", f"Associated W Transverse Mass; m_{{T}}^{{W}} [GeV]; Events", 50, 0, 200)

    # Style histograms
    for hist in h.values():
        hist.SetLineColor(dataset_color)
        hist.SetFillColor(dataset_color)
        hist.SetFillStyle(3004)
        hist.SetLineWidth(2)
        hist.SetDirectory(0)

    # 2. EVENT LOOP
    for event in tree:
        n_leps = getattr(event, 'lep_n', 0)
        if n_leps < 4:
            continue
            
        # Build Lorentz Vectors for all leptons (assuming pT, E are in MeV)
        leptons = []
        for i in range(n_leps):
            vec = ROOT.TLorentzVector()
            vec.SetPtEtaPhiE(event.lep_pt[i] / 1000.0, 
                             event.lep_eta[i], 
                             event.lep_phi[i], 
                             event.lep_E[i] / 1000.0)
            leptons.append({
                'idx': i,
                'vec': vec,
                'charge': event.lep_charge[i],
                'type': event.lep_type[i]
            })
            
        # Identify OSSF pairs to find Z candidates
        pairs = []
        for i in range(len(leptons)):
            for j in range(i + 1, len(leptons)):
                l1 = leptons[i]
                l2 = leptons[j]
                if (l1['type'] == l2['type']) and (l1['charge'] * l2['charge'] < 0):
                    mass = (l1['vec'] + l2['vec']).M()
                    pairs.append({
                        'idx1': l1['idx'], 
                        'idx2': l2['idx'],
                        'vec_sum': l1['vec'] + l2['vec'],
                        'mass': mass
                    })

        if len(pairs) < 2:
            continue

        best_z1_diff = 999999.0
        cand_z1 = None
        cand_z2 = None

        # Select the best combination of 4 leptons
        for p1_idx in range(len(pairs)):
            for p2_idx in range(p1_idx + 1, len(pairs)):
                p1 = pairs[p1_idx]
                p2 = pairs[p2_idx]

                if len({p1['idx1'], p1['idx2'], p2['idx1'], p2['idx2']}) == 4:
                    dist1 = abs(p1['mass'] - Z_MASS_GEV)
                    dist2 = abs(p2['mass'] - Z_MASS_GEV)

                    if dist1 < dist2:
                        current_z1_diff = dist1
                        tmp_z1, tmp_z2 = p1, p2
                    else:
                        current_z1_diff = dist2
                        tmp_z1, tmp_z2 = p2, p1

                    if current_z1_diff < best_z1_diff:
                        best_z1_diff = current_z1_diff
                        cand_z1 = tmp_z1
                        cand_z2 = tmp_z2

        if cand_z1 is None or cand_z2 is None:
            continue
            
        # Higgs Reconstruction
        higgs_vec = cand_z1['vec_sum'] + cand_z2['vec_sum']
        m_4l = higgs_vec.M()
        
        # Higgs mass window cut: Isolate the Higgs signal (e.g., 115 < m_4l < 130 GeV)
        if not (115.0 < m_4l < 130.0):
            continue

        # Extract extra associated leptons (not part of the Higgs candidate)
        higgs_lep_indices = {cand_z1['idx1'], cand_z1['idx2'], cand_z2['idx1'], cand_z2['idx2']}
        extra_leptons = [l['vec'] for l in leptons if l['idx'] not in higgs_lep_indices]
        extra_leptons.sort(key=lambda l: l.Pt(), reverse=True)

        # Process Jets (assuming ATLAS Open Data pT, E are in MeV)
        n_jets = getattr(event, 'jet_n', 0)
        jets = []
        for i in range(n_jets):
            pt = event.jet_pt[i] / 1000.0
            eta = event.jet_eta[i]
            # Standard jet selection threshold
            if pt > 30.0 and abs(eta) < 4.4:
                vec = ROOT.TLorentzVector()
                vec.SetPtEtaPhiE(pt, eta, event.jet_phi[i], event.jet_E[i] / 1000.0)
                jets.append(vec)
        
        jets.sort(key=lambda j: j.Pt(), reverse=True)

        # Process Missing Transverse Energy (MET)
        met_et = getattr(event, 'met_et', 0) / 1000.0
        met_phi = getattr(event, 'met_phi', 0)
        met_vec = ROOT.TLorentzVector()
        met_vec.SetPtEtaPhiM(met_et, 0, met_phi, 0)

        # 3. FILL HISTOGRAMS
        h["pT_H"].Fill(higgs_vec.Pt())
        h["eta_H"].Fill(higgs_vec.Eta())
        
        h["Njets"].Fill(len(jets))
        h["HT"].Fill(sum([j.Pt() for j in jets]))
        
        if len(jets) >= 2:
            j1, j2 = jets[0], jets[1]
            h["m_jj"].Fill((j1 + j2).M())
            h["dEta_jj"].Fill(abs(j1.Eta() - j2.Eta()))
            h["dPhi_jj"].Fill(abs(j1.DeltaPhi(j2)))

        h["N_l"].Fill(len(extra_leptons))
        h["Emiss_T"].Fill(met_vec.Pt())
        
        if len(extra_leptons) >= 2:
            h["m_ll"].Fill((extra_leptons[0] + extra_leptons[1]).M())
            
        if len(extra_leptons) >= 1 and met_vec.Pt() > 0:
            l1 = extra_leptons[0]
            dphi = abs(l1.DeltaPhi(met_vec))
            mt = math.sqrt(max(0.0, 2 * l1.Pt() * met_vec.Pt() * (1.0 - math.cos(dphi))))
            h["m_WT"].Fill(mt)

    root_file.Close()

    # 4. SAVE RESULTS
    # Write to a ROOT file for later overlay with MadGraph
    out_root_path = os.path.join(output_dir, f"{base_name}_kinematics.root")
    out_root_file = ROOT.TFile.Open(out_root_path, "RECREATE")
    for hist in h.values():
        hist.Write()
    out_root_file.Close()
    print(f"Saved ROOT file: {out_root_path}")

    # Generate and save PNG files
    canvas = ROOT.TCanvas("c", "", 800, 600)
    for var_name, hist in h.items():
        hist.Draw("HIST")
        png_path = os.path.join(output_dir, f"{base_name}_{var_name}.png")
        canvas.SaveAs(png_path)

def main():
    parser = argparse.ArgumentParser(description="Extract 11 kinematic variables from ATLAS Open Data focusing on associated Higgs production.")
    parser.add_argument("-i", "--input", nargs='+', required=True, help="List of input ROOT files (ATLAS data or MC background).")
    parser.add_argument("-o", "--output", default="output_atlas_kinematics", help="Directory to save the PNG and .root files.")
    
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Define specific color wheel for the 4 distinct datasets
    colors = [ROOT.kBlue+2, ROOT.kRed+1, ROOT.kGreen+2, ROOT.kOrange+7]

    for idx, file_path in enumerate(args.input):
        if not os.path.exists(file_path):
            print(f"Warning: File '{file_path}' not found. Skipping.")
            continue
        # Assign a color based on file index
        color = colors[idx % len(colors)]
        analyze_single_file(file_path, args.output, color)

if __name__ == "__main__":
    main()