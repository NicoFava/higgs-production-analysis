#!/usr/bin/env python3
import argparse
import os
import sys

try:
    import ROOT
except ImportError as e:
    print('ERROR: PyROOT is required.')
    print(e)
    sys.exit(1)

# Set ROOT to batch mode
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

def get_event_weight(tree, branches, is_data=False):
    """
    Computes the event weight based on the presence of 'mcWeight' and scale factors.
    If the event is from data, returns 1.0.
    """
    if is_data:
        return 1.0
    w = 1.0
    if 'mcWeight' in branches:
        w *= float(tree.mcWeight)
    for b_name in branches:
        if b_name.startswith('scaleFactor_'):
            w *= float(getattr(tree, b_name))
    return w

def branch_names(tree):
    return {b.GetName() for b in tree.GetListOfBranches()}

def analyze_single_file(file_path, dataset_color, tree_name="mini"):
    """
    Analyzes a single blind ROOT file, applies lepton pairing for >= 4 leptons,
    and returns histograms for m_4l, m_Z1, and m_Z2.
    """
    root_file = ROOT.TFile.Open(file_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"Error: Cannot open file {file_path}")
        return None
        
    tree = root_file.Get(tree_name)
    if not tree:
        print(f"Error: Cannot find tree '{tree_name}' in {file_path}")
        return None
        
    base_name = os.path.basename(file_path).replace(".root", "")
    print(f"Processing {base_name}... Total events: {tree.GetEntries()}")

    Z_MASS_GEV = 91.18 
    
    # Initialize histograms
    h_m4l = ROOT.TH1F(f"h_m4l_{base_name}", f"Higgs Candidate Mass ({base_name}); m_{{4l}} [GeV]; Events", 30, 110, 140)
    h_mZ1 = ROOT.TH1F(f"h_mZ1_{base_name}", f"Leading Z Boson Mass ({base_name}); m_{{Z1}} [GeV]; Events", 40, 40, 120)
    h_mZ2 = ROOT.TH1F(f"h_mZ2_{base_name}", f"Subleading Z Boson Mass ({base_name}); m_{{Z2}} [GeV]; Events", 50, 10, 110)
    
    # Apply color styling to histograms
    for hist in [h_m4l, h_mZ1, h_mZ2]:
        hist.SetLineColor(dataset_color)
        hist.SetFillColor(dataset_color)
        hist.SetFillStyle(3004)
        hist.SetLineWidth(2)
        hist.SetDirectory(0)

    branches = branch_names(tree)

    # EVENT LOOP
    for event in tree:
        n_leps = getattr(event, 'lep_n', 0)
        w = get_event_weight(tree, branches, is_data=False)
        if n_leps < 4:
            continue
            
        # Build Lorentz Vectors for all leptons in the event
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
            
        pairs = []
        # OSSF pair (Opposite Sign, Same Flavor)
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
        final_z1_vec = None
        final_z2_vec = None

        for p1_idx in range(len(pairs)):
            for p2_idx in range(p1_idx + 1, len(pairs)):
                p1 = pairs[p1_idx]
                p2 = pairs[p2_idx]

                if len({p1['idx1'], p1['idx2'], p2['idx1'], p2['idx2']}) == 4:
                    dist1 = abs(p1['mass'] - Z_MASS_GEV)
                    dist2 = abs(p2['mass'] - Z_MASS_GEV)

                    if dist1 < dist2:
                        current_z1_diff = dist1
                        cand_z1 = p1
                        cand_z2 = p2
                    else:
                        current_z1_diff = dist2
                        cand_z1 = p2
                        cand_z2 = p1

                    if current_z1_diff < best_z1_diff:
                        best_z1_diff = current_z1_diff
                        final_z1_vec = cand_z1['vec_sum']
                        final_z2_vec = cand_z2['vec_sum']

        if final_z1_vec is None or final_z2_vec is None:
            continue
            
        # Reconstruct Higgs candidate kinematics
        higgs_vec = final_z1_vec + final_z2_vec
        m_4l = higgs_vec.M()
        
        # Fill histograms
        h_m4l.Fill(m_4l, w)
        h_mZ1.Fill(final_z1_vec.M(), w)
        h_mZ2.Fill(final_z2_vec.M(), w)

    root_file.Close()
    return h_m4l, h_mZ1, h_mZ2


def main():
    parser = argparse.ArgumentParser(description="Independent analysis of 4 blind ATLAS Open Data samples.")
    parser.add_argument("-i", "--input", nargs='+', required=True, 
                        help="List of the 4 input ROOT files.")
    parser.add_argument("-o", "--output", default="output_plots", 
                        help="Directory to save the larger PNG plots.")
    
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Defined specific color wheel for the 4 distinct datasets
    colors = [ROOT.kBlue+2, ROOT.kRed+1, ROOT.kGreen+2, ROOT.kOrange+7]

    for idx, file_path in enumerate(args.input):
        if not os.path.exists(file_path):
            print(f"Warning: File '{file_path}' not found. Skipping.")
            continue
            
        # Assign a color based on file index loop
        color = colors[idx % len(colors)]
        hist_tuple = analyze_single_file(file_path, color)
        
        if hist_tuple is None:
            continue
            
        h_m4l, h_mZ1, h_mZ2 = hist_tuple
        base_name = os.path.basename(file_path).replace(".root", "")
        
        # Plotting block for each variable
        hist_dict = {"m4l": h_m4l, "mZ1": h_mZ1, "mZ2": h_mZ2}
        
        for var_name, hist in hist_dict.items():
            canvas = ROOT.TCanvas(f"c_{base_name}_{var_name}", "", 1200, 800)
            
            # Pad margin optimization to prevent vertical axis label clipping
            ROOT.gPad.SetLeftMargin(0.20)
            ROOT.gPad.SetBottomMargin(0.12)
            
            # Refine title and label offsets dynamically
            hist.GetYaxis().SetTitleOffset(1.5)
            hist.GetXaxis().SetTitleOffset(1.1)
            
            hist.Draw("HIST")
            
            output_file = os.path.join(args.output, f"{base_name}_{var_name}.png")
            canvas.SaveAs(output_file)
            print(f"Saved: {output_file}")
        print()

if __name__ == "__main__":
    main()