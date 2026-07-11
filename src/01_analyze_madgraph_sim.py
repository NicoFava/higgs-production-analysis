#!/usr/bin/env python3
import argparse
import math
import os
import sys

try:
    import ROOT
except ImportError as e:
    print('ERROR: PyROOT is required to run this script.')
    print(e)
    sys.exit(1)

# Disable graphic output on screen for faster batch processing
ROOT.gROOT.SetBatch(True)
# Ensure correct error tracking when scaling histograms
ROOT.TH1.SetDefaultSumw2(True)

def main():
    # ==============================================================================
    # 1. ARGUMENT PARSING
    # ==============================================================================
    parser = argparse.ArgumentParser(description="Analyze Higgs production at Parton and Shower levels.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Path to the input Delphes ROOT file")
    parser.add_argument('-o', '--outdir', type=str, required=True, help="Output directory for the generated plots")
    parser.add_argument('-p', '--process', type=str, required=True, help="Name of the physical process (e.g., ggF, VBF, WH, ZH)")
    args = parser.parse_args()

    # ==============================================================================
    # 2. LOAD DELPHES LIBRARIES
    # ==============================================================================
    DELPHES_DIR = os.environ.get("DELPHES_DIR")

    if DELPHES_DIR is None:
        # Try common locations
        possible_paths = [
            os.path.expanduser("~/Delphes"),
            "/usr/local/Delphes",
            "/opt/Delphes",
            os.path.expanduser("~/MG5_aMC_v3_5_3/Delphes") 
        ]

        for path in possible_paths:
            if os.path.isdir(path):
                DELPHES_DIR = path
                break

    if DELPHES_DIR is None:
        raise EnvironmentError(
            "Delphes not found. Please set the DELPHES_DIR environment variable."
        )

    DELPHES_LIB = os.path.join(DELPHES_DIR, "libDelphes.so")

    ROOT.gInterpreter.AddIncludePath(DELPHES_DIR)
    ROOT.gInterpreter.AddIncludePath(os.path.join(DELPHES_DIR, "external"))

    os.environ["ROOT_INCLUDE_PATH"] = (
        DELPHES_DIR + ":" + os.path.join(DELPHES_DIR, "external")
    )

    ROOT.gSystem.Load(DELPHES_LIB)

    ROOT.gInterpreter.Declare(r'''
    #include "classes/SortableObject.h"
    #include "classes/DelphesClasses.h"
    #include "ExRootAnalysis/ExRootTreeReader.h"
    ''')

    # ==============================================================================
    # 3. INPUT FILE AND TREE SETUP
    # ==============================================================================
    if not os.path.exists(args.input):
        print(f"ERROR: Input ROOT file not found at: {args.input}")
        sys.exit(1)

    f = ROOT.TFile.Open(args.input)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open input file: {args.input}")

    tree = f.Get("Delphes")
    if not tree:
        raise RuntimeError("Could not find tree 'Delphes' in the input file.")

    reader = ROOT.ExRootTreeReader(tree)
    branch_particle = reader.UseBranch("Particle")
    branch_genjet = reader.UseBranch("GenJet")

    n_entries = reader.GetEntries()
    print(f"Total events to process = {n_entries}")

    # ==============================================================================
    # 4. HISTOGRAM DEFINITIONS
    # ==============================================================================
    h = {"parton": {}, "shower": {}}

    for level in ["parton", "shower"]:
        prefix = f"{level}_{args.process}"
        
        # General observables
        h[level]["Njets"] = ROOT.TH1F(f"Njets_{prefix}", f"{args.process}: Number of Jets; N_{{jets}}; Events", 10, -0.5, 9.5)
        h[level]["pT_H"]  = ROOT.TH1F(f"pT_H_{prefix}", f"{args.process}: Higgs Boson p_{{T}}; p_{{T}}^{{H}} [GeV]; Events", 50, 0, 400)
        h[level]["eta_H"] = ROOT.TH1F(f"eta_H_{prefix}", f"{args.process}: Higgs Boson #eta; #eta_{{H}}; Events", 50, -5, 5)
        h[level]["HT"]    = ROOT.TH1F(f"HT_{prefix}", f"{args.process}: Scalar Sum of p_{{T}}; H_{{T}} [GeV]; Events", 50, 0, 800)
        
        # VBF Specific Observables (Requires >= 2 jets)
        h[level]["m_jj"]    = ROOT.TH1F(f"m_jj_{prefix}", f"{args.process}: Dijet Invariant Mass; m_{{jj}} [GeV]; Events", 50, 0, 2500)
        h[level]["dEta_jj"] = ROOT.TH1F(f"dEta_jj_{prefix}", f"{args.process}: Pseudorapidity Separation; |#Delta#eta_{{jj}}|; Events", 50, 0, 8)
        h[level]["dPhi_jj"] = ROOT.TH1F(f"dPhi_jj_{prefix}", f"{args.process}: Azimuthal Separation; #Delta#phi_{{jj}} [rad]; Events", 50, 0, ROOT.TMath.Pi())

        # VH Specific Observables (Requires Leptons / MET)
        h[level]["N_l"]     = ROOT.TH1F(f"N_l_{prefix}", f"{args.process}: Number of Leptons; N_{{#ell}}; Events", 6, -0.5, 5.5)
        h[level]["Emiss_T"] = ROOT.TH1F(f"Emiss_T_{prefix}", f"{args.process}: Missing Transverse Energy; E_{{T}}^{{miss}} [GeV]; Events", 50, 0, 400)
        h[level]["m_ll"]    = ROOT.TH1F(f"m_ll_{prefix}", f"{args.process}: Dilepton Invariant Mass; m_{{#ell#ell}} [GeV]; Events", 50, 0, 200)
        h[level]["m_WT"]    = ROOT.TH1F(f"m_WT_{prefix}", f"{args.process}: W Transverse Mass; m_{{T}}^{{W}} [GeV]; Events", 50, 0, 200)


    # ==============================================================================
    # 5. EVENT LOOP
    # ==============================================================================
    for ievt in range(n_entries):
        reader.ReadEntry(ievt)

        if ievt % 1000 == 0 and ievt > 0:
            print(f"Processed {ievt}/{n_entries} events...")

        # Collections to build physics objects
        higgs_candidates = []
        
        partonic_jets = []
        partonic_leptons = []
        partonic_met_x = 0.0
        partonic_met_y = 0.0
        
        shower_jets = []
        shower_leptons = []
        shower_met_x = 0.0
        shower_met_y = 0.0

        # --------------------------------------------------------------------------
        # A. PARSING GENERATOR PARTICLES (Reading 'Particle' branch)
        # --------------------------------------------------------------------------
        for ipart in range(branch_particle.GetEntries()):
            part = branch_particle.At(ipart)
            pid = abs(part.PID)

            # Store all Higgs copies to identify hard-process vs post-shower
            if part.PID == 25:
                higgs_candidates.append(part)
            
            # --- PARTON LEVEL (Status 23: Outgoing from hard process) ---
            if part.Status == 23:
                # Partonic Jets (Quarks 1-6 and Gluons 21)
                if pid <= 6 or pid == 21:
                    if part.PT > 25.0 and abs(part.Eta) < 4.5:
                        jet = ROOT.TLorentzVector()
                        jet.SetPtEtaPhiM(part.PT, part.Eta, part.Phi, part.Mass)
                        partonic_jets.append(jet)
                
                # Partonic Leptons (Electrons 11, Muons 13)
                elif pid in [11, 13]:
                    if part.PT > 10.0 and abs(part.Eta) < 2.5:
                        lep = ROOT.TLorentzVector()
                        lep.SetPtEtaPhiM(part.PT, part.Eta, part.Phi, part.Mass)
                        partonic_leptons.append(lep)
                
                # Partonic Neutrinos for MET (12, 14, 16)
                elif pid in [12, 14, 16]:
                    partonic_met_x += part.PT * math.cos(part.Phi)
                    partonic_met_y += part.PT * math.sin(part.Phi)

            # --- SHOWER LEVEL (Status 1: Final state particles before detector) ---
            if part.Status == 1:
                # Shower Leptons
                if pid in [11, 13]:
                    if part.PT > 10.0 and abs(part.Eta) < 2.5:
                        lep = ROOT.TLorentzVector()
                        lep.SetPtEtaPhiM(part.PT, part.Eta, part.Phi, part.Mass)
                        shower_leptons.append(lep)
                
                # Shower Neutrinos for MET
                elif pid in [12, 14, 16]:
                    shower_met_x += part.PT * math.cos(part.Phi)
                    shower_met_y += part.PT * math.sin(part.Phi)

        # --------------------------------------------------------------------------
        # B. SHOWER JETS (Reading 'GenJet' branch from Delphes)
        # --------------------------------------------------------------------------
        for ijet in range(branch_genjet.GetEntries()):
            genjet = branch_genjet.At(ijet)
            if genjet.PT > 25.0 and abs(genjet.Eta) < 4.5:
                jet = ROOT.TLorentzVector()
                jet.SetPtEtaPhiM(genjet.PT, genjet.Eta, genjet.Phi, genjet.Mass)
                shower_jets.append(jet)

        # --------------------------------------------------------------------------
        # C. HIGGS KINEMATICS EXTRACTION
        # --------------------------------------------------------------------------
        higgs_parton = None
        higgs_shower = None
        
        if len(higgs_candidates) > 0:
            # Parton level: First copy (Hard process, usually Status 22)
            h_part = higgs_candidates[0]
            higgs_parton = ROOT.TLorentzVector()
            higgs_parton.SetPtEtaPhiM(h_part.PT, h_part.Eta, h_part.Phi, h_part.Mass)
            
            # Shower level: Last copy before decay (Includes ISR recoil, usually Status 62)
            h_show = higgs_candidates[-1]
            higgs_shower = ROOT.TLorentzVector()
            higgs_shower.SetPtEtaPhiM(h_show.PT, h_show.Eta, h_show.Phi, h_show.Mass)

        # --------------------------------------------------------------------------
        # D. SORTING AND MET CALCULATION
        # --------------------------------------------------------------------------
        # Sort collections by descending pT
        partonic_jets.sort(key=lambda j: j.Pt(), reverse=True)
        shower_jets.sort(key=lambda j: j.Pt(), reverse=True)
        partonic_leptons.sort(key=lambda l: l.Pt(), reverse=True)
        shower_leptons.sort(key=lambda l: l.Pt(), reverse=True)

        # Reconstruct MET vectors
        part_met_vec = ROOT.TLorentzVector()
        part_met_vec.SetPxPyPzE(partonic_met_x, partonic_met_y, 0.0, math.sqrt(partonic_met_x**2 + partonic_met_y**2))
        
        show_met_vec = ROOT.TLorentzVector()
        show_met_vec.SetPxPyPzE(shower_met_x, shower_met_y, 0.0, math.sqrt(shower_met_x**2 + shower_met_y**2))

        # --------------------------------------------------------------------------
        # E. FILLING HISTOGRAMS
        # --------------------------------------------------------------------------
        
        # --- Helper function to avoid code duplication ---
        def fill_histograms(level, higgs, jets, leptons, met_vec):
            # 1. Higgs Observables
            if higgs:
                h[level]["pT_H"].Fill(higgs.Pt())
                h[level]["eta_H"].Fill(higgs.Eta())
            
            # 2. Jet Observables
            h[level]["Njets"].Fill(len(jets))
            h[level]["HT"].Fill(sum([j.Pt() for j in jets]))
            
            if len(jets) >= 2:
                j1, j2 = jets[0], jets[1]
                h[level]["m_jj"].Fill((j1 + j2).M())
                h[level]["dEta_jj"].Fill(abs(j1.Eta() - j2.Eta()))
                h[level]["dPhi_jj"].Fill(abs(j1.DeltaPhi(j2)))

            # 3. VH (Leptonic) Observables
            h[level]["N_l"].Fill(len(leptons))
            h[level]["Emiss_T"].Fill(met_vec.Pt())
            
            if len(leptons) >= 2:
                h[level]["m_ll"].Fill((leptons[0] + leptons[1]).M())
                
            if len(leptons) >= 1 and met_vec.Pt() > 0:
                l1 = leptons[0]
                dphi = abs(l1.DeltaPhi(met_vec))
                # mT = sqrt(2 * pT_l * MET * (1 - cos(DeltaPhi)))
                mt = math.sqrt(max(0.0, 2 * l1.Pt() * met_vec.Pt() * (1.0 - math.cos(dphi))))
                h[level]["m_WT"].Fill(mt)

        # Fill for both levels
        fill_histograms("parton", higgs_parton, partonic_jets, partonic_leptons, part_met_vec)
        fill_histograms("shower", higgs_shower, shower_jets, shower_leptons, show_met_vec)

    # ==============================================================================
    # 6. SAVING PLOTS WITH COMPARISON
    # ==============================================================================
    output_dir = os.path.join(args.outdir, args.process)
    os.makedirs(output_dir, exist_ok=True)
    canvas = ROOT.TCanvas("c", "Canvas", 800, 600)

    for var in h["parton"].keys():
        h["parton"][var].SetLineColor(ROOT.kBlue)
        h["parton"][var].SetLineWidth(2)
        
        h["shower"][var].SetLineColor(ROOT.kRed)
        h["shower"][var].SetLineWidth(2)
        
        # Normalize histograms to unity (Shape comparison)
        if h["parton"][var].Integral() > 0: 
            h["parton"][var].Scale(1.0 / h["parton"][var].Integral())
        if h["shower"][var].Integral() > 0: 
            h["shower"][var].Scale(1.0 / h["shower"][var].Integral())
        
        max_y = max(h["parton"][var].GetMaximum(), h["shower"][var].GetMaximum())
        h["parton"][var].SetMaximum(max_y * 1.3)
        
        # Plot styling
        h["parton"][var].GetYaxis().SetTitle("Normalized Units")
        h["parton"][var].SetStats(0)
        h["shower"][var].SetStats(0)

        h["parton"][var].Draw("HIST")
        h["shower"][var].Draw("HIST SAME")
        
        legend = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
        legend.AddEntry(h["parton"][var], f"{args.process}: Parton Level", "l")
        legend.AddEntry(h["shower"][var], f"{args.process}: Parton Shower", "l")
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.Draw()
        
        output_file = os.path.join(output_dir, f"{var}_{args.process}.png")
        canvas.SaveAs(output_file)
        
        # Clear canvas for next iteration
        canvas.Clear()

    print(f"Analysis complete! Plots for {args.process} saved in '{output_dir}'.")

if __name__ == "__main__":
    main()