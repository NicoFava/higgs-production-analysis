#!/usr/bin/env python3
import argparse
import math
import os
import sys
import gzip

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

# ==============================================================================
# LOAD DELPHES LIBRARY AND HEADERS
# ==============================================================================
DELPHES_DIR = os.environ.get("DELPHES_DIR")

if DELPHES_DIR is None:
    possible_paths = [
        os.path.expanduser("~/Delphes"),
        os.path.expanduser("~/lab/Delphes"),
        "/usr/local/Delphes",
        "/opt/Delphes",
    ]
    for path in possible_paths:
        if os.path.isdir(path):
            DELPHES_DIR = path
            break

if DELPHES_DIR is None:
    raise EnvironmentError("Delphes not found. Please set DELPHES_DIR environment variable.")

DELPHES_LIB = os.path.join(DELPHES_DIR, "libDelphes.so")

ROOT.gInterpreter.AddIncludePath(DELPHES_DIR)
ROOT.gInterpreter.AddIncludePath(os.path.join(DELPHES_DIR, "external"))

os.environ["ROOT_INCLUDE_PATH"] = DELPHES_DIR + ":" + os.path.join(DELPHES_DIR, "external")
ROOT.gSystem.Load(DELPHES_LIB)

ROOT.gInterpreter.Declare(r'''
#include "classes/SortableObject.h"
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
''')

# ==============================================================================
# PARSING FUNCTIONS FOR LHE
# ==============================================================================
def read_lhe_events(filepath):
    """Yields (event_particles, event_weight) from a zipped LHE file."""
    with gzip.open(filepath, 'rt') as f:
        in_event = False
        is_first_line = False
        event_particles = []
        event_weight = 1.0
        
        for line in f:
            line = line.strip()
            if not line: continue
            
            if line.startswith("<event>"):
                in_event = True
                is_first_line = True
                event_particles = []
                event_weight = 1.0
                continue
                
            if line.startswith("</event>"):
                in_event = False
                yield event_particles, event_weight
                continue
                
            if in_event:
                # The first line inside <event> contains: NUP IDPRUP XWGTUP SCALUP AQEDUP AQCDUP
                if is_first_line:
                    parts = line.split()
                    if len(parts) >= 3:
                        event_weight = float(parts[2])
                    is_first_line = False
                    continue
                
                parts = line.split()
                # Ensure it's a valid particle line
                if len(parts) >= 11 and not line.startswith('#'):
                    pdg = int(parts[0])
                    status = int(parts[1])
                    px, py, pz, e, mass = map(float, parts[6:11])
                    event_particles.append({
                        'pdg': pdg, 'status': status,
                        'px': px, 'py': py, 'pz': pz, 'e': e, 'm': mass
                    })

# ==============================================================================
# MAIN ANALYSIS LOOP
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(description="Analyze Higgs production from .lhe and Delphes .root directly.")
    parser.add_argument('--lhe', type=str, required=True, help="Path to the input .lhe.gz file (Parton Level)")
    parser.add_argument('--delphes', type=str, required=True, help="Path to the input Delphes .root file (Shower/Detector Level)")
    parser.add_argument('-o', '--outdir', type=str, required=True, help="Output directory for the generated plots")
    parser.add_argument('-p', '--process', type=str, required=True, help="Name of the physical process (e.g., ggF, VBF, WH, ZH)")
    args = parser.parse_args()

    if not os.path.exists(args.lhe):
        raise RuntimeError(f"Could not find LHE file: {args.lhe}")
    if not os.path.exists(args.delphes):
        raise RuntimeError(f"Could not find Delphes file: {args.delphes}")

    # --------------------------------------------------------------------------
    # 1. HISTOGRAM DEFINITIONS
    # --------------------------------------------------------------------------
    h = {"parton": {}, "shower": {}}

    for level in ["parton", "shower"]:
        prefix = f"{level}_{args.process}"
        
        # General observables
        h[level]["Njets"] = ROOT.TH1F(f"Njets_{prefix}", f"{args.process}: Number of Jets; N_{{jets}}; Events", 10, -0.5, 9.5)
        h[level]["pT_H"]  = ROOT.TH1F(f"pT_H_{prefix}", f"{args.process}: Higgs Boson p_{{T}}; p_{{T}}^{{H}} [GeV]; Events", 50, 0, 400)
        h[level]["eta_H"] = ROOT.TH1F(f"eta_H_{prefix}", f"{args.process}: Higgs Boson #eta; #eta_{{H}}; Events", 50, -5, 5)
        h[level]["HT"]    = ROOT.TH1F(f"HT_{prefix}", f"{args.process}: Scalar Sum of Jet p_{{T}}; H_{{T}} [GeV]; Events", 50, 0, 800)
        
        # VBF Specific Observables (Requires >= 2 jets)
        h[level]["m_jj"]    = ROOT.TH1F(f"m_jj_{prefix}", f"{args.process}: Dijet Invariant Mass; m_{{jj}} [GeV]; Events", 50, 0, 2500)
        h[level]["dEta_jj"] = ROOT.TH1F(f"dEta_jj_{prefix}", f"{args.process}: Pseudorapidity Separation; |#Delta#eta_{{jj}}|; Events", 50, 0, 8)
        h[level]["dPhi_jj"] = ROOT.TH1F(f"dPhi_jj_{prefix}", f"{args.process}: Azimuthal Separation; #Delta#phi_{{jj}} [rad]; Events", 50, 0, ROOT.TMath.Pi())

        # VH Specific Observables (Requires Leptons / MET)
        h[level]["N_l"]     = ROOT.TH1F(f"N_l_{prefix}", f"{args.process}: Number of Leptons (N_{{l}}); N_{{l}}; Events", 6, -0.5, 5.5)
        h[level]["Emiss_T"] = ROOT.TH1F(f"Emiss_T_{prefix}", f"{args.process}: Missing Transverse Energy; E_{{T}}^{{miss}} [GeV]; Events", 50, 0, 400)
        h[level]["m_ll"]    = ROOT.TH1F(f"m_ll_{prefix}", f"{args.process}: Dilepton Invariant Mass (m_{{ll}}); m_{{ll}} [GeV]; Events", 50, 0, 200)
        h[level]["m_WT"]    = ROOT.TH1F(f"m_WT_{prefix}", f"{args.process}: W Transverse Mass; m_{{T}}^{{W}} [GeV]; Events", 50, 0, 200)

    # Helper function to avoid repetition (Now includes weight)
    def fill_histograms(level, higgs, jets, leptons, met_vec, weight):
        if higgs:
            h[level]["pT_H"].Fill(higgs.Pt(), weight)
            h[level]["eta_H"].Fill(higgs.Eta(), weight)
        
        h[level]["Njets"].Fill(len(jets), weight)
        h[level]["HT"].Fill(sum([j.Pt() for j in jets]), weight)
        
        if len(jets) >= 2:
            j1, j2 = jets[0], jets[1]
            h[level]["m_jj"].Fill((j1 + j2).M(), weight)
            h[level]["dEta_jj"].Fill(abs(j1.Eta() - j2.Eta()), weight)
            h[level]["dPhi_jj"].Fill(abs(j1.DeltaPhi(j2)), weight)

        h[level]["N_l"].Fill(len(leptons), weight)
        h[level]["Emiss_T"].Fill(met_vec.Pt(), weight)
        
        if len(leptons) >= 2:
            h[level]["m_ll"].Fill((leptons[0] + leptons[1]).M(), weight)
            
        if len(leptons) >= 1 and met_vec.Pt() > 0:
            l1 = leptons[0]
            dphi = abs(l1.DeltaPhi(met_vec))
            mt = math.sqrt(max(0.0, 2 * l1.Pt() * met_vec.Pt() * (1.0 - math.cos(dphi))))
            h[level]["m_WT"].Fill(mt, weight)

    # --------------------------------------------------------------------------
    # 2. PROCESS LHE FILE (PARTON LEVEL)
    # --------------------------------------------------------------------------
    print("Processing LHE events (Parton Level)...")
    for ievt, (event, weight) in enumerate(read_lhe_events(args.lhe)):
        if ievt % 1000 == 0 and ievt > 0:
            print(f"  Processed {ievt} LHE events...")
            
        higgs_parton = None
        jets = []
        leptons = []
        met_x, met_y = 0.0, 0.0

        for p in event:
            # Higgs
            if abs(p['pdg']) == 25:
                higgs_parton = ROOT.TLorentzVector()
                higgs_parton.SetPxPyPzE(p['px'], p['py'], p['pz'], p['e'])
            
            # Status 1 = outgoing final state partons in LHE
            if p['status'] == 1:
                pid = abs(p['pdg'])
                vec = ROOT.TLorentzVector()
                vec.SetPxPyPzE(p['px'], p['py'], p['pz'], p['e'])
                
                # Partonic Jets (Quarks 1-6 and Gluons 21)
                if (pid <= 6 or pid == 21) and vec.Pt() > 25.0 and abs(vec.Eta()) < 4.5:
                    jets.append(vec)
                # Leptons
                elif pid in [11, 13] and vec.Pt() > 10.0 and abs(vec.Eta()) < 2.5:
                    leptons.append(vec)
                # Neutrinos (MET)
                elif pid in [12, 14, 16]:
                    met_x += p['px']
                    met_y += p['py']

        jets.sort(key=lambda j: j.Pt(), reverse=True)
        leptons.sort(key=lambda l: l.Pt(), reverse=True)
        met_vec = ROOT.TLorentzVector()
        met_vec.SetPxPyPzE(met_x, met_y, 0.0, math.sqrt(met_x**2 + met_y**2))

        fill_histograms("parton", higgs_parton, jets, leptons, met_vec, weight)

    # --------------------------------------------------------------------------
    # 3. PROCESS DELPHES ROOT FILE (SHOWER / DETECTOR LEVEL)
    # --------------------------------------------------------------------------
    print("Processing Delphes events (Shower/Detector Level)...")
    
    chain = ROOT.TChain("Delphes")
    chain.Add(args.delphes)
    
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()
    
    # Extract branches
    branchEvent     = treeReader.UseBranch("Event")
    branchParticle  = treeReader.UseBranch("Particle")
    branchJet       = treeReader.UseBranch("Jet")
    branchElectron  = treeReader.UseBranch("Electron")
    branchMuon      = treeReader.UseBranch("Muon")
    branchMissingET = treeReader.UseBranch("MissingET")

    for entry in range(numberOfEntries):
        treeReader.ReadEntry(entry)
        
        if entry % 1000 == 0 and entry > 0:
            print(f"  Processed {entry} Delphes events...")

        # Extract Event Weight
        event_weight = 1.0
        if branchEvent.GetEntries() > 0:
            event_weight = branchEvent.At(0).Weight

        higgs_shower = None
        shower_jets = []
        leptons = []
        met_vec = ROOT.TLorentzVector()

        # 3.1 Find the Higgs from the Generator Particles
        for i in range(branchParticle.GetEntries()):
            p = branchParticle.At(i)
            if abs(p.PID) == 25:
                higgs_shower = ROOT.TLorentzVector()
                higgs_shower.SetPtEtaPhiM(p.PT, p.Eta, p.Phi, p.Mass)

        # 3.2 Collect Reconstructed Leptons
        for i in range(branchElectron.GetEntries()):
            e = branchElectron.At(i)
            if e.PT > 10.0 and abs(e.Eta) < 2.5:
                vec = ROOT.TLorentzVector()
                vec.SetPtEtaPhiM(e.PT, e.Eta, e.Phi, 0.000511)
                leptons.append(vec)

        for i in range(branchMuon.GetEntries()):
            mu = branchMuon.At(i)
            if mu.PT > 10.0 and abs(mu.Eta) < 2.5:
                vec = ROOT.TLorentzVector()
                vec.SetPtEtaPhiM(mu.PT, mu.Eta, mu.Phi, 0.10566)
                leptons.append(vec)

        # 3.3 Collect Reconstructed Jets
        for i in range(branchJet.GetEntries()):
            j = branchJet.At(i)
            if j.PT > 25.0 and abs(j.Eta) < 4.5:
                vec = ROOT.TLorentzVector()
                vec.SetPtEtaPhiM(j.PT, j.Eta, j.Phi, j.Mass)
                shower_jets.append(vec)

        # 3.4 Collect Missing Transverse Energy (MET)
        if branchMissingET.GetEntries() > 0:
            met = branchMissingET.At(0)
            met_vec.SetPtEtaPhiM(met.MET, 0.0, met.Phi, 0.0)

        # Sort collections by pT descending
        shower_jets.sort(key=lambda j: j.Pt(), reverse=True)
        leptons.sort(key=lambda l: l.Pt(), reverse=True)

        fill_histograms("shower", higgs_shower, shower_jets, leptons, met_vec, event_weight)

    # --------------------------------------------------------------------------
    # 4. PLOTTING & SAVING
    # --------------------------------------------------------------------------
    output_dir = os.path.join(args.outdir, args.process)
    os.makedirs(output_dir, exist_ok=True)
    canvas = ROOT.TCanvas("c", "Canvas", 800, 600)

    for var in h["parton"].keys():
        h["parton"][var].SetLineColor(ROOT.kBlue)
        h["parton"][var].SetLineWidth(2)
        
        h["shower"][var].SetLineColor(ROOT.kRed)
        h["shower"][var].SetLineWidth(2)
        
        max_y = max(h["parton"][var].GetMaximum(), h["shower"][var].GetMaximum())
        h["parton"][var].SetMaximum(max_y * 1.3)
        
        # Plot styling
        h["parton"][var].SetStats(0)
        h["shower"][var].SetStats(0)

        h["parton"][var].Draw("HIST")
        h["shower"][var].Draw("HIST SAME")
        
        legend = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
        legend.AddEntry(h["parton"][var], f"{args.process}: Parton Level", "l")
        legend.AddEntry(h["shower"][var], f"{args.process}: Delphes Reco", "l")
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.Draw()
        
        output_file = os.path.join(output_dir, f"{var}_{args.process}.png")
        canvas.SaveAs(output_file)
        canvas.Clear()

    # Save all histograms to a ROOT file
    root_output_file = os.path.join(output_dir, f"{args.process}_histograms.root")
    out_root = ROOT.TFile.Open(root_output_file, "RECREATE")
    for level in ["parton", "shower"]:
        for hist in h[level].values():
            hist.Write()
    out_root.Close()
    print(f"Saved ROOT file with all histograms: {root_output_file}")

    print(f"Analysis complete! Plots for {args.process} saved in '{output_dir}'.")

if __name__ == "__main__":
    main()