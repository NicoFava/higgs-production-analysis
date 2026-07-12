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

try:
    import numpy as np
    from pyjet import cluster
except ImportError as e:
    print('ERROR: numpy and pyjet are required for HepMC jet clustering.')
    print('Please run: pip install numpy pyjet')
    sys.exit(1)

# Disable graphic output on screen for faster batch processing
ROOT.gROOT.SetBatch(True)
# Ensure correct error tracking when scaling histograms
ROOT.TH1.SetDefaultSumw2(True)

# ==============================================================================
# PARSING FUNCTIONS FOR LHE AND HEPMC
# ==============================================================================
def read_lhe_events(filepath):
    """Yields events from a zipped LHE file. Each event is a list of particles."""
    with gzip.open(filepath, 'rt') as f:
        in_event = False
        is_first_line = False
        event_particles = []
        for line in f:
            line = line.strip()
            if not line: continue
            
            if line.startswith("<event>"):
                in_event = True
                is_first_line = True
                event_particles = []
                continue
                
            if line.startswith("</event>"):
                in_event = False
                yield event_particles
                continue
                
            if in_event:
                if is_first_line:
                    is_first_line = False
                    continue  # Skip event weight/scale line
                
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

def read_hepmc_events(filepath):
    """Yields events from a zipped HepMC file. Each event is a list of particles."""
    with gzip.open(filepath, 'rt') as f:
        event_particles = []
        for line in f:
            line = line.strip()
            if line.startswith('E '):
                if event_particles:
                    yield event_particles
                event_particles = []
            elif line.startswith('P '):
                parts = line.split()
                if len(parts) >= 8:
                    pdg = int(parts[2])
                    px, py, pz, e, mass = map(float, parts[3:8])
                    status = int(parts[8])
                    event_particles.append({
                        'pdg': pdg, 'status': status,
                        'px': px, 'py': py, 'pz': pz, 'e': e, 'm': mass
                    })
        if event_particles:
            yield event_particles

# ==============================================================================
# MAIN ANALYSIS LOOP
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(description="Analyze Higgs production from .lhe and .hepmc directly.")
    parser.add_argument('--lhe', type=str, required=True, help="Path to the input .lhe.gz file (Parton Level)")
    parser.add_argument('--hepmc', type=str, required=True, help="Path to the input .hepmc.gz file (Shower Level)")
    parser.add_argument('-o', '--outdir', type=str, required=True, help="Output directory for the generated plots")
    parser.add_argument('-p', '--process', type=str, required=True, help="Name of the physical process (e.g., ggF, VBF, WH, ZH)")
    args = parser.parse_args()

    if not os.path.exists(args.lhe):
        raise RuntimeError(f"Could not find LHE file: {args.lhe}")
    if not os.path.exists(args.hepmc):
        raise RuntimeError(f"Could not find HepMC file: {args.hepmc}")

    # --------------------------------------------------------------------------
    # 1. HISTOGRAM DEFINITIONS
    # --------------------------------------------------------------------------
    h = {"parton": {}, "shower": {}}

    for level in ["parton", "shower"]:
        prefix = f"{level}_{args.process}"
        
        # General observables
        h[level]["Njets"] = ROOT.TH1F(f"Njets_{prefix}", f"{args.process}: Number of Jets; N_{{jets}}; Events", 10, -0.5, 9.5)
        h[level]["pT_H"]  = ROOT.TH1F(f"pT_H_{prefix}", f"{args.process}: Higgs Boson p_{{T}}; p_{{T}}^{{H}} [GeV]; Events", 50, 0, 400)
        h[level]["HT"]    = ROOT.TH1F(f"HT_{prefix}", f"{args.process}: Scalar Sum of Jet p_{{T}}; H_{{T}} [GeV]; Events", 50, 0, 800)
        
        # VBF Specific Observables (Requires >= 2 jets)
        h[level]["m_jj"]    = ROOT.TH1F(f"m_jj_{prefix}", f"{args.process}: Dijet Invariant Mass; m_{{jj}} [GeV]; Events", 50, 0, 2500)
        h[level]["dEta_jj"] = ROOT.TH1F(f"dEta_jj_{prefix}", f"{args.process}: Pseudorapidity Separation; |#Delta#eta_{{jj}}|; Events", 50, 0, 8)
        h[level]["dPhi_jj"] = ROOT.TH1F(f"dPhi_jj_{prefix}", f"{args.process}: Azimuthal Separation; #Delta#phi_{{jj}} [rad]; Events", 50, 0, ROOT.TMath.Pi())

        # VH Specific Observables (Requires Leptons / MET)
        h[level]["N_l"]     = ROOT.TH1F(f"N_l_{prefix}", f"{args.process}: Number of Leptons (N_{{f}}); N_{{#ell}}; Events", 6, -0.5, 5.5)
        h[level]["Emiss_T"] = ROOT.TH1F(f"Emiss_T_{prefix}", f"{args.process}: Missing Transverse Energy; E_{{T}}^{{miss}} [GeV]; Events", 50, 0, 400)
        h[level]["m_ll"]    = ROOT.TH1F(f"m_ll_{prefix}", f"{args.process}: Dilepton Invariant Mass (m_{{ee}}); m_{{#ell#ell}} [GeV]; Events", 50, 0, 200)
        h[level]["m_WT"]    = ROOT.TH1F(f"m_WT_{prefix}", f"{args.process}: W Transverse Mass; m_{{T}}^{{W}} [GeV]; Events", 50, 0, 200)

    # Helper function to avoid repetition
    def fill_histograms(level, higgs, jets, leptons, met_vec):
        if higgs:
            h[level]["pT_H"].Fill(higgs.Pt())
        
        h[level]["Njets"].Fill(len(jets))
        h[level]["HT"].Fill(sum([j.Pt() for j in jets]))
        
        if len(jets) >= 2:
            j1, j2 = jets[0], jets[1]
            h[level]["m_jj"].Fill((j1 + j2).M())
            h[level]["dEta_jj"].Fill(abs(j1.Eta() - j2.Eta()))
            h[level]["dPhi_jj"].Fill(abs(j1.DeltaPhi(j2)))

        h[level]["N_l"].Fill(len(leptons))
        h[level]["Emiss_T"].Fill(met_vec.Pt())
        
        if len(leptons) >= 2:
            h[level]["m_ll"].Fill((leptons[0] + leptons[1]).M())
            
        if len(leptons) >= 1 and met_vec.Pt() > 0:
            l1 = leptons[0]
            dphi = abs(l1.DeltaPhi(met_vec))
            mt = math.sqrt(max(0.0, 2 * l1.Pt() * met_vec.Pt() * (1.0 - math.cos(dphi))))
            h[level]["m_WT"].Fill(mt)

    # --------------------------------------------------------------------------
    # 2. PROCESS LHE FILE (PARTON LEVEL)
    # --------------------------------------------------------------------------
    print("Processing LHE events (Parton Level)...")
    for ievt, event in enumerate(read_lhe_events(args.lhe)):
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

        fill_histograms("parton", higgs_parton, jets, leptons, met_vec)

    # --------------------------------------------------------------------------
    # 3. PROCESS HEPMC FILE (SHOWER LEVEL)
    # --------------------------------------------------------------------------
    print("Processing HepMC events (Shower Level)...")
    for ievt, event in enumerate(read_hepmc_events(args.hepmc)):
        if ievt % 1000 == 0 and ievt > 0:
            print(f"  Processed {ievt} HepMC events...")

        higgs_shower = None
        leptons = []
        met_x, met_y = 0.0, 0.0
        
        visible_particles = []

        for p in event:
            pid = abs(p['pdg'])
            
            # Track the last copy of the Higgs (usually post-ISR radiation)
            if pid == 25:
                higgs_shower = ROOT.TLorentzVector()
                higgs_shower.SetPxPyPzE(p['px'], p['py'], p['pz'], p['e'])
            
            # Status 1 = final state particles after shower & hadronization
            if p['status'] == 1:
                vec = ROOT.TLorentzVector()
                vec.SetPxPyPzE(p['px'], p['py'], p['pz'], p['e'])
                
                # Leptons
                if pid in [11, 13] and vec.Pt() > 10.0 and abs(vec.Eta()) < 2.5:
                    leptons.append(vec)
                # Neutrinos (MET)
                elif pid in [12, 14, 16]:
                    met_x += p['px']
                    met_y += p['py']
                
                # Collect all visible particles for Jet Clustering (exclude neutrinos/prompt leptons)
                if pid not in [11, 12, 13, 14, 16, 25]:
                    vec = ROOT.TLorentzVector()
                    vec.SetPxPyPzE(p['px'], p['py'], p['pz'], p['e'])
                    # Require a tiny minimum pT to avoid infinite pseudo-rapidity for beamline remnants
                    if vec.Pt() > 0.01:
                        visible_particles.append((vec.Pt(), vec.Eta(), vec.Phi(), vec.M()))

        # Jet Clustering using PyJet (Anti-kT, R=0.4)
        shower_jets = []
        if visible_particles:
            # Use the default kinematic variables that pyjet expects
            vectors = np.array(visible_particles, dtype=np.dtype([
                ('pT', 'f8'), ('eta', 'f8'), ('phi', 'f8'), ('mass', 'f8')
            ]))
            sequence = cluster(vectors, R=0.4, p=-1) # p=-1 is Anti-kT
            for jet in sequence.inclusive_jets(ptmin=25.0):
                if abs(jet.eta) < 4.5:
                    vec = ROOT.TLorentzVector()
                    vec.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass)
                    shower_jets.append(vec)

        # Jet Clustering using PyJet (Anti-kT, R=0.4)
        shower_jets = []
        if visible_particles:
            vectors = np.array(visible_particles, dtype=np.dtype([('px', 'f8'), ('py', 'f8'), ('pz', 'f8'), ('E', 'f8')]))
            sequence = cluster(vectors, R=0.4, p=-1) # p=-1 is Anti-kT
            for jet in sequence.inclusive_jets(ptmin=25.0):
                if abs(jet.eta) < 4.5:
                    vec = ROOT.TLorentzVector()
                    vec.SetPxPyPzE(jet.px, jet.py, jet.pz, jet.e)
                    shower_jets.append(vec)

        shower_jets.sort(key=lambda j: j.Pt(), reverse=True)
        leptons.sort(key=lambda l: l.Pt(), reverse=True)
        met_vec = ROOT.TLorentzVector()
        met_vec.SetPxPyPzE(met_x, met_y, 0.0, math.sqrt(met_x**2 + met_y**2))

        fill_histograms("shower", higgs_shower, shower_jets, leptons, met_vec)

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
        canvas.Clear()

    print(f"Analysis complete! Plots for {args.process} saved in '{output_dir}'.")

if __name__ == "__main__":
    main()