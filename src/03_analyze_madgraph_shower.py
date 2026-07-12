#!/usr/bin/env python3
import argparse
import math
import os
import sys
import gzip

try:
    import ROOT
except ImportError:
    print('ERROR: PyROOT is required.')
    sys.exit(1)
try:
    import numpy as np
    from pyjet import cluster
except ImportError:
    print('ERROR: numpy and pyjet are required.')
    sys.exit(1)

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

def read_hepmc_events(filepath):
    with gzip.open(filepath, 'rt') as f:
        event_particles = []
        for line in f:
            line = line.strip()
            if line.startswith('E '):
                if event_particles: yield event_particles
                event_particles = []
            elif line.startswith('P '):
                parts = line.split()
                if len(parts) >= 8:
                    event_particles.append({
                        'pdg': int(parts[2]), 'status': int(parts[8]),
                        'px': float(parts[3]), 'py': float(parts[4]), 'pz': float(parts[5]), 
                        'e': float(parts[6]), 'm': float(parts[7])
                    })
        if event_particles: yield event_particles

def main():
    parser = argparse.ArgumentParser(description="Compare 4 processes at Shower Level (HepMC).")
    parser.add_argument('--ggf', required=True, help="Path to ggF .hepmc.gz")
    parser.add_argument('--vbf', required=True, help="Path to VBF .hepmc.gz")
    parser.add_argument('--wh',  required=True, help="Path to WH .hepmc.gz")
    parser.add_argument('--zh',  required=True, help="Path to ZH .hepmc.gz")
    parser.add_argument('-o', '--outdir', required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    processes = {
        "ggF": {"file": args.ggf, "color": ROOT.kBlue},
        "VBF": {"file": args.vbf, "color": ROOT.kRed},
        "WH":  {"file": args.wh,  "color": ROOT.kGreen+2},
        "ZH":  {"file": args.zh,  "color": ROOT.kMagenta}
    }

    hist_configs = {
        "Njets":   {"title": "Number of Jets; N_{jets}; Normalized Events", "b": 10, "xmin": -0.5, "xmax": 9.5},
        "pT_H":    {"title": "Higgs Boson p_{T}; p_{T}^{H} [GeV]; Normalized Events", "b": 50, "xmin": 0, "xmax": 400},
        "eta_H":   {"title": "Higgs Boson #eta; #eta_{H}; Normalized Events", "b": 50, "xmin": -5, "xmax": 5},
        "HT":      {"title": "Scalar Sum of Jet p_{T}; H_{T} [GeV]; Normalized Events", "b": 50, "xmin": 0, "xmax": 800},
        "m_jj":    {"title": "Dijet Invariant Mass; m_{jj} [GeV]; Normalized Events", "b": 50, "xmin": 0, "xmax": 2500},
        "dEta_jj": {"title": "Pseudorapidity Separation; |#Delta#eta_{jj}|; Normalized Events", "b": 50, "xmin": 0, "xmax": 8},
        "dPhi_jj": {"title": "Azimuthal Separation; #Delta#phi_{jj} [rad]; Normalized Events", "b": 50, "xmin": 0, "xmax": ROOT.TMath.Pi()},
        "N_l":     {"title": "Number of Leptons; N_{l}; Normalized Events", "b": 6, "xmin": -0.5, "xmax": 5.5},
        "Emiss_T": {"title": "Missing Transverse Energy; E_{T}^{miss} [GeV]; Normalized Events", "b": 50, "xmin": 0, "xmax": 400},
        "m_ll":    {"title": "Dilepton Invariant Mass; m_{ll} [GeV]; Normalized Events", "b": 50, "xmin": 0, "xmax": 200},
        "m_WT":    {"title": "W Transverse Mass; m_{T}^{W} [GeV]; Normalized Events", "b": 50, "xmin": 0, "xmax": 200}
    }

    hists = {var: {} for var in hist_configs.keys()}

    for proc_name, proc_info in processes.items():
        print(f"Processing Shower Level for: {proc_name}")
        for var, cfg in hist_configs.items():
            hists[var][proc_name] = ROOT.TH1F(f"{var}_{proc_name}_shower", cfg["title"], cfg["b"], cfg["xmin"], cfg["xmax"])
            hists[var][proc_name].SetLineColor(proc_info["color"])
            hists[var][proc_name].SetLineWidth(2)
            hists[var][proc_name].SetStats(0)

        if not os.path.exists(proc_info["file"]):
            print(f"  WARNING: File {proc_info['file']} not found.")
            continue

        for ievt, event in enumerate(read_hepmc_events(proc_info["file"])):
            if ievt % 1000 == 0 and ievt > 0: print(f"  ... {ievt} events")
            
            higgs = None
            leptons = []
            met_x, met_y = 0.0, 0.0
            visible_particles = []

            for p in event:
                pid = abs(p['pdg'])
                if pid == 25:
                    higgs = ROOT.TLorentzVector()
                    higgs.SetPxPyPzE(p['px'], p['py'], p['pz'], p['e'])
                
                if p['status'] == 1:
                    vec = ROOT.TLorentzVector()
                    vec.SetPxPyPzE(p['px'], p['py'], p['pz'], p['e'])
                    if pid in [11, 13] and vec.Pt() > 10.0 and abs(vec.Eta()) < 2.5:
                        leptons.append(vec)
                    elif pid in [12, 14, 16]:
                        met_x += p['px']
                        met_y += p['py']
                    
                    if pid not in [11, 12, 13, 14, 16, 25] and vec.Pt() > 0.01:
                        visible_particles.append((vec.Pt(), vec.Eta(), vec.Phi(), vec.M()))

            shower_jets = []
            if visible_particles:
                vectors = np.array(visible_particles, dtype=np.dtype([
                    ('pT', 'f8'), ('eta', 'f8'), ('phi', 'f8'), ('mass', 'f8')
                ]))
                sequence = cluster(vectors, R=0.4, p=-1)
                for jet in sequence.inclusive_jets(ptmin=25.0):
                    if abs(jet.eta) < 4.5:
                        vec = ROOT.TLorentzVector()
                        vec.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass)
                        shower_jets.append(vec)

            shower_jets.sort(key=lambda j: j.Pt(), reverse=True)
            leptons.sort(key=lambda l: l.Pt(), reverse=True)
            met_vec = ROOT.TLorentzVector()
            met_vec.SetPxPyPzE(met_x, met_y, 0.0, math.sqrt(met_x**2 + met_y**2))

            if higgs:
                hists["pT_H"][proc_name].Fill(higgs.Pt())
                hists["eta_H"][proc_name].Fill(higgs.Eta())
            hists["Njets"][proc_name].Fill(len(shower_jets))
            hists["HT"][proc_name].Fill(sum([j.Pt() for j in shower_jets]))
            
            if len(shower_jets) >= 2:
                j1, j2 = shower_jets[0], shower_jets[1]
                hists["m_jj"][proc_name].Fill((j1 + j2).M())
                hists["dEta_jj"][proc_name].Fill(abs(j1.Eta() - j2.Eta()))
                hists["dPhi_jj"][proc_name].Fill(abs(j1.DeltaPhi(j2)))

            hists["N_l"][proc_name].Fill(len(leptons))
            hists["Emiss_T"][proc_name].Fill(met_vec.Pt())
            if len(leptons) >= 2:
                hists["m_ll"][proc_name].Fill((leptons[0] + leptons[1]).M())
            if len(leptons) >= 1 and met_vec.Pt() > 0:
                dphi = abs(leptons[0].DeltaPhi(met_vec))
                mt = math.sqrt(max(0.0, 2 * leptons[0].Pt() * met_vec.Pt() * (1.0 - math.cos(dphi))))
                hists["m_WT"][proc_name].Fill(mt)

    canvas = ROOT.TCanvas("c", "Canvas", 800, 600)
    for var in hists.keys():
        max_y = 0.0
        for proc in processes.keys():
            if hists[var][proc].Integral() > 0:
                hists[var][proc].Scale(1.0 / hists[var][proc].Integral())
            if hists[var][proc].GetMaximum() > max_y:
                max_y = hists[var][proc].GetMaximum()

        legend = ROOT.TLegend(0.70, 0.70, 0.88, 0.88)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)

        first = True
        for proc in processes.keys():
            hists[var][proc].SetMaximum(max_y * 1.3)
            hists[var][proc].GetYaxis().SetTitle("Normalized Units")
            opt = "HIST" if first else "HIST SAME"
            hists[var][proc].Draw(opt)
            legend.AddEntry(hists[var][proc], proc, "l")
            first = False
        
        legend.Draw()
        canvas.SaveAs(os.path.join(args.outdir, f"{var}_shower_comparison.png"))
        canvas.Clear()

    print("Shower level comparison complete!")

if __name__ == "__main__":
    main()