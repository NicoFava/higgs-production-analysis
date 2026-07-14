#!/usr/bin/env python3
import argparse
import math
import os
import sys

try:
    import ROOT
except ImportError:
    print('ERROR: PyROOT is required.')
    sys.exit(1)

ROOT.gROOT.SetBatch(True)
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

def main():
    parser = argparse.ArgumentParser(description="Compare 4 processes at Shower Level (Delphes).")
    parser.add_argument('--ggf', required=True, help="Path to ggF Delphes ROOT file")
    parser.add_argument('--vbf', required=True, help="Path to VBF Delphes ROOT file")
    parser.add_argument('--wh',  required=True, help="Path to WH Delphes ROOT file")
    parser.add_argument('--zh',  required=True, help="Path to ZH Delphes ROOT file")
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

        chain = ROOT.TChain("Delphes")
        chain.Add(proc_info["file"])
        treeReader = ROOT.ExRootTreeReader(chain)
        numberOfEntries = treeReader.GetEntries()

        branchEvent     = treeReader.UseBranch("Event")
        branchParticle  = treeReader.UseBranch("Particle")
        branchJet       = treeReader.UseBranch("Jet")
        branchElectron  = treeReader.UseBranch("Electron")
        branchMuon      = treeReader.UseBranch("Muon")
        branchMissingET = treeReader.UseBranch("MissingET")

        for entry in range(numberOfEntries):
            treeReader.ReadEntry(entry)
            if entry % 1000 == 0 and entry > 0:
                print(f"  ... {entry} events")

            event_weight = 1.0
            if branchEvent.GetEntries() > 0:
                event_weight = branchEvent.At(0).Weight

            higgs = None
            shower_jets = []
            leptons = []
            met_vec = ROOT.TLorentzVector()

            for i in range(branchParticle.GetEntries()):
                p = branchParticle.At(i)
                if abs(p.PID) == 25:
                    higgs = ROOT.TLorentzVector()
                    higgs.SetPtEtaPhiM(p.PT, p.Eta, p.Phi, p.Mass)

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

            for i in range(branchJet.GetEntries()):
                j = branchJet.At(i)
                if j.PT > 25.0 and abs(j.Eta) < 4.5:
                    vec = ROOT.TLorentzVector()
                    vec.SetPtEtaPhiM(j.PT, j.Eta, j.Phi, j.Mass)
                    shower_jets.append(vec)

            if branchMissingET.GetEntries() > 0:
                met = branchMissingET.At(0)
                met_vec.SetPtEtaPhiM(met.MET, 0.0, met.Phi, 0.0)

            shower_jets.sort(key=lambda j: j.Pt(), reverse=True)
            leptons.sort(key=lambda l: l.Pt(), reverse=True)

            if higgs:
                hists["pT_H"][proc_name].Fill(higgs.Pt(), event_weight)
                hists["eta_H"][proc_name].Fill(higgs.Eta(), event_weight)
            hists["Njets"][proc_name].Fill(len(shower_jets), event_weight)
            hists["HT"][proc_name].Fill(sum([j.Pt() for j in shower_jets]), event_weight)

            if len(shower_jets) >= 2:
                j1, j2 = shower_jets[0], shower_jets[1]
                hists["m_jj"][proc_name].Fill((j1 + j2).M(), event_weight)
                hists["dEta_jj"][proc_name].Fill(abs(j1.Eta() - j2.Eta()), event_weight)
                hists["dPhi_jj"][proc_name].Fill(abs(j1.DeltaPhi(j2)), event_weight)

            hists["N_l"][proc_name].Fill(len(leptons), event_weight)
            hists["Emiss_T"][proc_name].Fill(met_vec.Pt(), event_weight)
            if len(leptons) >= 2:
                hists["m_ll"][proc_name].Fill((leptons[0] + leptons[1]).M(), event_weight)
            if len(leptons) >= 1 and met_vec.Pt() > 0:
                dphi = abs(leptons[0].DeltaPhi(met_vec))
                mt = math.sqrt(max(0.0, 2 * leptons[0].Pt() * met_vec.Pt() * (1.0 - math.cos(dphi))))
                hists["m_WT"][proc_name].Fill(mt, event_weight)

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