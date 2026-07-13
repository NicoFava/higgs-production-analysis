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

# Global ROOT settings for a clean visual output
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.gStyle.SetOptStat(0) # Remove statistics box

# Exact list of the 11 kinematic variables
VARIABLES = [
    "Njets", "pT_H", "eta_H", "HT", 
    "m_jj", "dEta_jj", "dPhi_jj", 
    "N_l", "Emiss_T", "m_ll", "m_WT"
]

def get_histogram(tfile, hist_name, fallback_path=None):
    """Fetches a histogram from the ROOT file. Tries direct path, then fallback."""
    hist = tfile.Get(hist_name)
    if not hist and fallback_path:
        hist = tfile.Get(fallback_path)
    return hist

def normalize_hist(hist):
    """Normalizes the histogram to Area = 1 for shape comparison."""
    if not hist.GetSumw2N():
        hist.Sumw2()
    integral = hist.Integral()
    if integral > 0:
        hist.Scale(1.0 / integral)
    return hist

def format_data_hist(hist):
    """Style for ATLAS data: Black markers with error bars."""
    hist.SetLineColor(ROOT.kBlack)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetMarkerStyle(20) # Full Circle
    hist.SetMarkerSize(1.0)
    hist.SetLineWidth(2)

def format_mc_hist(hist, color, style):
    """Style for Monte Carlo: Solid or dashed lines."""
    hist.SetLineColor(color)
    hist.SetLineWidth(3)
    hist.SetLineStyle(style)
    hist.SetFillStyle(0) # No fill to prevent covering other curves

def main():
    parser = argparse.ArgumentParser(description="Overlay ATLAS Data vs MadGraph (Parton & Shower).")
    parser.add_argument("--data", nargs='+', required=True, help="List of ATLAS blind sample ROOT files")
    parser.add_argument("--mc", nargs='+', required=True, help="List of MadGraph ROOT files")
    parser.add_argument("-o", "--outdir", default="results/comparisons", help="Base output directory")
    args = parser.parse_args()

    # Create the canvas exactly ONCE outside the loops
    canvas = ROOT.TCanvas("c_compare", "Kinematics Comparison", 800, 600)

    # 1. Loop over each ATLAS Data file
    for data_path in args.data:
        if not os.path.exists(data_path):
            print(f"[!] Error: Data file {data_path} not found. Skipping.")
            continue
            
        f_data = ROOT.TFile.Open(data_path, "READ")
        data_name = os.path.basename(data_path).replace("_kinematics.root", "")
        print(f"\n[*] Processing ATLAS Data Sample: {data_name}")

        # 2. Loop over each Monte Carlo process
        for mc_path in args.mc:
            if not os.path.exists(mc_path):
                print(f"  [!] Warning: MC file {mc_path} not found. Skipping.")
                continue
                
            f_mc = ROOT.TFile.Open(mc_path, "READ")
            process_name = os.path.basename(mc_path).split('_')[0]
            print(f"  -> Comparing with MadGraph Process: {process_name}")

            # Create specific sub-directory for this Data vs MC comparison
            process_outdir = os.path.join(args.outdir, data_name, process_name)
            os.makedirs(process_outdir, exist_ok=True)

            # 3. Loop over the 11 variables
            for var in VARIABLES:
                # Fetch Histograms
                h_data = get_histogram(f_data, f"{var}_{data_name}")
                h_parton = get_histogram(f_mc, f"{var}_parton_{process_name}", f"parton/{var}_parton_{process_name}")
                h_shower = get_histogram(f_mc, f"{var}_shower_{process_name}", f"shower/{var}_shower_{process_name}")

                if not h_data or not h_parton or not h_shower:
                    continue # Skip this variable if any of the 3 histograms are missing

                # Clone to avoid modifying original objects in memory
                hd = h_data.Clone(f"plot_{var}_data_{data_name}_{process_name}")
                hp = h_parton.Clone(f"plot_{var}_parton_{data_name}_{process_name}")
                hs = h_shower.Clone(f"plot_{var}_shower_{data_name}_{process_name}")

                # Detach clones from ROOT directory management
                hd.SetDirectory(0)
                hp.SetDirectory(0)
                hs.SetDirectory(0)

                # Normalize all to Area = 1
                hd = normalize_hist(hd)
                hp = normalize_hist(hp)
                hs = normalize_hist(hs)

                # Apply visual styles
                format_data_hist(hd)
                format_mc_hist(hp, ROOT.kBlue, 2)  # Parton: Dashed Blue
                format_mc_hist(hs, ROOT.kRed, 1)   # Shower: Solid Red

                # DYNAMIC Y-AXIS MAXIMUM CALCULATION (Safe check using Integrals)
                max_y = 0.0
                if hd.Integral() > 0:
                    for b in range(1, hd.GetNbinsX() + 1):
                        data_upper_edge = hd.GetBinContent(b) + hd.GetBinError(b)
                        if data_upper_edge > max_y:
                            max_y = data_upper_edge
                        
                if hp.Integral() > 0:
                    max_y = max(max_y, hp.GetMaximum())
                if hs.Integral() > 0:
                    max_y = max(max_y, hs.GetMaximum())
                
                # Fallback value if every single histogram is completely empty
                if max_y <= 0:
                    max_y = 1.0
                
                # --- THE ULTIMATE FIX: USE A DUMMY FRAME ---
                # Create a completely empty histogram to act ONLY as the canvas axis frame.
                # This decouples the drawing of the axes from the actual data distributions.
                h_frame = hs.Clone(f"frame_{var}_{data_name}_{process_name}")
                h_frame.Reset() # Clears all data and errors, leaving only the binning structure
                h_frame.SetDirectory(0)

                # Configure the dummy frame with our dynamically calculated maximum
                h_frame.SetMaximum(max_y * 1.35) 
                h_frame.SetMinimum(0.0)
                h_frame.SetTitle(f"Shape Comparison: {var} ({data_name} vs {process_name})")
                h_frame.GetYaxis().SetTitle("Normalized Units")
                h_frame.GetYaxis().SetTitleOffset(1.4)

                # Clear the canvas and reset margins
                canvas.Clear()
                ROOT.gPad.SetLeftMargin(0.15)
                ROOT.gPad.SetBottomMargin(0.12)

                # 1. Draw ONLY the axes framework first
                h_frame.Draw("AXIS")

                # 2. Overlay the actual histograms on top safely.
                # If any of these are mathematically empty, ROOT will gracefully draw nothing
                # for them without breaking the axis scales we established above.
                hs.Draw("HIST SAME")
                hp.Draw("HIST SAME")
                hd.Draw("PE SAME")

                # Create Legend
                leg = ROOT.TLegend(0.50, 0.70, 0.88, 0.88)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)
                leg.AddEntry(hd, f"ATLAS {data_name}", "pe")
                leg.AddEntry(hp, f"MadGraph {process_name} (Parton)", "l")
                leg.AddEntry(hs, f"MadGraph {process_name} (Shower)", "l")
                leg.Draw()

                # Save the plot
                out_plot_path = os.path.join(process_outdir, f"compare_{var}_{data_name}_vs_{process_name}.png")
                canvas.SaveAs(out_plot_path)

                # Free memory completely
                h_frame.Delete()
                hd.Delete()
                hp.Delete()
                hs.Delete()

            f_mc.Close()

        f_data.Close()

    print("\n[+] Comparison completed successfully! All plots saved in:", args.outdir)

if __name__ == "__main__":
    main()