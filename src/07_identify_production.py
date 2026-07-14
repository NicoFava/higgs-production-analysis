#!/usr/bin/env python3
"""
Script for Higgs Production Mechanism Identification (Signal Region Approach).

This script identifies the production mechanism by defining a signal region
(value window) for a single highly discriminating variable for each process.
"""

import argparse
import os
import sys
import array

try:
    import ROOT
except ImportError as e:
    print('ERROR: PyROOT is required to run this script.')
    print(e)
    sys.exit(1)

# ROOT global settings
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

# ==============================================================================
# DEFINITION OF SIGNAL REGIONS
# ==============================================================================
# For each process define the golden variable and the range (min, max) that
# identifies the signal region.
# NOTE: Tune these values by looking at your histograms!
SIGNAL_REGIONS = {
    "VBF": ("m_jj", 400.0, 1100.0),    # Remove the empty tail to avoid Chi2 warnings
    "WH":  ("m_WT", 0.0, 120.0),      # W peak (already good)
    "ZH":  ("m_ll", 80.0, 100.0),      # Use the associated Z mass peak!
    "ggF": ("pT_H", 0.0, 100.0)        # Low pT (already good)
}

def get_histogram_from_file(tfile, var_name, sample_name, is_mc=False):
    """
    Retrieve the histogram from the ROOT file using naming patterns.
    """
    patterns = []
    if is_mc:
        patterns.append(f"{var_name}_shower_{sample_name}")
        patterns.append(f"{var_name}_{sample_name}_shower")
        patterns.append(f"{var_name}_shower")
    else:
        patterns.append(f"{var_name}_{sample_name}")
        patterns.append(f"{var_name}_{sample_name}_kinematics")
    
    patterns.append(var_name)
    
    for pattern in patterns:
        hist = tfile.Get(pattern)
        if hist:
            return hist
            
    return None

def extract_signal_region(hist, xmin, xmax, name_suffix):
    """
    Create a new TH1D that contains only the bins in the signal region [xmin, xmax].
    Returns the new histogram and the fraction of events in this region.
    """
    if not hist:
        return None, 0.0
        
    xaxis = hist.GetXaxis()
    bin_min = xaxis.FindBin(xmin)
    bin_max = xaxis.FindBin(xmax)
    
    # Se il max va oltre l'ultimo bin, fermiamoci alla fine dell'asse
    if bin_max > hist.GetNbinsX():
        bin_max = hist.GetNbinsX()
        
    n_bins = bin_max - bin_min + 1
    
    # Raccogliamo i bordi dei bin per ricreare fedelmente l'asse X
    edges = []
    for i in range(bin_min, bin_max + 2):
        edges.append(xaxis.GetBinLowEdge(i))
        
    new_name = f"{hist.GetName()}_{name_suffix}"
    h_sr = ROOT.TH1D(new_name, hist.GetTitle(), n_bins, array.array('d', edges))
    h_sr.Sumw2()
    
    # Copy contents and errors
    for i in range(1, n_bins + 1):
        orig_bin = bin_min + i - 1
        h_sr.SetBinContent(i, hist.GetBinContent(orig_bin))
        h_sr.SetBinError(i, hist.GetBinError(orig_bin))
        
    # Calculate the event fraction in the SR
    total_integral = hist.Integral()
    sr_integral = h_sr.Integral()
    fraction = (sr_integral / total_integral) if total_integral > 0 else 0.0
    
    return h_sr, fraction

def normalize_hist(hist):
    """Normalize the histogram to unit area for shape-only comparison."""
    if not hist:
        return None
    integral = hist.Integral()
    if integral > 0:
        hist.Scale(1.0 / integral)
    return hist

def main():
    parser = argparse.ArgumentParser(description="Signal region Chi-square test for Higgs production.")
    parser.add_argument("--data", nargs='+', required=True, help="List of processed blind sample ROOT files")
    parser.add_argument("--mc", nargs='+', required=True, help="List of MC ROOT files containing shower histograms")
    args = parser.parse_args()

    # Struttura dei risultati: { blind_sample: { mc_process: (var, chi2_ndf, p_val, frac_data, frac_mc) } }
    results = {}

    for data_path in args.data:
        if not os.path.exists(data_path):
            continue
            
        f_data = ROOT.TFile.Open(data_path, "READ")
        data_file_base = os.path.basename(data_path)
        blind_name = data_file_base.replace(".root", "").replace("_kinematics", "")
        results[blind_name] = {}

        for mc_path in args.mc:
            if not os.path.exists(mc_path):
                continue
                
            f_mc = ROOT.TFile.Open(mc_path, "READ")
            mc_process = os.path.basename(mc_path).split('_')[0]
            
            if mc_process not in SIGNAL_REGIONS:
                f_mc.Close()
                continue
                
            var, xmin, xmax = SIGNAL_REGIONS[mc_process]

            # 1. Retrieve full histograms
            hist_data = get_histogram_from_file(f_data, var, blind_name, is_mc=False)
            hist_mc = get_histogram_from_file(f_mc, var, mc_process, is_mc=True)

            if not hist_data or not hist_mc:
                print(f" [!] Error: Variable '{var}' missing for {blind_name} or {mc_process}")
                results[blind_name][mc_process] = (var, xmin, xmax, 999.9, 0.0, 0.0, 0.0)
                f_mc.Close()
                continue

            # 2. Extract the signal region (and event fraction)
            h_data_sr, frac_data = extract_signal_region(hist_data, xmin, xmax, "data_sr")
            h_mc_sr, frac_mc = extract_signal_region(hist_mc, xmin, xmax, "mc_sr")

            # --- ADD THESE TWO LINES ---
            h_data_sr.Rebin(2)
            h_mc_sr.Rebin(2)
            # ---------------------------------

            # 3. Normalize only the SR to unit area (for pure shape testing)
            if h_data_sr.Integral() == 0 or h_mc_sr.Integral() == 0:
                chi2_ndf = 999.9
                p_val = 0.0
            else:
                normalize_hist(h_data_sr)
                normalize_hist(h_mc_sr)
                
                # 4. Chi-square calculation on the SR
                p_val = h_data_sr.Chi2Test(h_mc_sr, "WW")
                chi2_ndf = h_data_sr.Chi2Test(h_mc_sr, "WW CHI2/NDF")

            results[blind_name][mc_process] = (var, xmin, xmax, chi2_ndf, p_val, frac_data, frac_mc)

            # Cleanup
            h_data_sr.Delete()
            h_mc_sr.Delete()
            f_mc.Close()
            
        f_data.Close()

    # ==============================================================================
    # PRINT RESULTS
    # ==============================================================================
    output_lines = []
    output_lines.append("\n" + "="*95)
    output_lines.append("      HIGGS PRODUCTION IDENTIFICATION MATRIX (SIGNAL REGION APPROACH)")
    output_lines.append("="*95)

    for blind_name in sorted(results.keys()):
        output_lines.append(f"\n>>> TEST RESULTS FOR BLIND SAMPLE: {blind_name}")
        output_lines.append("-" * 95)
        
        mc_hypotheses = results[blind_name]
        for mc_process in sorted(mc_hypotheses.keys()):
            var, xmin, xmax, chi2, pval, f_data, f_mc = mc_hypotheses[mc_process]
            
            # Formatting to make results clear
            output_lines.append(f" -> HYPOTHESIS: {mc_process:<4} | Key variable: {var} (in [{xmin}, {xmax}])")
            
            # Event fraction: this is very useful!
            f_data_str = f"{f_data*100:.1f}%"
            f_mc_str = f"{f_mc*100:.1f}%"
            
            # Color if the shape is compact
            chi2_str = f"{chi2:.2f}"
            if chi2 < 2.5 and chi2 > 0.0:
                chi2_display = f"{chi2_str:<6} (Shape Match)"
            else:
                chi2_display = f"{chi2_str:<19}"

            output_lines.append(f"    Events in SR: Data = {f_data_str:<7} | MC expectation = {f_mc_str:<7}")
            output_lines.append(f"    Shape test : Chi2/NDF = {chi2_display} | p-value = {pval:.4e}")
            output_lines.append("    " + "-"*90)
            
        output_lines.append("="*95)

    results_dir = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
    os.makedirs(results_dir, exist_ok=True)
    output_path = os.path.join(results_dir, "output.txt")

    with open(output_path, "w") as outfile:
        outfile.write("\n".join(output_lines))

    for line in output_lines:
        print(line)

    print(f"\nDetailed results also written to: {output_path}")

if __name__ == "__main__":
    main()