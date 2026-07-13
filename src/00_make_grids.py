#!/usr/bin/env python3
import os
import argparse
import sys

# Require Pillow for fast image manipulation
try:
    from PIL import Image
except ImportError:
    print("ERROR: The 'Pillow' library is required to run this script.")
    print("Please install it by running: pip install Pillow")
    sys.exit(1)

# Exact list of the 11 kinematic variables
VARIABLES = [
    "Njets", "pT_H", "eta_H", "HT", 
    "m_jj", "dEta_jj", "dPhi_jj", 
    "N_l", "Emiss_T", "m_ll", "m_WT"
]

# The 4 production processes in the exact order you requested:
# Top-Left: ggF, Top-Right: VBF
# Bottom-Left: WH, Bottom-Right: ZH
PROCESSES = ["ggF", "VBF", "WH", "ZH"]

def main():
    parser = argparse.ArgumentParser(description="Combine 4 comparison plots into a 2x2 grid.")
    parser.add_argument("--data", nargs='+', required=True, help="List of ATLAS data names (e.g., H4l_A H4l_B)")
    parser.add_argument("-i", "--indir", default="results/comparisons", help="Input directory with the individual plots")
    parser.add_argument("-o", "--outdir", default="results/grids", help="Output directory for the 2x2 grids")
    args = parser.parse_args()

    # Create the base output directory
    os.makedirs(args.outdir, exist_ok=True)

    # 1. Loop over each Data sample (e.g., H4l_A)
    for data_name in args.data:
        # Clean up data_name just in case the user passed the file path instead of the string name
        data_name = os.path.basename(data_name).replace("_kinematics.root", "").replace(".root", "")
        print(f"\n[*] Creating grids for Data Sample: {data_name}")

        # Create a specific output folder for this data sample's grids
        sample_outdir = os.path.join(args.outdir, data_name)
        os.makedirs(sample_outdir, exist_ok=True)

        # 2. Loop over each of the 11 variables
        for var in VARIABLES:
            images = []
            
            # Fetch the 4 corresponding plots (ggF, VBF, WH, ZH)
            for process in PROCESSES:
                plot_path = os.path.join(args.indir, data_name, process, f"compare_{var}_{data_name}_vs_{process}.png")
                
                if os.path.exists(plot_path):
                    img = Image.open(plot_path)
                    images.append(img)
                else:
                    print(f"  [!] Warning: Missing plot {plot_path}. Grid for {var} might be incomplete.")
                    images.append(None) # Keep None to maintain the 2x2 positioning if one is missing
            
            # If all 4 images are completely missing, skip this variable
            if all(img is None for img in images):
                continue

            # Assume all generated plots have the same dimensions based on the first valid one
            valid_img = next(img for img in images if img is not None)
            width, height = valid_img.size

            # Create a blank white canvas big enough to hold 2x2 images
            grid_img = Image.new('RGB', (width * 2, height * 2), color='white')

            # Paste the images into the grid
            # Index 0: ggF (Top-Left)
            if images[0]: grid_img.paste(images[0], (0, 0))
            # Index 1: VBF (Top-Right)
            if images[1]: grid_img.paste(images[1], (width, 0))
            # Index 2: WH (Bottom-Left)
            if images[2]: grid_img.paste(images[2], (0, height))
            # Index 3: ZH (Bottom-Right)
            if images[3]: grid_img.paste(images[3], (width, height))

            # Save the final grid
            out_grid_path = os.path.join(sample_outdir, f"grid_{var}_{data_name}.png")
            grid_img.save(out_grid_path)
            print(f"  -> Saved {out_grid_path}")

    print("\n[+] All 2x2 grids generated successfully!")

if __name__ == "__main__":
    main()