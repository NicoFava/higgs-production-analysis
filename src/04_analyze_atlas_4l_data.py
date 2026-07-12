#!/usr/bin/env python3
import argparse
import math
import os
import re
import shutil
import sys

try:
    import ROOT
except Exception as e:
    print('ERROR: PyROOT is required.')
    print(e)
    sys.exit(e)

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

# Load Delphes
DELPHES_DIR = os.environ.get("DELPHES_DIR")

if DELPHES_DIR is None:
    # Try common locations
    possible_paths = [
        os.path.expanduser("~/Delphes"),
        "/usr/local/Delphes",
        "/opt/Delphes",
    ]

    for path in possible_paths:
        if os.path.isdir(path):
            DELPHES_DIR = path
            break

if DELPHES_DIR is None:
    raise EnvironmentError(
        "Delphes not found. Set DELPHES_DIR environment variable."
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



