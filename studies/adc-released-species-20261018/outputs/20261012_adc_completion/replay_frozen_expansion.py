#!/usr/bin/env python3
"""Reproduce the immutable experiment in a fresh result directory."""
import argparse
import importlib.util
from pathlib import Path

here = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output-dir", default="scoring_reproduction")
args = parser.parse_args()
destination = (here / args.output_dir).resolve()
if not destination.is_relative_to(here) or destination == here / "scoring":
    raise ValueError("Choose a fresh directory within the experiment; scoring/ is protected.")
spec = importlib.util.spec_from_file_location("frozen_adc_completion", here / "run_frozen_expansion.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.OUT = destination
module.main()
