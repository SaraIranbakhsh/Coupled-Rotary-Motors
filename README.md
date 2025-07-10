
# Coupled Rotary Molecular Motors Simulation

This repository contains the code used to simulate the dynamics of coupled rotary molecular motors with symmetry mismatch under constant and scaling driving forces. The model is based on the Fokker–Planck formulation and numerically solves the resulting partial differential equations to analyze fluxes, power output, efficiency, and disruption regimes.
## Overview

The simulations explore how the symmetry mismatch between two coupled rotary motors—modeled after components of ATP synthase (Fo and F1)—affects system performance, including output power and energy transduction efficiency. Both constant and symmetry-scaled driving schemes are considered.

Key features:
	•	Finite-difference time and space discretization of the Fokker–Planck equation
	•	Steady-state flux and power computation
	•	Visualization of flux fields and power curves
	•	Reproduction of analytical results in simplified barrierless cases

## Requirements
	•	Python 3.x
	•	NumPy
	•	SciPy
	•	Matplotlib


## Acknowledgment

This project builds upon the codebase from jnlucero96/ATP_response, which provided the original finite-difference implementation for a similar model. We adapted and extended their framework to explore coupled rotary systems with symmetry mismatch.

## Reference

If you use this code, please cite:

Iranbakhsh & Sivak, “Effects of Symmetry Mismatch on the Performance of Coupled Rotary Molecular Motors”, [].

