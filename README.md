
# Coupled Rotary Molecular Motors Simulation
This repository contains the simulation and analysis code accompanying a peer-reviewed physics publication.
From an applied perspective, it demonstrates large-scale data generation and quantitative analysis to detect regime changes in complex systems using statistical and information-theoretic metrics.
Specifically, the code simulates coupled rotary molecular motors with symmetry mismatch using a Fokker–Planck–based model, numerically solving partial differential equations to analyze fluxes, power output, efficiency, and disruption regimes.

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

Iranbakhsh & Sivak, “Effects of Symmetry Mismatch on the Performance of Coupled Rotary Molecular Motors”, [Phys. Rev. E 112, L052103].

