# H → ττ Analysis — Uproot / Awkward Implementation

Analysis of the Higgs boson decay

\[
H \rightarrow \tau\tau
\]

in the \(\mu\tau_h\) final state using CMS Open Data NanoAOD ROOT files.

This repository contains the development of a vectorized version of an original PyROOT event-loop analysis, together with the subsequent debugging and validation of the physics selection, trigger matching, event selection, and Monte Carlo normalization.

---

# 1. Project Overview

The analysis studies the final state

\[
H\rightarrow\tau\tau
\rightarrow \mu\nu_\mu\nu_\tau + \tau_h\nu_\tau
\]

where one tau lepton decays leptonically to a muon and the other tau decays hadronically.

The original analysis was implemented using a conventional Python/PyROOT event loop. Although functional, processing the full Monte Carlo sample was very slow.

The main computational objective of this part of the project was therefore to replace the event-by-event Python loop with vectorized operations using:

- Uproot
- Awkward Array
- NumPy
- ROOT / PyROOT

The physics selection was subsequently checked and modified against the cut-based \(H\rightarrow\tau\tau\) analysis described in the CMS reference publication.

**Reference:**

CMS Collaboration, *Measurements of Higgs boson production in the decay channel with a pair of τ leptons*, Eur. Phys. J. C 83, 562 (2023).

DOI:

```text
10.1140/epjc/s10052-023-11452-8
