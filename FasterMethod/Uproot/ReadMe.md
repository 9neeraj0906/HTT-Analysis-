# H → ττ Analysis — μ + τh Channel

## Overview

This repository contains the development of a cut-based analysis of the Higgs boson decay

\[
H \rightarrow \tau\tau
\]

in the **muon + hadronic tau (\(\mu+\tau_h\)) final state**.

The analysis is being developed using CMS Open Data in NanoAOD ROOT format. The event processing is performed with **Python, Uproot, and Awkward Array**, with ROOT used for histogram storage and visualization.

The main goal is to implement the event selection, trigger requirements, muon–tau candidate construction, Monte Carlo normalization, and mass-related observables needed for an \(H\rightarrow\tau\tau\) analysis.

---

## Analysis Channel

The channel studied in this repository is

\[
H\rightarrow\tau\tau
\rightarrow \mu\nu_\mu\nu_\tau + \tau_h\nu_\tau
\]

where:

- one tau decays leptonically to a muon;
- the other tau decays hadronically;
- the reconstructed final state contains a muon and a hadronic tau candidate.

The analysis therefore focuses on the

\[
\boxed{\mu+\tau_h}
\]

final state.

---

## Software

The analysis uses:

- Python 2.7
- Uproot
- Awkward Array
- NumPy
- ROOT
- Python ROOT bindings

The current analysis code uses the older Awkward Array API because it was developed in the corresponding analysis environment.

---

## Analysis Workflow

The analysis follows the approximate workflow:

```text
CMS NanoAOD ROOT files
        │
        ▼
  Read event branches
        │
        ▼
    MET filters
        │
        ▼
     b-jet veto
        │
        ▼
    lepton vetoes
        │
        ▼
   Muon selection
        │
        ▼
   Trigger selection
        │
        ▼
    Tau selection
        │
        ▼
 Muon–tau trigger matching
        │
        ▼
   μ–τ combinations
        │
        ▼
 Pair selection
        │
        ▼
   Event weights
        │
        ▼
 ROOT histograms
        │
        ▼
 Data / MC comparison
```

The main event-selection and histogram-production function is contained in:

```text
analysis/analysisfxn.py
```

---

## Repository Structure

```text
Uproot/
│
├── README.md
│
├── analysis/
│   ├── analysisfxn.py
│   ├── run_data_part1.py
│   ├── run_data_part2.py
│   ├── run_background.py
│   ├── run_signal.py
│   └── make_plots.py
│
├── plots/
│   ├── visible_mass/
│   ├── transverse_mass/
│   └── reconstructed_mass/
│
├── results/
│   ├── data/
│   ├── background/
│   └── signal/
│
├── docs/
│   ├── analysis_notes.md
│   └── debugging.md
│
├── requirements.txt
│
└── .gitignore
```

Large ROOT input datasets are **not stored in this repository**.

---

# Analysis Code

## `analysisfxn.py`

This is the central physics-analysis file.

It contains the event-level selection and histogram production, including:

- reading NanoAOD branches;
- MET filter selection;
- b-jet veto;
- electron veto;
- muon selection;
- trigger requirements;
- trigger-object matching;
- tau selection;
- muon–tau pairing;
- transverse mass calculation;
- Monte Carlo event weighting;
- histogram filling.

The function used to process an individual ROOT file is:

```python
Analysis_up(...)
```

---

## Data Processing

The data drivers process the available data samples separately:

```text
run_data_part1.py
run_data_part2.py
```

The two outputs can subsequently be combined into a single data histogram.

Data events are assigned unit event weight:

\[
w_{\mathrm{data}}=1.
\]

---

## Background MC

Background Monte Carlo is processed using:

```text
run_background.py
```

The current background sample used in the analysis is a Drell–Yan sample.

For Monte Carlo events, the event weight is based on the generator weight and the sample cross section:

\[
w =
\frac{
w_{\mathrm{gen}}\,
\sigma\,
\mathcal{L}
}{
\sum w_{\mathrm{gen}}
}.
\]

The total generator-weight sum is calculated over the relevant MC files before processing the individual files.

---

## Signal MC

Signal Monte Carlo is processed using:

```text
run_signal.py
```

The signal files are processed using the same central analysis function so that the event-selection procedure remains consistent between data, background, and signal.

---

# Muon Selection

The analysis applies a strict reconstructed-muon selection based on:

- transverse momentum;
- pseudorapidity;
- transverse and longitudinal impact parameters;
- medium muon identification;
- relative isolation;
- primary-vertex requirements.

Trigger-object matching is subsequently applied to associate the reconstructed muon with the appropriate trigger object.

---

# Tau Selection

The hadronic tau candidate selection includes requirements on:

- tau pseudorapidity;
- longitudinal impact parameter;
- DeepTau anti-electron identification;
- DeepTau anti-muon identification;
- DeepTau jet discrimination;
- allowed tau decay modes.

The selected tau information includes quantities such as:

```text
pT
eta
phi
charge
DeepTau VSjet discriminator
```

The tau candidate is subsequently matched to the corresponding tau trigger object.

---

# Muon–Tau Pairing

For each event, all possible reconstructed muon–tau combinations are constructed.

Conceptually:

```text
muons × taus
```

produces all possible pairs in the event.

Pair-level requirements can then be applied to quantities such as:

\[
\Delta R(\mu,\tau),
\]

the relative charge,

\[
Q_\mu Q_\tau,
\]

and the tau transverse momentum.

The selected pair is used to construct the analysis observables.

---

# Visible Mass

One of the current observables is the visible invariant mass of the reconstructed muon and hadronic tau:

\[
m_{\mathrm{vis}} =
m(\mu,\tau_h).
\]

Because neutrinos are not directly reconstructed, this is **not** the full \(\tau\tau\) invariant mass.

A current data/MC comparison of the visible-mass distribution is stored in the `plots/` directory.

---

# Transverse Mass

The transverse mass between the muon and missing transverse momentum is also calculated:

\[
m_T =
\sqrt{
2p_T^\mu E_T^{\mathrm{miss}}
\left[
1-\cos\left(\Delta\phi(\mu,E_T^{\mathrm{miss}})\right)
\right]
}.
\]

This observable is used to study backgrounds containing genuine missing transverse momentum, particularly \(W\)-related backgrounds.

---

# Current Results

The current analysis produces data and Monte Carlo distributions for observables including the visible mass.

Example:

```text
plots/
└── visible_mass/
    └── combined_visible_mass.pdf
```

The current visible-mass comparison is a useful intermediate result, but it should not be interpreted as a complete reproduction of the final CMS \(H\rightarrow\tau\tau\) background model.

---

# Current Status

### Implemented

- [x] NanoAOD ROOT event reading
- [x] Awkward Array event processing
- [x] MET filter selection
- [x] b-jet veto
- [x] electron veto
- [x] muon selection
- [x] trigger-object matching
- [x] tau selection
- [x] muon–tau candidate construction
- [x] event weighting
- [x] ROOT histogram production
- [x] data/MC comparison plots
- [ ] Complete background model
- [ ] Final \(H\rightarrow\tau\tau\) mass reconstruction
- [ ] Final statistical interpretation

---

# Known Limitations

This repository represents a developing analysis rather than a full reproduction of the CMS \(H\rightarrow\tau\tau\) analysis.

In particular:

1. The current background treatment does not yet reproduce the complete data-driven background estimation used in a full CMS analysis.

2. The visible mass

\[
m_{\mathrm{vis}}(\mu,\tau_h)
\]

is only an intermediate observable and is not equivalent to the reconstructed Higgs-to-\(\tau\tau\) mass.

3. Trigger and object selections are being validated against the relevant dataset and analysis requirements.

4. The analysis is currently being developed as a cut-based selection rather than a complete statistical measurement.

These limitations are kept explicit so that intermediate results are not presented as a final measurement.

---

# Running the Analysis

The analysis drivers are located in:

```text
analysis/
```

The general workflow is:

```bash
python2 run_data_part1.py
python2 run_data_part2.py
python2 run_background.py
python2 run_signal.py
python2 make_plots.py
```

The exact input ROOT file locations are defined in the corresponding driver scripts and are not included in the GitHub repository.

---

# References

The analysis is based on the methodology and event-selection requirements of the corresponding CMS \(H\rightarrow\tau\tau\) analysis and the available CMS Open Data samples.

Relevant papers and dataset documentation should be added here as the analysis is finalized.

---

## Development Note

This repository contains the actual development of the analysis, including intermediate selections, validation studies, debugging, and ongoing improvements.

The intention is to keep the analysis code transparent and reproducible rather than presenting only the final plots.
