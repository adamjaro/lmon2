# Example on using LPS for low-Q2 reconstruction

## Scope

Full chain from physics event generation to reconstructed tracks

## Sequence of steps:

### (1) Electrons from quasi-real photoproduction by GETaLM event generator

- Command: <code> GETaLM getalm_quasireal.ini </code>

- Produces: qr_17p8x275.hepmc3.tree.root

### (2) Geant4 simulation

- Command: <code> run_lmon2 -m lmon2_batch.mac </code>

- Produces: lmon2.root

### (3) Create LPS response

- Command: <code> run_TagTpix4Reco.py --config make_lps.ini </code>

- Produces: lps_response.root

### (4) Reconstruction

- Command: <code> run_TagTpix4Reco.py --config rec_lps.ini </code>

- Produces: tracks.root

