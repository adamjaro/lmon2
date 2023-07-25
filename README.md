# Luminosity monitor and electron tagger for the EIC, version 2

## Dependencies

- Geant4
- ROOT6
- boost

## Optional dependency

- HecMC3

## Steps to checkout the repository, compile and run

### Checkout and build

<pre><code> git clone https://github.com/adamjaro/lmon2.git </pre></code>
<pre><code> cd lmon2 </pre></code>
<pre><code> mkdir build </pre></code>
<pre><code> cd build </pre></code>
<pre><code> cmake ../ </pre></code>
<pre><code> make </pre></code>

### Initialize for terminal

In lmon2 directory, do:

- bash: <code> . init-lmon2.sh </code>
- tcsh: <code> source init-lmon2.csh </code>

### Run in batch

<pre><code> run_lmon2 -m run.mac </pre></code>

### Run with visualization

<pre><code> run_lmon2 --vis run.mac </pre></code>

## All command line options

### Compulsory (one of)

- <code> --mac, -m [name of batch macro] </code> : run in batch mode with the batch macro
- <code> --vis [name of visualization macro] </code> : run with visualization with the visualization macro

### Optional

- <code> --seed, -s [seed value] </code> : seed for Geant4, default: 0
- <code> --phys [name of physics list] </code> : physics list for Geant4, default: FTFP_BERT
- <code> --optics [1 or 0] </code> : turn optics processes on when set to 1, turn off when 0, default: 0

## Source code directories

- <code> base </code> : Basic package functionality
- <code> beamline </code> : Implementation of beam and vacuum elements
- <code> calorimeters </code> : Calorimeter detectors
- <code> reconstruction </code> : Reconstruction implementation
- <code> trackers </code> : Tracker detectors

## Configuration and macro directories

- <code> geometry </code> : Examle geometry configurations
- <code> macro </code> : Run and visualization macros
- <code> plots </code> : Plotting routines

## Source files in lmon2 top directory

- <code> CMakeLists.txt </code> : Top level cmake
- <code> init-lmon2.csh </code> : Initialization for tcsh
- <code> init-lmon2.sh </code> : Initialization for bash shell
- <code> README.md </code> : Top level readme
- <code> run_lmon2.cxx </code> : Source for run_lmon2 executable








