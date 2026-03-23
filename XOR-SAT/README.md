# XOR-SAT

A small research-oriented pipeline for formulating triangulation distance constraints as a CNF/XOR-SAT instance and solving them with CryptoMiniSat.

The project converts triangulation data from JSON into a custom `.in` format, generates CNF chunks with a C++ OpenMP binary, merges them into a single formula, and runs `cryptominisat5` to determine satisfiability.

This repository is structured around running experiments locally or on the ETH Euler cluster.

## What the pipeline does

Given an input instance containing points and multiple triangulations, the pipeline:

1. converts the JSON input into the internal `.in` format,
2. reads per-triangulation distance bounds from `config/distances.txt`,
3. optionally restricts which triangulations are used via `config/selection.txt`,
4. generates CNF/XOR constraints with the C++ executable,
5. merges all generated `.cnf` fragments into one `total.cnf`,
6. runs `cryptominisat5`,
7. writes the final verdict to `verdict.txt`.

## Repository layout

```text
XOR-SAT/
├── euler/
│   ├── jobs/
│   │   ├── cmsat.sh
│   │   ├── exact.sh
│   │   ├── linear_pairwise.sh
│   │   └── pairwise.sh
│   └── scripts/
│       ├── CMakeLists.txt
│       ├── Utils.hpp
│       ├── all_pairs.hpp
│       ├── convert_format.py
│       ├── exact_run.py
│       ├── main.cpp
│       ├── merge.py
│       ├── pairwise_linear_run.py
│       ├── pairwise_run.py
│       └── run_all.py
└── README.md
```

## Prerequisites

You need:

- Linux or macOS shell environment
- Python 3.10+
- CMake 3.26+
- A C++17 compiler with OpenMP support
- `cryptominisat5`

Tested assumptions from the code:

- the solver binary is expected at `euler/cgshop26/XOR-SAT/cryptominisat5`
- the compiled generator binary is expected at `euler/scripts/build/2sat_linux`
- the scripts are run from within the `euler/` directory or from paths consistent with the relative paths used in the code

## Expected directories

The code uses several directories that are not included in the zip but are expected to exist next to `scripts/`.
Create them before running the pipeline:

```bash
cd euler
mkdir -p data config formulations edges validate quads scripts/build cgshop26/XOR-SAT
```

If `cryptominisat5` is installed globally, you can either:

- copy or symlink it to `euler/cgshop26/XOR-SAT/cryptominisat5`, or
- modify `scripts/run_all.py` to point to your local installation

Example using a symlink:

```bash
ln -s "$(which cryptominisat5)" cgshop26/XOR-SAT/cryptominisat5
```

## Build instructions

Build the C++ binary from `euler/scripts`:

```bash
cd euler/scripts
cmake -S . -B build
cmake --build build -j
```

This should produce:

```text
euler/scripts/build/2sat_linux
```

## Input format

The main input is expected at:

```text
euler/data/<INSTANCE_NAME>.json
```

The JSON must contain:

```json
{
  "points_x": [0, 1, 2],
  "points_y": [0, 1, 0],
  "triangulations": [
    [[0, 1], [1, 2], [0, 2]],
    [[0, 1], [1, 2], [0, 2]]
  ]
}
```

`convert_format.py` transforms this into:

```text
euler/data/<INSTANCE_NAME>.in
```

## Configuration files

### `config/distances.txt`

A space-separated list of distance bounds, one per triangulation.

Example:

```text
3 3 4 2
```

### `config/selection.txt`

This file controls which triangulations are included.

Format:

- first number = number of selected triangulations
- followed by the selected triangulation indices

Examples:

Use all triangulations:

```text
0
```

Use triangulations `1` and `3` only:

```text
2 1 3
```

## How to run

### Standard pipeline

From the `euler/` directory:

```bash
cd euler
python3 scripts/run_all.py <INSTANCE_NAME>
```

Example:

```bash
python3 scripts/run_all.py random_instance_459_320_3
```

This runs:

1. `scripts/convert_format.py`
2. `scripts/build/2sat_linux`
3. `scripts/merge.py`
4. `cryptominisat5`

### Output files

After a successful run, the important outputs are:

- `data/<INSTANCE_NAME>.in` - converted input file
- `formulations/*.cnf` - generated CNF fragments
- `total.cnf` - merged formula
- `scripts/output.txt` - raw CryptoMiniSat output
- `verdict.txt` - final status (`SAT` or `UNSAT`)
- `validate/validate_path.txt` - remembers which instance was run

## Additional scripts

### `exact_run.py`

Runs the full pipeline for multiple distance vectors listed in `../distances.txt` and stops as soon as one is satisfiable.

Run:

```bash
cd euler/jobs
python3 ../scripts/exact_run.py <INSTANCE_NAME>
```

### `pairwise_run.py`

Computes pairwise distances using repeated SAT calls and a binary-search-like process. This script currently contains hardcoded initialization data and is meant for experiment-specific runs.

Run:

```bash
cd euler/jobs
python3 ../scripts/pairwise_run.py <INSTANCE_NAME> <m> <initial_distance>
```

### `pairwise_linear_run.py`

Similar to `pairwise_run.py`, but checks candidate distances in a more linear fashion. Also contains hardcoded experiment logic.

Run:

```bash
cd euler/jobs
python3 ../scripts/pairwise_linear_run.py <INSTANCE_NAME> <m> <initial_distance>
```

## Running on ETH Euler

The `euler/jobs/` directory contains Slurm job scripts:

- `cmsat.sh`
- `exact.sh`
- `pairwise.sh`
- `linear_pairwise.sh`

Example:

```bash
cd euler/jobs
sbatch cmsat.sh
```

Before submitting, make sure the instance names and paths inside the job scripts match your data.

## Quick start

A minimal end-to-end setup looks like this:

```bash
# 1) enter project
cd XOR-SAT/euler

# 2) create missing directories
mkdir -p data config formulations edges validate quads scripts/build cgshop26/XOR-SAT

# 3) make cryptominisat visible at the expected path
ln -s "$(which cryptominisat5)" cgshop26/XOR-SAT/cryptominisat5

# 4) build the generator
cd scripts
cmake -S . -B build
cmake --build build -j
cd ..

# 5) prepare configuration
printf "0\n" > config/selection.txt
printf "3 3\n" > config/distances.txt

# 6) place your input JSON
# expected path: data/example.json

# 7) run
python3 scripts/run_all.py example

# 8) inspect result
cat verdict.txt
```

## Common issues

### `Could not open input file`

Make sure this file exists:

```text
euler/data/<INSTANCE_NAME>.in
```

Normally it is created automatically from `data/<INSTANCE_NAME>.json` by `convert_format.py`.

### `cryptominisat5` not found

The script currently expects the solver at:

```text
euler/cgshop26/XOR-SAT/cryptominisat5
```

Create a symlink there or edit `scripts/run_all.py`.

### Missing directories

The generator writes into directories like `formulations/`, `edges/`, `validate/`, and `quads/`.
If they do not exist, the pipeline may fail.

### Wrong working directory

Many paths in the code are relative. The safest way is to run the main pipeline from:

```bash
cd euler
python3 scripts/run_all.py <INSTANCE_NAME>
```

## Notes

- Some scripts contain hardcoded instance names or hardcoded distance matrices and may need adaptation for new experiments.
- `main.cpp` currently forces `omp_set_num_threads(48)`. Adjust this in the source if your machine has fewer or more cores.

## Suggested future improvements

- make all paths configurable via CLI arguments
- remove hardcoded solver path
- replace hardcoded experiment parameters in the pairwise scripts
- add a sample input instance to the repository
- add dependency installation instructions for Euler specifically
- add validation for malformed `distances.txt` and `selection.txt`

