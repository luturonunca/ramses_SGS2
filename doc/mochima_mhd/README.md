# Mochima MHD run

MHD rerun of Mochima (noAGN) with the same ICs and physics as the hydro run.
The SGS turbulence model used by star formation now also runs with `SOLVER=mhd`.

## 1. Get the code
```bash
git clone git@github.com:luturonunca/ramses_SGS2.git   # first time only
cd ramses_SGS2
git fetch origin
git checkout mhd-mochima
git pull
```

## 2. Compile
```bash
cd bin2
make -f MakefileMHD clean     # required: NVAR/SOLVER differ from the hydro build
make -f MakefileMHD
```
This produces `bin2/ramses_MHD3d`.
Build settings: `NVAR=12`, `SOLVER=mhd`, `PATCH=../patch/hydro/MeanBondi`, `-DQUADHILBERT`.

## 3. Run
```bash
mkdir run_mochima_mhd && cd run_mochima_mhd
cp /path/to/ramses_SGS2/doc/mochima_mhd/mochima_MHD.nml .
cp /path/to/ramses_SGS2/bin2/ramses_MHD3d .
mpirun -np <N> ./ramses_MHD3d mochima_MHD.nml
```
- `initfile(...)` points to `../../IC_Mochima/level_0XX`: edit the paths to your IC location.
- The IC folders must contain `ic_pvar_00004` (the zoom refinement scalar, read into variable 12).
- Adjust `walltime_hrs` / `minutes_dump` in `&OUTPUT_PARAMS` to your queue.

## Namelist changes vs `mochima_noAGN.nml`
- `nrestart=0`: hydro outputs cannot be restarted in MHD.
- `B_ave=2.07d-14`: uniform initial Bz seed of 3e-20 G comoving, as in Martin-Alvarez et al. 2021 (MNRAS, "Towards convergence of turbulent dynamo amplification..."). Code units for this box (25 Mpc/h): B[G] = 1.447e-6 · B_ave / a².
- `ivar_refine=12`: passive scalars move up by 3 in MHD (was 9).
- `riemann='hlld'`, `riemann2d='hlld'`: `hllc` does not exist in the MHD solver.

## Quick sanity check
Look at the run log header: "Hydro var indices" should show `imetal=9`, `idelay=10`, `ivirial1=11`.
In the first outputs, variable 11 (turbulent energy) should be non-zero in dense, sheared gas.
