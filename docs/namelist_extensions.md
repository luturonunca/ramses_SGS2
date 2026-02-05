Local Namelist Notes (bin2 build)
=================================

Scope
-----
- This note documents the namelist blocks used by the bin2 build and the
  example run in ../namelist_0.nml.
- It is an add-on to the existing LaTeX user guide (doc/src/*.tex). It does
  not replace or modify those files.

Build context (bin2/Makefile)
-----------------------------
- SOLVER=hydro, GRACKLE=0, NDIM=3, NVAR=10.
- Physics modules linked: amr, hydro, poisson, pm (particles), clump finder,
  sink/AGN feedback.

Where each namelist block is read
---------------------------------
- &RUN_PARAMS, &OUTPUT_PARAMS, &AMR_PARAMS, &POISSON_PARAMS,
  &LIGHTCONE_PARAMS, &MOVIE_PARAMS
  - Read in amr/read_params.f90.
- &INIT_PARAMS, &HYDRO_PARAMS, &REFINE_PARAMS, &BOUNDARY_PARAMS,
  &COOLING_PARAMS, &SF_PARAMS, &FEEDBACK_PARAMS, &UNITS_PARAMS
  - Read in hydro/read_hydro_params.f90 (called during init_hydro).
  - Note: &PHYSICS_PARAMS is explicitly rejected there as deprecated.
- &SINK_PARAMS
  - Read in pm/sink_particle.f90 (read_sink_params).
- &CLUMPFIND_PARAMS
  - Read in pm/clump_finder.f90 (read_clumpfind_params).

Blocks used by ../namelist_0.nml (and their code paths)
------------------------------------------------------
- &RUN_PARAMS
  - Enables cosmo, hydro, poisson, pic, sink, clumpfind.
  - Controls coarse stepping and output frequency used in amr/adaptive_loop.f90
    and amr/amr_step.f90.
- &INIT_PARAMS
  - filetype='grafic' and initfile(*) select IC reading in hydro/init_hydro.f90
    and pm/init_part.f90.
- &AMR_PARAMS / &REFINE_PARAMS
  - Define the grid hierarchy and refinement triggers used in
    amr/refine_utils.f90 and amr/flag_utils.f90.
- &POISSON_PARAMS
  - Sets gravity solver tolerances used in poisson/* routines.
- &HYDRO_PARAMS
  - Selects scheme/riemann/limiter for hydro/godunov_fine.f90.
- &COOLING_PARAMS
  - Cooling + UV background used in hydro/cooling_module.f90.
- &SF_PARAMS
  - Star formation knobs used in pm/star_formation.f90.
- &FEEDBACK_PARAMS
  - Supernova/mechanical feedback settings used in pm/feedback.f90 and
    pm/mechanical_fine.f90.
- &SINK_PARAMS
  - SMBH/AGN sink settings used in pm/sink_particle.f90.
- &CLUMPFIND_PARAMS
  - Clump finder thresholds used in pm/clump_finder.f90.
- &MOVIE_PARAMS
  - Output movie frames in amr/movie.f90.

Documentation gaps vs. current run
----------------------------------
- The LaTeX guide (doc/src/30-runtimeparams.tex) still documents
  &PHYSICS_PARAMS, but this build rejects it in hydro/read_hydro_params.f90.
- Newer blocks used here (&COOLING_PARAMS, &SF_PARAMS, &FEEDBACK_PARAMS,
  &SINK_PARAMS, &CLUMPFIND_PARAMS) are not described in doc/src/*.tex.
- The guide assumes compiling from bin/, while this run uses bin2/Makefile.

Suggested next steps (optional)
-------------------------------
- Add a short section in doc/src/30-runtimeparams.tex that documents the
  new block names and points to the relevant modules.
- Add a one-page build note describing bin2/Makefile as the active build.
