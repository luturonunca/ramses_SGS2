---
name: Ramses_developing
description: Workflow guardrails for RAMSES feature updates (passive scalars, particle families, I/O, MPI packing). Use to agree on steps before any code changes.
---

# Ramses_developing

## Core rule
- Do not edit code until the user explicitly agrees on the approach for the current task.

## Workflow
1) Restate the task and ask clarifying questions before any edits.
2) Identify impacted areas and confirm ordering with the user.
3) Propose a minimal change set and wait for approval.
4) Implement in small steps with checkpoints (commit/push if requested).

## Common impact checklist (use only if relevant)
- Passive scalars: index chain in `hydro/read_hydro_params.f90`, `hydro/hydro_parameters.f90`, `NVAR` in build files.
- Star inheritance: `pm/star_formation.f90`, `pm/pm_commons.f90`, `pm/init_part.f90`.
- Restart/output: `pm/output_part.f90`, `pm/init_part.f90`, `hydro/output_hydro.f90`.
- MPI packing: `pm/particle_tree.f90` data width and packing order.
- Particle families: `pm/pm_commons.f90` constants/functions, filters in `pm/move_fine.f90`, `pm/synchro_fine.f90`, clump/feedback filters.

## Communication style
- Confirm before each code edit phase.
- Keep changes focused and ordered; avoid bundling unrelated edits.
