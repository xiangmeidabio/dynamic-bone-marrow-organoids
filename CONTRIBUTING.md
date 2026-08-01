# Contributing

This repository is intended to preserve the analysis associated with a scientific manuscript. Changes should improve reproducibility without silently altering biological decisions.

## Before proposing a change

1. Open an issue describing the affected analysis stage and scientific rationale.
2. State whether the change alters a threshold, sample set, cell label, statistical test, or figure.
3. Link the relevant manuscript section or figure panel.

## Code conventions

- Use English comments and descriptive object names.
- Keep steps in execution order.
- Use repository-relative paths through `project_root`.
- Never install packages or change the working directory inside analysis scripts.
- Write generated files only under `data/processed/` or `results/`.
- Do not commit raw data, large single-cell objects, credentials, or personal identifiers.
- Set and report random seeds for stochastic operations.
- Save `sessionInfo()` or equivalent environment provenance for final runs.

## Validation

Run the changed stage from a clean session, inspect warnings, and compare regenerated numerical results and figures with the expected outputs. Document intentional changes in the pull request and update the reproducibility checklist when appropriate.
