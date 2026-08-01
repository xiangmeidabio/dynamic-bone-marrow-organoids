# Software environment

The original notebook metadata records:

- R 4.2.3 for all R notebooks; and
- Python 3.10.18 for the CellOracle notebook.

The source files did not contain a complete package-version report. Inventing exact versions would create false reproducibility, so this repository separates environment reconstruction from version locking.

## R environment

From the repository root:

```r
install.packages("renv")
renv::init(bare = TRUE)
source("environment/install_r_packages.R")
```

After the analysis runs successfully and the figures match the manuscript, record the validated versions:

```r
renv::snapshot()
writeLines(capture.output(sessionInfo()), "environment/r_session_info.txt")
```

Commit `renv.lock` and `environment/r_session_info.txt` with the paper release. Do not commit `renv/library/`.

## Python environment

Create an isolated Python 3.10 environment and install the listed packages:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r environment/requirements.txt
python -m ipykernel install --user --name bmo-celloracle --display-name "BMO CellOracle"
```

On Windows PowerShell, activate the environment with `.venv\Scripts\Activate.ps1`.

After CellOracle has been rerun and validated, pin exact installed versions in `requirements.txt` or generate a lock file with the environment manager used for the final analysis.

## pySCENIC

pySCENIC is called from `analysis/08_gene_regulatory_networks/01_run_pyscenic.sh`. Record the exact pySCENIC version, ranking database release, motif-annotation release, and transcription-factor list in `docs/reproducibility_checklist.md` before publication.
