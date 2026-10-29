MetaHopper Conda environments
=============================

Create the environments:

    mamba env create -f metahopper.yml
    mamba env create -f checkm.yml
    mamba env create -f unicycler.yml

or replace `mamba` with `conda`.

Run MetaHopper from the main environment:

    conda activate metahopper
    python MetaHopper.v3.py ...

MetaHopper searches for CheckM in environments named `checkm` / `checkm_env`
and Unicycler in environments named `unicycler` / `unicycler_env`, so the
environment names in these YAML files match its defaults.

CheckM database
---------------

The CheckM software package does not include its reference database. After
creating the `checkm` environment, configure/download the CheckM data separately,
for example using CheckM's normal `checkm data setRoot ...` workflow, or pass the
database location to MetaHopper with `--checkm-data-path`.

DIAMOND database
----------------

MetaHopper also requires a taxonomy-enabled DIAMOND NR database built with
--taxonmap, --taxonnodes, and --taxonnames. That database is external data and
is not installed by Conda.
