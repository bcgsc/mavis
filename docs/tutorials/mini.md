# MAVIS (Mini) Tutorial

This tutorial is based on the data included in the tests folder of
MAVIS. The data files are very small and this tutorial is really only
intended for testing a MAVIS install. The data here is simulated and
results are not representative of the typical events you would see
reported from MAVIS. For a more complete tutorial with actual fusion
gene examples, please see the [full tutorial](../../tutorials/full/).

The first step is to clone or download a zip of the MAVIS repository
(<https://github.com/bcgsc/mavis>). You will need the tests directory.
The tag you check out should correspond to the MAVIS version you have
installed

```bash
git clone https://github.com/bcgsc/mavis.git
git checkout <VERSION_TAG>
mv mavis/tests .
mv mavis/Snakefile .
rm -r mavis
```

Now you should have a folder called `tests` in your current directory. Since this is a trivial
example, it can easily be run locally. However in order to run the snakemake file you will need
to have the config validation module `mavis_config` installed which has minimal dependencies.

```bash
pip install mavis_config
```

Now you are ready to run MAVIS. This can be done in a single command using snakemake.

```bash
snakemake -j 1 --configfile=tests/mini-tutorial.config.json -s Snakefile
```

Which will run the mini tutorial version and output files into a folder called `output_dir` in the
current directory
