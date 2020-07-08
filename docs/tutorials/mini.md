# MAVIS (Mini) Tutorial

This tutorial is based on the data included in the tests folder of
MAVIS. The data files are very small and this tutorial is really only
intended for testing a MAVIS install. The data here is simulated and
results are not representitive of the typical events you would see
reported from MAVIS. For a more complete tutorial with actual fusion
gene examples, please see the [full tutorial](../../tutorials/full/).

The first step is to clone or download a zip of the MAVIS repository
(<https://github.com/bcgsc/mavis>). You will need the tests directory.
The tag you check out should correspond to the MAVIS version you have
installed

```bash
git clone https://github.com/bcgsc/mavis.git
git checkout v2.0.0
mv mavis/tests .
rm -r mavis
```

Now you should have a folder called `tests` in your current directory.
You will need to specify the scheduler if you want to test one that is
not the default. For example

```bash
export MAVIS_SCHEDULER=LOCAL
```

Since this is a trivial example, it can easily be run locally. By
default MAVIS in local mode will run a maximum of 1 less than the
current cpu count processes. If you are running other things on the same
machine you may find it useful to set this directly.

```bash
export MAVIS_CONCURRENCY_LIMIT=2
```

The above will limit mavis to running 2 processes concurrently.

Now you are ready to run MAVIS itself. This can be done in two commands
(since the config file we are going to use is already built). First set
up the pipeline

```bash
mavis setup tests/data/pipeline_config.cfg -o output_dir
```

Now if you run the schedule step (without the submit flag, schedule acts
as a checker) you should see something like

```bash
mavis schedule -o output_dir/
```

```text
                        MAVIS: 1.8.4
                        hostname: gphost08.bcgsc.ca
[2018-06-01 12:19:31] arguments
                        command = 'schedule'
                        log = None
                        log_level = 'INFO'
                        output = 'output_dir/'
                        resubmit = False
                        submit = False
[2018-06-01 12:19:31] validate
                        MV_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-1 is NOT SUBMITTED
                        MV_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-2 is NOT SUBMITTED
                        MV_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-1 is NOT SUBMITTED
                        MV_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-2 is NOT SUBMITTED
                        MV_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-3 is NOT SUBMITTED
[2018-06-01 12:19:31] annotate
                        MA_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-1 is NOT SUBMITTED
                        MA_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-2 is NOT SUBMITTED
                        MA_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-1 is NOT SUBMITTED
                        MA_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-2 is NOT SUBMITTED
                        MA_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-3 is NOT SUBMITTED
[2018-06-01 12:19:31] pairing
                        MP_batch-s4W2Go4tinn49nkhSuusrE is NOT SUBMITTED
[2018-06-01 12:19:31] summary
                        MS_batch-s4W2Go4tinn49nkhSuusrE is NOT SUBMITTED
                        rewriting: output_dir/build.cfg
```

Adding the submit argument will start the pipeline

```bash
mavis schedule -o output_dir/ --submit
```

After this completes, run schedule without the submit flag again and you
should see something like

```text
                        MAVIS: 1.8.4
                        hostname: gphost08.bcgsc.ca
[2018-06-01 13:15:28] arguments
                        command = 'schedule'
                        log = None
                        log_level = 'INFO'
                        output = 'output_dir/'
                        resubmit = False
                        submit = False
[2018-06-01 13:15:28] validate
                        MV_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-1 (zQJYndSMimaoALwcSSiYwi) is COMPLETED
                        MV_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-2 (BHFVf3BmXVrDUA5X4GGSki) is COMPLETED
                        MV_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-1 (tUpx3iabCrpR9iKu9rJtES) is COMPLETED
                        MV_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-2 (hgmH7nqPXZ49a8yTsxSUWZ) is COMPLETED
                        MV_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-3 (cEoRN582An3eAGALaSKmpJ) is COMPLETED
[2018-06-01 13:15:28] annotate
                        MA_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-1 (tMHiVR8ueNokhBDnghXYo6) is COMPLETED
                        MA_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-2 (AsNpNdvUyhNtKmRZqRSPpR) is COMPLETED
                        MA_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-1 (k7qQiAzxfC2dnZwsGH7BzD) is COMPLETED
                        MA_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-2 (dqAuhhcVKejDvHGBXn22xb) is COMPLETED
                        MA_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-3 (eB69Ghed2xAdp2VRdaCJBf) is COMPLETED
[2018-06-01 13:15:28] pairing
                        MP_batch-s4W2Go4tinn49nkhSuusrE (6LfEgBtBsmGhQpLQp9rXmi) is COMPLETED
[2018-06-01 13:15:28] summary
                        MS_batch-s4W2Go4tinn49nkhSuusrE (HDJhXgKjRmseahcQ7mgNoD) is COMPLETED
                        rewriting: output_dir/build.cfg
                        run time (hh/mm/ss): 0:00:00
                        run time (s): 0
```

If you see the above, then MAVIS has completed correctly!
