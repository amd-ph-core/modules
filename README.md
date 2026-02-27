```md
Org: NCEZID
Contact Email: ncezid_shareit@cdc.gov
Status: Active
Keywords: bioinformatics
Version: N/A
Contract#: 47QFCA23F0058
```

# ![amd-ph-core/modules](docs/images/amd_logo.png)

[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

A repository for hosting [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) module files containing AMD Platform tool-specific process definitions and their associated documentation.

## Table of contents

- [Using existing modules](#using-existing-modules)
- [Adding new modules](#adding-new-modules)
- [Help](#help)
- [Citation](#citation)

## Using existing modules

The module files hosted in this repository define a set of processes for software tools such as `varpipe/coverageanalysis`, `picard/buildbamindex`, `irma` etc. This allows you to share and add common functionality across multiple pipelines in a modular fashion.

We have written a helper command in the `nf-core/tools` package that uses the GitHub API to obtain the relevant information for the module files present in the [`modules/`](modules/) directory of this repository. This includes using `git` commit hashes to track changes for reproducibility purposes, and to download and install all of the relevant module files.

1. Install the latest version of [`nf-core/tools`](https://github.com/nf-core/tools#installation) (`>=2.0`)
2. List the available modules:

```console
$ nf-core modules list remote

                                      ,--./,-.
      ___     __   __   __   ___     /,-._.--~\
|\ | |__  __ /  ` /  \ |__) |__         }  {
| \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                      `._,._,'

nf-core/tools version 2.0

INFO     Modules available from nf-core/modules (master):                       pipeline_modules.py:164

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Module Name                    ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ bandage/image                  │
│ bcftools/consensus             │
│ bcftools/filter                │
│ bcftools/isec                  │
..truncated..
```

3. Install the module in your pipeline directory:

```console
$ nf-core modules --git-remote https://github.com/amd-ph-core/modules.git -b install varpipe/coverageanalysis

                                      ,--./,-.
      ___     __   __   __   ___     /,-._.--~\
|\ | |__  __ /  ` /  \ |__) |__         }  {
| \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                      `._,._,'

nf-core/tools version 2.0

INFO     Installing varpipe/coverageanalysis

INFO     Downloaded 3 files to ./modules/nf-core/modules/varpipe/coverageanalysis
```

4. Import the module in your Nextflow script:

```nextflow
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARPIPE_COVERAGEANALYSIS } from './modules/nf-core/modules/varpipe/coverageanalysis/main'
```

5. Remove the module from the pipeline repository if required:

```console
$ nf-core modules remove varpipe/coverageanalysis

                                      ,--./,-.
      ___     __   __   __   ___     /,-._.--~\
|\ | |__  __ /  ` /  \ |__) |__         }  {
| \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                      `._,._,'

nf-core/tools version 2.0

INFO     Removing varpipe/coverageanalysis

INFO     Successfully removed varpipe/coverageanalysis
```

6. Check that a locally installed nf-core module is up-to-date compared to the one hosted in this repo:

```console
$ nf-core modules lint varpipe/coverageanalysis

                                      ,--./,-.
      ___     __   __   __   ___     /,-._.--~\
|\ | |__  __ /  ` /  \ |__) |__         }  {
| \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                      `._,._,'

nf-core/tools version 2.0

INFO     Linting pipeline: .

INFO     Linting module: varpipe/coverageanalysis

╭─────────────────────────────────────────────────────────────────────────────────╮
│ [!] 1 Test Warning                                                              │
╰─────────────────────────────────────────────────────────────────────────────────╯
╭──────────────┬───────────────────────────────┬──────────────────────────────────╮
│ Module name  │ Test message                  │ File path                        │
├──────────────┼───────────────────────────────┼──────────────────────────────────┤
│ varpipe/coverageanalysis       │ Local copy of module outdated │ modules/amd-ph-core/modules/varpipe/coverageanalysis/  │
╰──────────────┴────────────────────────────── ┴──────────────────────────────────╯
╭──────────────────────╮
│ LINT RESULTS SUMMARY │
├──────────────────────┤
│ [✔]  15 Tests Passed │
│ [!]   1 Test Warning │
│ [✗]   0 Test Failed  │
╰──────────────────────╯
```

<!---

### Offline usage

If you want to use an existing module file available in `nf-core/modules`, and you're running on a system that has no internet connection, you'll need to download the repository (e.g. `git clone https://github.com/amd-ph-core/modules.git`) and place it in a location that is visible to the file system on which you are running the pipeline. Then run the pipeline by creating a custom config file called e.g. `custom_module.conf` containing the following information:

```bash
include /path/to/downloaded/modules/directory/
```

Then you can run the pipeline by directly passing the additional config file with the `-c` parameter:

```bash
nextflow run /path/to/pipeline/ -c /path/to/custom_module.conf
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs + modules in one go for you, to make this process easier.

# New test data created for the module- sequenzautils/bam2seqz
The new test data is an output from another module- sequenzautils/bcwiggle- (which uses sarscov2 genome fasta file as an input).
-->

# CI Runners

We are using self-hosted runners for the CI tests.
