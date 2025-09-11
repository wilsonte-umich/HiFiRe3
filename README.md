# HiFiRe3 Tools Suite

The **HiFiRe3** tool suite
carries pipelines and apps for analyzing sequencing data generated 
by the HiFiRe3 family of long read library methods, 
which achieve 
**Hi**gh **Fi**delity error-corrected scoring of both 
single-nucleotide variants (SNVs) and structural variants (SVs)
by controlling input fragment sizes and coverage through a
**Re**striction enzyme-mediated **Re**duced **Re**presentation
approach (hence, HiFiRe3, pronounces "high fire 3").

The steps to using HiFiRe3 are to:
- install the codebase
- build the conda enviroment
- download a reference genome
- PENDING

## Single-suite installation (recommended)

HiFiRe3 is implemented in the
[Michigan Data Interface](https://midataint.github.io/) (MDI),
a framework for developing, installing and running 
Stage 1 HPC **pipelines** and Stage 2 interactive web **apps**.
Because HiFiRe3 is packaged as standalone software, 
we recommend a single-suite installation (see the 
[MDI documentation](https://midataint.github.io/) 
for multi-suite installations), which is accomplished by:
- cloning this tool suite repository
- running _install.sh_ to create a suite-specific MDI installation
- calling _alias.pl_ to create a `hf3` alias to the suite's _run_ utility

### Install this tool suite

```bash
git clone https://github.com/wilsontelab/HiFiRe3.git
cd HiFiRe3
./install.sh
```

To start, answer 'y' (yes) to install the Stage 1 Pipelines, then after a
minute, answer 'n' (no) to skip installation of the Stage 2 Apps.

### Create an alias to the suite's _run_ utility (optional)

```bash
# you can use a different alias if you'd like, e.g., replace hf3 with HiFiRe3
perl alias.pl hf3 
```

Answer 'y' to add the alias to your bash profile.

Reload a new shell to activate the alias for use.

If you prefer not to use an alias, 
you can add the installation directory to your PATH variable,
or `cd` into the directory prior to calling `./run`.

### Test the command line interface (CLI)

For help, call the _run_ utility with no arguments, which describes the format for pipeline calls. 

```bash
# use the alias, if you created it as described above
hf3 
hf3 --help                    # tool suite help

# or call the run utility directly without an alias
cd HiFiRe3
./run
./run --help # etc.
```

## Recommended HiFiRe3 workspace organization

Sections below discuss required input files and optional job files,
and, of course, pipelines generate output files.
The following is an optional but time-tested strategy for organizing 
these files in your HiFiRe3 workspace (*** marks optional folders you
create manually as needed).

The organization uses the best practice of separating folders with your 
job files and input and output data files, with each folder type
using a parallel organization for different transgenes.

```sh
HiFiRe3                    # root folder you created above
├── input                  # *** folder for your input data files
│   └── project1           # *** subfolder for a related set of samples
├── jobFiles               # *** folder for your job configuration files
│   └── project1
├── mdi                    # MDI codebase folder created by HiFiRe3 installation
│   ├── config             # folder with configuration files you may need
│   └── resources/genomes  # folder where genome files are placed by default
└── output                 # *** folder for your output data files
    └── project1
```

## Build the required Conda environment

HiFiRe3 pipelines use 3rd-party software built into a conda environment.
Build that environment in your HiFiRe3 installation as follows:

```sh
hf3 genome conda --create
```

You must have `conda` available in your server environment. If you need to run
a command to make conda available, follow the instructions in
`.../mdi/config/stage1-pipelines.yml` (this is pre-configured to work on
the University of Michigan Great Lakes cluster).

**PENDING**: once HiFiRe3 stabilizes, we will release Singularity containers
that can be used instead of building the conda environment yourself.

## Download a reference genome

HiFiRe3 requires a reference genome and support files. The fastest
and best way to obtain these is using the built-in `genome download` pipeline action.
Download only takes a few minutes and does not need to be repeated.


```sh
hf3 genome download --help
hf3 genome download --genome hs1 --output-dir $PWD/_download_hs1 --data-name hs1
```

In addition, you can use the 
[templates/genome.yml](https://github.com/wilsontelab/HiFiRe3/blob/main/templates/genome.yml)
job file template. 

## Input read data

TODO: update this

HiFiRe3 pipelines accept either PacBio CCS/HiFi or ONT read data in either FASTQ 
or unaligned BAM file formats. Read data files must be found in folder `--read-file-dir`.
If needed, set option `--read-file-prefix` to select only specific files in `--read-file-dir`.

Because tGenLR automatically tracks transgene 5mC methylation at CpG dinucleotides and 
6mA adenylation as used in Fiber-seq, it is recommended to use unaligned bam files with 
MM and ML tags as input.

## Execute a Stage 1 pipeline from the command line

### Job files

tGenLR pipelines can be called entirely using the CLI introduced above. However, you 
are strongly encouraged to create YAML-format job configuration files that define the
parameters for your job and coordinate execution steps.

See [the templates folder](https://github.com/wilsontelab/HiFiRe3/tree/main/templates)
for job file templates for all tGenLR pipelines
and actions, and <https://midataint.github.io/mdi/docs/job_config_files.html>
for extended help on using job files. Job file templates can also be generated with 
command `hf3 transgene template`.

The tGenLR CLI and job files can run pipeline actions either 
inline in the calling shell or by submitting jobs to your server job scheduler,
which is recommended for most use cases. Thus, our most common usage pattern is:

```sh
hf3 submit --dry-run myJob.yml # test the job file to see what will happen
hf3 submit myJob.yml           # submit the job to Slurm or your scheduler
hf3 myJob.yml status           # show the state of all submitted jobs
hf3 myJob.yml top              # monitor a running job
hf3 myJob.yml report           # show a job log report
```

### Workflow sequence

HiFiRe3 has XXX main pipeline actions (once you have downloaded a genome assembly), 
listed here in execution order:
- `xxxx xxxx` = xxxx

See **PENDING** for additional details on what happens in these workflows, and be sure
to use commands such as `hf3 <pipeline> <action> --help` 
or `hf3 <pipeline> template` for complete option information.

## Launch the Stage 2 web apps server

To install and launch the HiFiRe3 apps web server, we recommend using the 
[MDI Desktop app](https://midataint.github.io/mdi-desktop-app),
which allows you to control both local and remote MDI web servers.

After following the instructions to run the Desktop on your local machine
or server, load a HiFiRe3 data package file ending in `.mdi.package.zip`,
into the app interface. 
