# HiFiRe3 Tool Suite

The **HiFiRe3** ("high fire") tool suite
carries pipelines and apps for analyzing sequencing data 
from library methods that
achieve **Hi**gh **Fi**delity error-corrected scoring of both 
single nucleotide variants (SNVs) and structural variants (SVs)
by controlling input fragment sizes and/or by
**Re**striction enzyme-mediated **Re**duced **Re**presentation.

**IMPORTANT**: HiFiRe3 tools assume that any
insert size selection during library preparation was performed 
_before adapter ligation_, targeting molecules from 1N to 2N bp. 
They also assume that input libararies are either PCR-free
or used pooled PCR with unique dual indices to minimize chimeric 
template switching.

The steps to using HiFiRe3 are to:
- install the codebase
- build the conda enviroment
- download and digest a reference genome using `prepare genome`
- if needed, basecall reads using `basecall PacBio` or `basecall ONT`
- align and analyze basecalled reads using:
    - `analyze fragments`
    - `analyze SVs` and/or `analyze SNVs`
- compare samples using `compare SVs` and/or `compare SNVs`
- visualize results in the R Shiny apps

## Single-suite installation (recommended)

HiFiRe3 is implemented in the
[Michigan Data Interface](https://midataint.github.io/) (MDI)
for developing, installing and running 
Stage 1 HPC **pipelines** and Stage 2 interactive web **apps**.
Because HiFiRe3 is packaged as standalone software, 
we recommend a single-suite installation (see the 
[MDI documentation](https://midataint.github.io/) 
for multi-suite installation), which is accomplished by:
- cloning this tool suite repository
- running _install.sh_ to create a suite-specific MDI installation
- calling _alias.pl_ to create a `hf3` alias to the suite's command line interface (CLI)

### Install this tool suite

```bash
git clone https://github.com/wilsontelab/HiFiRe3.git
cd HiFiRe3
./install.sh
```

To start, answer 'y' (yes) to install the Stage 1 Pipelines, then after a
minute, answer 'n' (no) to skip installation of the Stage 2 Apps (for now).

### Create an alias to the suite's command line interface (optional)

```bash
# you can use a different alias if you'd like, e.g., replace hf3 with HiFiRe3
perl alias.pl hf3 
```

Answer 'y' to add the alias to your bash profile, then 
reload a new shell to activate the alias for use.

If you prefer not to use an alias, 
add the installation directory to your PATH variable
or `cd` into the directory prior to calling `./hf3`.

### Test the command line interface (CLI)

For help, call the CLI with no arguments, which describes the format for pipeline calls. 

```bash
# use the alias, if you created it as described above
hf3 
hf3 --help

# or call the CLI directly without an alias
cd HiFiRe3
./hf3
./hf3 --help
```

## Recommended HiFiRe3 workspace organization

Sections below discuss required input files and optional job files,
and, of course, pipelines generate output files.
The following is an optional but time-tested strategy for organizing 
these files in your HiFiRe3 workspace (*** marks optional folders you
create manually as needed). The parallel organization uses the best 
practice of separating folders with your job files and input and output 
data files.

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

HiFiRe3 pipelines use version-controlled 3rd-party software built into a 
[conda](https://docs.conda.io/)
environment. Build that environment in your HiFiRe3 installation as follows:

```sh
hf3 analyze conda --create
```

You must have `conda` available in your server environment. If you need to run
a command to make conda available, follow the instructions in
`.../mdi/config/stage1-pipelines.yml` (HiFiRe3 is pre-configured to work on
the University of Michigan Great Lakes cluster).

In a shared server environment, the conda build command may get killed by the host.
If that happens, run the command on a cluster worker node with sufficient resources,
e.g., 4 CPU with 4G RAM per CPU.

**PENDING**: once HiFiRe3 stabilizes, we will release Singularity containers
that can be used instead of building the conda environment yourself.

## Execute a pipeline from the command line

### Job files

HiFiRe3 pipelines can be called entirely using the CLI introduced above. However, you 
are encouraged to create YAML-format job configuration files that define the
parameters for your job and coordinate execution steps.

See [the templates folder](https://github.com/wilsontelab/HiFiRe3/tree/main/templates)
for job file templates for all HiFiRe3 pipelines
and actions, and <https://midataint.github.io/mdi/docs/job_config_files.html>
for extended help on using job files. Job file templates can also be generated with 
command `hf3 <pipeline> template`, e.g., `hf3 basecall template`.

The HiFiRe3 CLI and job files can run pipeline actions either 
inline in the calling shell or by submitting jobs to your server job scheduler;
the latter is recommended for most use cases. Thus, our most common usage pattern is:

```sh
hf3 inspect myJob.yml          # check the formatting of your job file
hf3 mkdir myJob.yml            # create any missing output directories
hf3 submit --dry-run myJob.yml # test the job file to see what will happen
hf3 submit myJob.yml           # submit the job to Slurm or your scheduler
hf3 myJob.yml status           # show the state of all submitted jobs
hf3 myJob.yml top              # monitor a running job
hf3 myJob.yml report           # show a job log report
hf3 myJob.yml ls               # show the contents of a job's output diretory
```

### Workflow sequence

HiFiRe3 has _pipelines_ each with associated _actions_ with descriptive names, 
listed here in execution order:
- `prepare genome` (where `prepare` is the _pipeline_ and `genome` is the _action_)
- `basecall ONT` or `basecall PacBio`
- `analyze fragments`
- `analyze SVs`
- `analyze SNVs`
- `compare SVs`
- `compare SNVs`

Required/common options are described below; use 
`hf3 <pipeline> <action> --help` or `hf3 <pipeline> template` 
for complete option information, or see the output of all action help commands
[here](https://github.com/wilsontelab/HiFiRe3/tree/main/options).

### Universally required options

Options `--output-dir/-O` and `--data-name/-N` are required by all pipeline actions.
Sample-level output files are placed into directory `<--output-dir>/<--data-name>`.

## Download and digest a reference genome

HiFiRe3 requires a reference genome assembly and supporting annotation files. Obtain and create 
them using the `prepare genome` pipeline action.

```sh
hf3 prepare genome --help
hf3 prepare genome --genome hs1 --output-dir $PWD/_prepare_hs1 --data-name hs1
```

Even better, use the 
[templates/prepare_genome.yml](https://github.com/wilsontelab/HiFiRe3/blob/main/templates/prepare_genome.yml)
job file template. 

The action downloads all required genome assembly and annotation files and
tabulates blunt restriction enzyme (RE) sites in the assembly sequences.
These steps only need to be performed once per genome. 
Many blunt REs suitable for HiFiRe3 are tabulated by _in silico_ digestion in a single pass.
A given library preparation will only use one of these enzymes (or none at all).

### Optionally create a composite genome for mixed species analysis

Especially during development, quality assessment of SV artifact rates is
enhanced by mixing samples from two species prior to libary preparation and sequencing.
The `combine genomes` pipeline action will create the required composite reference,
which should then be provided as option `--genome` in later analysis steps while
also setting option `--is-composite-genome`.

```sh
hf3 combine genomes --help
hf3 combine genomes --genome1 hs1 --genome2 dm6 --output-dir $PWD/_prepare_hs1_dm6 --data-name hs1_dm6
```

Combining genomes requires that each individual genome was previously prepared as above.
Even better, use the 
[templates/combine_genomes.yml](https://github.com/wilsontelab/HiFiRe3/blob/main/templates/combine_genomes.yml)
job file template to coordinate all work.

The resulting composite reference assembly carries the original chromosome names
suffixed with the respective source genome, e.g., `chr1_hs1` and `chr4_dm6`.

## Input read data

HiFiRe3 accepts read input files from the following sequencing platforms,
read configurations, and library types, communicated via options 
`--sequencing-platform` and `--library-type`:

| --sequencing-platform | --library-type |
|---|---|
| Illumina_2x150 | Nextera |
| Illumina_2x150 | TruSeq |
| Aviti_2x150 | Nextera |
| Aviti_2x150 | TruSeq |
| Aviti_2x150 | Elevate |
| Aviti_1x300 | Nextera |
| Aviti_1x300 | TruSeq |
| Aviti_1x300 | Elevate |
| Ultima | Ultima |
| ONT | Rapid |
| ONT | Ligation |
| PacBio | HiFi |

All tools assume you performed a valid HiFiRe3 library preparation,
with any insert size selection performed _before adapter ligation_. However,
the pipeline is flexible with regard to which error correction methods you exploit, 
e.g., you can forego size selection or RE-mediated DNA fragmentation as suits your needs
(see option details below).

## Basecalling long reads

The `basecall` pipeline finishes processing PacBio or ONT
reads into the unaligned BAM format required for later pipeline actions. 
However, what happens during PacBio and ONT basecalling is different.

```sh
hf3 basecall --help
hf3 basecall PacBio --help
hf3 basecall ONT --help
hf3 basecall ...
```

As throughout, using the 
[templates/basecall.yml](https://github.com/wilsontelab/HiFiRe3/blob/main/templates/genome.yml)
job file template can simplify use of the `basecall` pipeline.

For short-read libraries, input reads will have been basecalled by the sequencing platform
and you should proceed directly to fragment analysis.

### PacBio, from unaligned single-strand consensus bam files

For PacBio HiFiRe3 with sufficiently small inserts, initial basecalling is 
performed during your Revio or Vega run, which must be configured with the following settings:
- Library Type = Standard (e.g., WGS)
- Insert Size = the middle of your selected size range, i.e., (N + 2N) / 2 = 1.5N
- Movie Acquisition Time = at least 24 hours (or the maximum allowed)
- Data Options:
    - Include Base Kinetics = YES (or NO, if you do not intend to look for damaged bases)
    - Consensus Mode = Strand **<<<< IMPORTANT!**

These settings ensure that you will have:
- sufficient time to get around shorter HiFiRe3 reads _more times than a normal HiFi library_ (ideally at least 4 subreads on each strand)
- consensus output files generated for each strand independently

Point `basecall PacBio` at your strand consensus unaligned BAM files using
option `--pacbio-dir`. HiFiRe3 compares the strand consensuses to each other and assigns a final 
base call and quality score. When heteroduplex DNA is detected, it might result from:
- a base mismatch characterized by two high-quality but non-Watson-Crick paired bases
- a damaged base on one of the two strands opposite an undamaged base

HiFiRe3 uses base quality and kinetic information to determine which case is more likely.
For mismatched high-quality bases, an N base and low quality are reported, with the values
of the mismatched bases reported in the XX tag. 
For base pairs with one high-quality and one low confidence/kinetically unusual base, 
the high-quality base is reported at lower confidence than fully validated base pairs
and information about the low-quality bases is reported in the XX tag.

PacBio reads processed in this way can be used for high-accuracy, error-corrected calling of 
both SVs and SNVs. Alternatively, larger PacBio reads can be subjected to normal basecalling
and used only for error-corrected SV calling.

### ONT, from POD5 files

Unlike PacBio, HiFiRe3 must perform basecalling of ONT reads from primary trace data 
to make best use of the fixed bases at RE-cleaved DNA ends. ONT reads cannot support
high-fidelity SNV calling, but fully support error-corrected calling of rare mosaic 
SVs (often using longer reads than for combined SNV and SV calling with PacBio).

The required input for ONT is one or more POD5 format files obtained from your
nanopore sequencing run. Use option `--pod5-dir` to point at your files in your call to
`basecall ONT`. HiFiRe3 uses 
[Dorado](https://github.com/nanoporetech/dorado)
for high accuracy basecalling and custom scripts
for RE-aware adapter trimming that preserves ONT channel information. See
the options help for information on how to tune ONT basecalling and/or to 
invoke modified basecalling.

### HiFiRe3 tags in unaligned bam files

The final output of either HiFiRe3 `basecall` action is one or more unaligned
BAM files in folder `<--output-dir>/<--data-name>/ubam`. Important custom tags specific to HiFiRe3
or derived from vendor files are (all propagate into aligned BAM files, below):
- ch:i: = ONT channel (absent for other platforms)
- tl:Z: = ONT adapter trim lengths in format `<5' trim>,<3' trim>` (absent for other platforms)

## Shared fragment analysis upstream of SV and SNV calling

The `analyze` pipeline aligns and processes reads with the goal of finding
rare mosaic variants, especially singleton variants, with extremely low baseline artifact rates
for SVs and, mainly for PacBio, SNVs/indels.

The `analyze fragments` action filters and indexes reads to prepare for variant analysis.
Briefly, various quality filters ensure that only high-quality conforming reads and bases are used for
variant calling downstream. Indexing allows for efficient comparison of individual reads
to RE fragment locations and read recovery during visualization.

```sh
hf3 analyze --help
hf3 analyze fragments --help
hf3 analyze fragments ...
```

### Read alignment using minimap2

The first step in `analyze fragments` aligns properly basecalled reads to the reference genome using 
[minimap2](https://github.com/lh3/minimap2), 
for both short and long reads. Use option --`read-file-dir` and `--read-file-prefix` to communicate 
the path to your input read files if they were not generated with `basecall PacBio` or `basecall ONT`.

The following custom BAM tags are added by HiFiRe3
alignment, in addition to the basecalling tags added above:
- fm:Z: = metadata describing merging of paired-end reads as mergeLevel:nRead1:nRead2 (absent for single-read platforms)
- hv:i: = bit-encoded flag of alignment and read variant status

### Inferring RE site locations in reduced representation reads

When option `--enzyme-name` is set to something other than NA, HiFiRe3 expects reads
to arise from a subset of all possible blunt RE fragments in a sample's genome, with each
fragment sequenced sufficient times to establish its consensus genotype without external information 
This is known as a _reduced representation_ libary. While most fragments can be predicted from the
_in silico_ RE site analysis performed by `prepare genome`, others will be genotype-specific due to RFLPs. 
Therefore, when applicable, `analyze fragments` next uses the ends of sequenced inserts to locate recurring sample-specific RE sites. 

Clonal SNV-containing fragments will be among the set of located fragments if they are within the 
(selected) insert size range. Only these reads are used for SNV discovery, where reference and 
SNV/indel-containing reads will be ~the same size.

In contrast, rare mosaic SV-containing fragments may not end at located RE sites depending on the 
relative sizes of clonal and SV fragments, i.e.,  an SV fragment may be in a library's insert 
size range even though  the parental fragment(s) it bridges are not.
For this reason, both _in silico_ and located RE sites are used in subsequent steps when
matching read outer ends to RE sites. Any fragment ending at expected RE sites of either type
are used for SV discovery, which, unlike clonal fragments, can be located ~anywhere in the genome 
or target regions.

These steps do not apply and are skipped for randomly sheared library inputs.

## Mosaic variant calling

Depending on your application and sequencing platform, 
you may wish to call SVs, SNVs/indels, or both from your processed reads, which
is done in distinct actions of the `analyze` pipeline:

```sh
hf3 analyze SVs --help
hf3 analyze SVs ...
hf3 analyze SNVs --help
hf3 analyze SNVs ...
```

The error correction methods available during SV and SNV calling depend on the nature
of the input libraries. Error-corrected SNV calling is mainly used
for sufficiently small PacBio reads processed using `basecall PacBio` as described above.
SV error correction enforces filters against end-to-end chimeras that may variably include:
- rejecting SVs with adapters inserted at junctions (always applied)
- rejecting junctions flanked by low quality bases or alignments (always applied)
- rejecting junctions that match a known RE fragment endpoint (if a blunt RE was used to fragment the genomic DNA)
- requiring junctions to be within 1N bp of insert ends (if size selection was performed before adapter ligation)

Use options `--sequencing-platform` and `--library-type` to tell the pipeline how
to perform adapter and alignment quality checks based on the nature of the reads.

Use option `--enzyme-name` to communicate the blunt RE used to prepare a HiFiRe3 library
upstream of adapter ligation or tagmentation. Leave `--enzyme-name`as the 
default value 'NA' if your library is a sheared DNA library with no RE cleavage.

Use options `--min-selected-size`, `--selected-size-cv`, and `--min-allowed-size` 
to communicate the selected DNA insert size range. Leave both `--min-selected-size` and 
`--min-allowed-size` at the default value of 0 if size selection was not performed. Otherwise,
you will most commonly set `--min-selected-size` to communicate the lower insert size cutoff,
which should be strictly established, e.g., using a BluePippin device. If needed, the values
calculated from `--min-selected-size` and `--selected-size-cv` can be overridden by setting
`--min-allowed-size` to a non-zero value. Whether you set `--min-selected-size` or 
`--min-allowed-size`, these options establish the 1N insert size value used to enforce
size-based SV error correction.

RE-based and size-based SV error correction methods are substantially redundant but 
independent, so both can be used to achieve maximal SV error correction. There is little 
point in running `analyze SVs` pipelines if neither method was used during library preparation. 

Importantly, each SV error correction method only allows the pipeline to reliably reject 
end-to-end chimeras. Many middle-to-middle chimeras arising from DNA fragmentation after 
size-selection cannot be error corrected, so care must be taken to suppress their formation 
by aggressively limiting DNA fragmentation that occurs between size selection and adapter ligation.

The final output of both variant calling actions is a data package suitable for
loading into the interactive visualization app.

## Comparing variants across samples

HiFiRe3 applications often seek to distinguish single-molecule variants,
i.e., SVs and SNVs detected only once in a sample, from those detected
multiple times. Multiply-sequenced variants must have arisen from different 
input DNA molecules in the PCR-free libraries HiFiRe3 is mainly intended to 
process, whereas single-molecule variants are consistent with mutations
that occurred during an experiment or recently in a tissue or cell lineage.

Sensitivity for detecting multiply-sequenced variants increases if data from 
many samples from the same DNA source can be combined, which is the
job of `compare SVs` and `compare SNVs`. These actions are not needed if you 
are only analyzing one sample.

The process for comparing variants across samples is essentially the same as
that used to count the number of reads that called a variant in a single
sample, now tracking sample sources whan applying fuzzy-logic variant comparisons.

## Targeted genomic analysis

Many error-corrected sequencing approaches benefit from focusing reads to specific
genomic regions. The details differ by platform and library type, but commonly
include hybridization capture and, most useful for HiFiRe3, ONT adaptive sampling.
If applicable, use options `--targets-bed` and `--region-padding` to communicate 
the genomic regions to which your reads were targeted. 

During SV calling, HiFiRe3 pipelines only consider reads where the 5' end 
aligned to a target region. This approach is critical for ONT adaptive sampling
where off-target 5' ends will be present in the library in a truncated form that
should not be used for equivalent SV calling to on-target reads.

Note that region targeting is compounded on top of RE-mediated reduced representation
sampling to allow considerable focusing of reads onto specific genomic fragments.

## Launch the interactive apps server

Once all pipelines have finished running, you will want to view and interact
with your data.

To install and launch the HiFiRe3 apps server, we recommend using the 
[MDI Desktop app](https://midataint.github.io/mdi-desktop-app),
which allows you to control both local and remote MDI web servers.

After following the instructions to run the Desktop on your local machine
or server, load a HiFiRe3 data package file ending in `.mdi.package.zip`,
into the app interface. 
