---
published: false
---

## compare SVs options

|Tool Suite|Pipeline|Action|
|---|---|---|
|HiFiRe3|compare|SVs|

```
hf3 compare SVs --help
```

```
compare: Compare variants called by the analyze pipeline across multiple related samples

SVs: compare structural variant junctions between multiple related samples

library-properties:
  -P,--sequencing-platform  <string> platform used to sequence reads (Illumina_2x150|Aviti_2x150|Aviti_1x300|Ultima|ONT|PacBio) *REQUIRED*
  -L,--library-type    <string> how inserts were prepared for --sequencing-platform (Nextera|TruSeq|Elevate|Ultima|Ligation|Rapid|HiFi) *REQUIRED*
  -I,--min-selected-size    <integer> smallest insert in bp selected BEFORE adapter ligation [0]
  -V,--selected-size-cv     <double> estimated coefficient of variation for --min-selected-size used to calculate min-allowed-size and max-allowed-size [0.15]
  -z,--min-allowed-size     <integer> use this value for min-allowed-size instead of calculating from --min-selected-size and --selected-size-cv [0]

genome:
  -g,--genome          <string> UCSC-compatible name of the reference genome to use (e.g., hg38) *REQUIRED*
  -G,--genome-dir      <string> directory with indexed <--genome>.fa or genome.fa file [.../mdi/resources/genomes/<--genome>] 
  --use-all-chroms     <boolean> use all chromosomes as they are found in genome fasta file without filtering 
  --is-composite-genome     <boolean> <--genome> is a composite of >1 reference (e.g., hs1_dm6 with cross-species spike-in) 

restriction-enzyme:
  -e,--enzyme-name     <string> name of the restriction enzyme used to cleave genomic DNA (from shared/modules/REs/blunt_enzymes.csv if not NA) [NA]

targets:
  -y,--targets-bed     <string> path to a BED file definining genomic regions targeted during sequencing [NA]
  -Y,--region-padding  <integer> bp of adjacency padding applied to target regions in --targets-bed-file [0]

compare:
  -P,--compare-project-dir  <string> directory with multiple sample-level subfolders whose `analyze` pipeline files will be compared [TASK_DIR] 
  -S,--compare-sample-dirs  <string> override --project-dir with a comma-delimited list of `analyze` output directories [--compare-project-dir/*] 

output:
  -O,--output-dir      <string> the directory where output files will be placed; must already exist *REQUIRED*
  -N,--data-name       <string> simple name for the data (e.g., sample) being analyzed (no spaces or periods) *REQUIRED*

push:
  --push-server        <string> external server domain name, e.g, on AWS, to which data packages should be pushed with scp 
  --push-dir           <string> directory on --push-server to which data packages will be pushed [/srv/data]
  --push-user          <string> valid user name on --push-server, authorized by --push-key [ubuntu]
  --push-key           <string> path to an ssh key file for which --push-user has a public key on --push-server [~/.ssh/mdi-push-key.pem]

version:
  -v,--version         <string> the version to use of the tool suite that provides the requested pipeline [latest]

resources:
  -m,--runtime         <string> execution environment: one of direct, conda, container, singularity, auto [auto]
  -p,--n-cpu           <integer> number of CPUs used for parallel processing [1]
  -u,--n-gpu           <character> number [and type] of GPUs used for data processing, as [gpu-type:]n-gpu [0]
  -r,--ram-per-cpu     <string> RAM allocated per CPU (e.g., 500M, 4G) [4G]
  -t,--tmp-dir         <string> directory used for small temporary files (recommend SSD) [/tmp]
  -T,--tmp-dir-large   <string> directory used for large temporary files (generally >10GB) [/tmp]

job-manager:
  --email              <string> email address of the user submitting the job [nobody@nowhere.edu]
  --account            <string> name of the account used to run a job on the server [NA]
  --time-limit         <string> time limit for the running job (e.g., dd-hh:mm:ss for slurm --time) [10:00]
  --partition          <string> slurm --partition (standard, gpu, largemem, viz, standard-oc) [standard]
  --exclusive          <boolean> ensure that only your Slurm job runs on a node; sets --ram-per-cpu to 0 

workflow:
  -f,--force           <boolean> execute certain actions and outcomes without prompting (create, rollback, etc.) 
  -R,--rollback        <integer> revert to this pipeline step number before beginning at the next step (implies --force) 
  -q,--quiet           <boolean> suppress the configuration feedback in the output log stream 

help:
  -h,--help            <boolean> show pipeline help 
  -d,--dry-run         <boolean> only show parsed variable values; do not execute the action 

```
