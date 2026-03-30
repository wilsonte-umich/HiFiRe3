---
published: false
---

## basecall ONT options

|Tool Suite|Pipeline|Action|
|---|---|---|
|HiFiRe3|basecall|ONT|

```
hf3 basecall ONT --help
```

```
basecall: Perform basecalling on input read data from the sequencing platform

ONT: use Dorado to basecall reads in ONT POD5 file(s), with RE-specific adapter trimming

library-properties:
  -P,--sequencing-platform  <string> platform used to sequence reads (Illumina_2x150|Aviti_2x150|Aviti_1x300|Ultima|ONT|PacBio) *REQUIRED*
  -L,--library-type    <string> how inserts were prepared for --sequencing-platform (Nextera|TruSeq|Elevate|Ultima|Ligation|Rapid|HiFi) *REQUIRED*
  -I,--min-selected-size    <integer> smallest insert in bp selected BEFORE adapter ligation [0]
  -V,--selected-size-cv     <double> estimated coefficient of variation for --min-selected-size used to calculate min-allowed-size and max-allowed-size [0.15]
  -z,--min-allowed-size     <integer> use this value for min-allowed-size instead of calculating from --min-selected-size and --selected-size-cv [0]

restriction-enzyme:
  -e,--enzyme-name     <string> name of the restriction enzyme used to cleave genomic DNA (from shared/modules/REs/blunt_enzymes.csv if not NA) [NA]

pod5:
  -i,--pod5-dir        <string> directory or tar archive with one or more *.pod5 files from a single ONT run *REQUIRED*

dorado:
  -A,--dorado-version  <string> Dorado version to be downloaded and used, e.g., 1.3.0-linux-x64 [1.3.0-linux-x64]
  -M,--ont-model-complex    <string> ONT main basecalling model complex = (fast|hac|sup)[__AT_SYMBOL__(version|latest)] [sup]
  -b,--modified-base-models <string> optional comma-separated list of modified base models = modification[__AT_SYMBOL__(version|latest)]] 
  -5,--pod5-buffer     <string> whether to use shared memory (shm) or --tmp-dir (tmp) for pod5 buffering [shm]
  -Z,--pod5-buffer-size     <string> maximum size of pod5-buffer in Gb used when auto-batching POD5 files, e.g., 50G [60G]
  -y,--min-pod5-size   <string> minimum size of a pod5 file in Mb to use when auto-batching POD5 files, e.g., 100M [250M]
  -D,--dorado-options  <string> additional options passed directly to `dorado basecaller` [none] 
  -F,--force-basecalling    <boolean> force basecalling of all reads, ignoring any existing ubam files 

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
