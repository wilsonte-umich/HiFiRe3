## basecall PacBio options

|Tool Suite|Pipeline|Action|
|---|---|---|
|HiFiRe3|basecall|PacBio|

```
hf3 basecall PacBio --help
```

```
basecall: Perform basecalling on input read data from the sequencing platform

PacBio: derive a read consensus from PacBioStrand strand consensuses generated in by-strand mode

pacbio:
  -H,--hifi-bam-file   <string> path to a PacBioStrand BAM file with hifi reads generated in by-strand mode *REQUIRED*
  -F,--fail-bam-file   <string> path to a PacBioStrand BAM file with fail reads generated in by-strand mode *REQUIRED*
  -X,--force-basecalling    <boolean> force basecalling of all reads, ignoring any existing ubam files 

genome:
  -g,--genome          <string> name of the reference genome to use (e.g., hg38) *REQUIRED*
  -G,--genome-dir      <string> directory with indexed <--genome>.fa or genome.fa file [.../mdi/resources/genomes/<--genome>] 
  --use-all-chroms     <boolean> use all chromosomes as they are found in genome fasta file without filtering 
  --is-composite-genome     <boolean> <--genome> is a composite of >1 reference (e.g., hs1_dm6 with cross-species spike-in) 

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
