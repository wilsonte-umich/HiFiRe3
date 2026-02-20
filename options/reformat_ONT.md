## reformat ONT options

|Tool Suite|Pipeline|Action|
|---|---|---|
|HiFiRe3|reformat|ONT|

```
hf3 reformat ONT --help
```

```
reformat: Utilities for reformatting input read files (not required by most users)

ONT: reformat ONT unaligned bam files to compress them and to convert legacy files

read-file:
  -i,--read-file-dir   <string> directory with one or more FASTQ or unaligned BAM files with required tags [TASK_DIR/ubam] 
  -x,--read-file-prefix     <string> file name prefix identifying files in --read-file-dir to align [use all files] 
  -n,--read-number-format   <string> portion of file name that identifies paired read files ('x' is replaced with 1 or 2) [_Rx_]

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
