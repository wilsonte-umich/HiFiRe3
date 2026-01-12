## Job file templates

The `templates` folder carries YAML job file templates
you can adapt for running HiFiRe3 pipelines on your data.

Template names reflect the pipeline and action they
correspond to, e.g., `<pipeline>_<action>.yml`.

The HiFiRe3 CLI can run job files either 
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
