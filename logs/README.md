## Job log reports

The `logs` folder carries text files containing the output
of the follow commands applied to example job files.

```sh
hf3 myJob.yml status # show the state of all submitted jobs
hf3 myJob.yml report # show a job log report
hf3 myJob.yml ls     # show the contents of a job's output diretory
```

Log file names reflect the pipeline and action they
correspond to, e.g., `<pipeline>_<action>.yml`.

Compare these log reports to your own job logs when
troubleshooting your installation and data analysis.
