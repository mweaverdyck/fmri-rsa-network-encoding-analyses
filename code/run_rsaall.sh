#!/bin/bash
#
# runs qsub scripts

# load modules and functions
source funcs
run_cmd=""
script="qsub_rsaall.sh"

while test $# -gt 0; do
  case "$1" in
    -h|--help)
      echo "run_script.sh launches individual jobs for each subject for each script (or runs in sequence in interactive session)"
      echo "Usage: bash run_script.sh [-r] <scriptname> <sub> [<sub>]..."
      echo "       [-r]           runs script in interactive session (should qrsh first)"
      echo "       <scriptname>   script to run or 'full' for all first-level scripts"
      echo "                      must match format for a script that begins with 'qsub_'"
      echo "                          qsub_scriptname.sh"
      echo "                          scriptname=qsub_*.sh"
      echo "                          scriptname=full (will run full firstlevel pipeline)"
      echo "       sub [sub]...   subject number or subject id (may include as many as you'd like or 'all')"
      exit 0
      ;;
    -r)
      run_cmd='-r'
      shift
      ;;
    *)
      break
      ;;
  esac
done

# get subjects
convert_sub_args -f numid -c "$@" RSA
subs=( "${SUBS[@]}" )

# go through each subject and launch each script dependent on the previous
for s in "${subs[@]}"; do
    # determine if launching job or running interactively
    if [[ $run_cmd  == '-r' ]]; then
        # run job interactively
        echo "Running interactive job: . $script $s"
        . "$script" $s
    else
      for i in {1..2}; do
        echo "Launching job: qsub $script $s"
        # set jobname to subject, run, dependency order (i), script name
        job_name="s${s}.${i}_${script:5}"
        # launch job
        if [[ $i -eq 1 ]]; then
            # if first script, launch without dependency
            qsub -N "${job_name}" "$script" $s
        else
            # if not first script, launch with dependency
            qsub -N "${job_name}" -hold_jid ${prev_job_name} "$script" $s
        fi
        # save previous job for dependency
        prev_job_name=${job_name}
      done
    fi
done
