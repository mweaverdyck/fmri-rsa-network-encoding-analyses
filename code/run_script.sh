#!/bin/bash
#
# runs qsub scripts

# load modules and functions
source funcs
run_cmd=""
scripts=()
subs=()

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

for arg in "$@"; do
  # test if this argument is in a subject format
  is_sub $arg
  if [[ ${isSub} -ne 0 ]]; then
    # is a subject
    subs=( ${subs[@]} $arg )
    echo "subject added to list: $arg"
  else
    if [[ $arg == $ALL ]] || [[ $arg == $NEW ]]; then
      # run all subjects
      subs=( $arg $2)
    elif [[ $arg == "full" ]]; then
      # run full first level analysis
      scripts=( qsub_reconall.sh qsub_fmriprep.sh qsub_glm.sh qsub_rsa.sh qsub_rsaregs.sh )
    else
      # not a subject
      if [[ $arg != qsub* ]] && [[ $arg != /*.sh ]]; then arg=qsub_"$arg".sh; fi
      if [[ ! -f $arg ]]; then
        echo "ERROR: $arg is not a subject or a script. Skipping"
      else
        scripts=( $scripts $arg )
        echo $arg 'script added to list:' "${scripts[@]}"
      fi
    fi
  fi
  shift
done
if [[ ${#scripts} -eq 0 ]]; then
  echo "Must include a script name!"
fi

# get first step
get_step_from_script ${scripts[0]}
echo $STEP

is_step_second_level ${STEP}
if [[ ${#subs} -eq 0 ]] && [[ $IS_LEVEL2 -eq 1 ]]; then
  subs=${ALL}
else
  # get subjects
  convert_sub_args -f numid -c "${subs[@]}" $STEP
  subs=( "${SUBS[@]}" )
fi

# go through each subject and launch each script dependent on the previous
for s in "${subs[@]}"; do
    i=1
    for script in "${scripts[@]}"; do
        # rename script to full qsub_*.sh name if not already
        if [[ $script != qsub* ]]; then script=qsub_"$script".sh; fi
        # if this script doesn't exist, throw error.
        if [[ ! -f $script ]]; then echo "ERROR: $script does not exist."; exit 1; fi
        # if script is feat, add runs
        if [[ $script == *feat* ]]; then runs=( {1..8} ); else runs=""; fi
        # go through each run (will run once if empty)
        for r in "${runs[@]}"; do
            # determine if launching job or running interactively
            if [[ $run_cmd  == '-r' ]]; then
                # run job interactively
                echo "Running interactive job: . $script $s"
                . "$script" $s $r
            else
                # if r is not empty, make string for jobname
                if [[ ! -z $r ]]; then r_str="-$r"; fi
                echo "Launching job: qsub $script $s"
                # set jobname to subject, run, dependency order (i), script name
                job_name="s${s}${r_str}.${i}_${script:5}"
                # launch job
                if [[ $i -eq 1 ]]; then
                    # if first script, launch without dependency
                    qsub -N "${job_name}" "$script" $s $r
                else
                    # if not first script, launch with dependency
                    qsub -N "${job_name}" -hold_jid ${prev_job_name} "$script" $s $r
                fi
            fi
        done
        # save previous job for dependency
        prev_job_name=$job_name
        # increase dependency order by 1
        ((i++))
    done
done
