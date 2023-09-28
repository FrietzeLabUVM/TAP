#!/bin/bash

SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
#echo $SCRIPT_PATH

#parameters parsing
#1 get config param if specified
#2 if config, parse params
#3 override any config params with remaining params

in_param="$@"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help) cat "$SCRIPT_PATH"/help_msg.txt; exit 0; shift ;;
    esac
    shift
done

#check if config specified
set -- ${in_param}
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c|--config) cfg="$2"; if [ ! -f $cfg ]; then echo "cannnot find config file: \"$cfg\". quit!"; exit 1; fi; shift ;;
    esac
    shift
done
if [ -z "$cfg" ]; then echo "config file (-c|--config) is required."; exit 1; fi

jh=job_headers
jt=job_templates
dt=../deployed_job_scripts
force=false
paths_relative=true

#parse args specified in config file by lines starting with #CFG
if [ -n "$cfg" ]; then
  echo gathering parameters from config file "$cfg"
  args=$(cat $cfg | awk -v cfg_prefix="#CFG" -v ORS=" " '{if ($1 == cfg_prefix){$1 = ""; print $0}}')
  args="${args//\~/$HOME}"
  if [ ! -z "$args" ]; then
    echo $args
    set -- $args
    while [[ "$#" -gt 0 ]]; do
      case $1 in
        -c|--config) echo ignoring config file specified in config file.; shift ;;
        -f|--force) force=true; ;;
        --job_headers) jh="$2"; shift ;;
        --job_templates) jt="$2"; shift ;;
        --deployment_target) dt="$2"; shift ;;
        --paths_relative) paths_relative=$2; shift;;
        *) echo "Unknown parameter passed: $1"; cat "$SCRIPT_PATH"/help_msg.txt; exit 1 ;;
        esac
      shift
    done
  fi
fi

#parse remaining args from command
set -- ${in_param}
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c|--config) shift ;;
        -f|--force) force=true; ;;
        --job_headers) jh="$2"; shift ;;
        --job_templates) jt="$2"; shift ;;
        --deployment_target) dt="$2"; shift ;;
        --paths_relative) paths_relative=$2; shift;;
        *) echo "Unknown parameter passed: $1"; cat "$SCRIPT_PATH"/help_msg.txt; exit 1 ;;
    esac
    shift
done

if [ "$paths_relative" = "true" ]; then
    jh="$SCRIPT_PATH/$jh"
    jt="$SCRIPT_PATH/$jt"
    dt="$SCRIPT_PATH/$dt"
fi

echo deployment config file is "$cfg"
echo job_headers directory is "$jh"
echo job_templates directory is "$jt"
echo deployment_target directory is "$dt"
echo force is $force



if [ ! -d "$jh" ]; then echo "job_headers directory \($jh\) not found!"; exit 1; fi
if [ ! -d "$jh" ]; then echo "job_templates directory \($jt\) not found!"; exit 1; fi
if [ -d "$dt" ]; then 
  if [ $force != "true" ]; then
    echo "deployment_target directory \($dt\) already exists! Remove or rerun this script with -f|--force to automatically replace."; 
    echo "rm -r $dt"
    exit 1; 
  fi
fi

mkdir -p "$dt"

todo=$(cat "$cfg" | awk '/^[^#]/ {print $0 }')


for f1 in $todo; do
  #template_script,job_name,header,processors
  temp_script=$(echo $f1 | awk -v FS="," '{print $1}');
  job_name=$(echo $f1 | awk -v FS="," '{print $2}');
  header_file=$(echo $f1 | awk -v FS="," '{print $3}');
  n_proc=$(echo $f1 | awk -v FS="," '{print $4}');
  echo "  temp_script $temp_script"
  echo "  job_name $job_name"
  echo "  header_file $header_file"
  echo "  n_proc $n_proc"

done
