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
        -h|--help) cat $SCRIPT_PATH/help_msg.txt; exit 0; shift ;;
    esac
    shift
done

#check if config specified
set -- ${in_param}
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c|--config) cfg="$2"; if [ ! -f $cfg ]; then echo cannnot find config file $cfg. quit!; exit 1; fi; shift ;;
    esac
    shift
done

jh=$SCRIPT_PATH/job_headers
jt=$SCRIPT_PATH/job_templates
dt=$SCRIPT_PATH/../deployed_job_scripts
force=false

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
        --job_headers) jh=job_headers; shift ;;
        --job_templates) jt=job_templates; shift ;;
        --deployment_target) dt=../deployed_job_scripts; shift ;;
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
        --job_headers) jh=job_headers; shift ;;
        --job_templates) jt=job_templates; shift ;;
        --deployment_target) dt=../deployed_job_scripts; shift ;;
        *) echo "Unknown parameter passed: $1"; cat "$SCRIPT_PATH"/help_msg.txt; exit 1 ;;
    esac
    shift
done

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