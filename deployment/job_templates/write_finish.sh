#!/bin/bash

date > ${1}.finish
echo if this file exists and is newer than start, all jobs have finished, possibly with errors >> ${1}.finish
