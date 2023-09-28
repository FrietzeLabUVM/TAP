#!/bin/bash

date > ${1}.complete
echo if this file exists and is newer than start, then all jobs have finished without any errors.  the output is complete and should not be overwritten.  >> ${1}.complete
