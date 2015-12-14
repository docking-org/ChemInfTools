#!/bin/bash --norc
#$ -S /bin/bash
#$ -cwd

set -e


if [ -z "${JOB_ID}" ]; then
	export JOB_ID=$$
fi

if [ -z "${SGE_TASK_ID}" ]; then
	export SGE_TASK_ID=1
fi

if [ -z "${__MR_SGE_JOB_DIR}" ]; then
	JOB_DIR=`pwd`
	if [ -e "${JOB_DIR}/.qsub-map-logs" ] ; then
		JOB_DIR="$( dirname "${JOB_DIR}" )"
	fi
else
	JOB_DIR="${__MR_SGE_JOB_DIR}"
fi

INPUTS="${JOB_DIR}/outputs"
OUTPUT="${JOB_DIR}/merged"

if [ -z "${SCRATCH_DIR}" ]; then
	SCRATCH_DIR=/scratch
	if [ ! -d "${SCRATCH_DIR}" ]; then
	    SCRATCH_DIR=/tmp
	fi
fi

if [ -z "${TASK_DIR}" ]; then
	TASK_DIR="${SCRATCH_DIR}/$( whoami )/${JOB_ID}"
fi

if [ ! -z "${__MR_SGE_TASK_SETUP_SCRIPT}" ]; then
	source "${__MR_SGE_TASK_SETUP_SCRIPT}" 1>&2
fi

echo "Input:   ${INPUTS}" 1>&2
echo "Working: $( hostname ):${TASK_DIR}" 1>&2
echo "Output:  ${OUTPUT}" 1>&2

mkdir -pv "${TASK_DIR}" 1>&2
pushd "${TASK_DIR}" 1>&2
set +e
"$@" "${INPUTS}"
set -e
popd 1>&2
mv -v "${TASK_DIR}" "${OUTPUT}" 1>&2
