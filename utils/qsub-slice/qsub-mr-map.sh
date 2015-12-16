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

if [ -z "${__MR_SGE_OUTPUT_METHOD}" ] ; then
	OUTPUT_METHOD="/"
else
	OUTPUT_METHOD="${__MR_SGE_OUTPUT_METHOD}"
fi

INPUTS="${JOB_DIR}/inputs"
OUTPUTS="${JOB_DIR}/outputs"
TASK_FILE=$( ls "${INPUTS}" | awk "NR == ${SGE_TASK_ID}" )
TASK_NAME="$( basename "${TASK_FILE}" )"
TASK_INPUT="${INPUTS}/${TASK_FILE}"
TASK_OUTPUT="${OUTPUTS}/${TASK_FILE}"


if [ ! -z "${__MR_SGE_MAP_ITERATE}" ] ; then
	TASK_COMMAND=( xargs -L 1 -a "${TASK_INPUT}" "$@" )
else
	TASK_COMMAND=( "$@" )
fi

if [ -x "${SCRATCH_DIR}" ] ; then
	SCRATCH_DIR="$( "${SCRATCH_DIR}" )"
elif [ -z "${SCRATCH_DIR}" ] ; then
	SCRATCH_DIR=/scratch
	if [ ! -d "${SCRATCH_DIR}" ]; then
	    SCRATCH_DIR=/tmp
	fi
	SCRATCH_DIR="${SCRATCH_DIR}/$( whoami )/$$"
fi

if [ -z "${TASK_DIR}" ]; then
	TASK_DIR="${SCRATCH_DIR}/${TASK_NAME}"
fi

if [ "${OUTPUT_METHOD}" == "-" ] ; then
	TASK_RESULT="${TASK_DIR}/__QSUB_MR_RESULT"
elif [ "${OUTPUT_METHOD}" == "/" ] ; then
	TASK_RESULT="${TASK_DIR}"
elif [ -z "${OUTPUT_METHOD}" ] ; then
	TASK_RESULT=""
else 
	TASK_RESULT="${TASK_DIR}/${OUTPUT_METHOD}"
fi

if [ ! -z "${__MR_SGE_TASK_SETUP_SCRIPT}" ]; then
	source "${__MR_SGE_TASK_SETUP_SCRIPT}" 1>&2
fi

echo "Input:   ${TASK_INPUT}" 1>&2
echo "Working: $( hostname ):${TASK_DIR}" 1>&2
echo "Result:  ${TASK_RESULT}" 1>&2
echo "Output:  ${TASK_OUTPUT}" 1>&2
echo "Task:    ${TASK_COMMAND[@]}" 1>&2

mkdir -pv "${TASK_DIR}" 1>&2
pushd "${TASK_DIR}" 1>&2
set +e
if [ "${OUTPUT_METHOD}" == "-" ] ; then
	"${TASK_COMMAND[@]}" "${TASK_INPUT}" > "${TASK_RESULT}"
else
	"${TASK_COMMAND[@]}" "${TASK_INPUT}"
fi
set -e
popd 1>&2

mv -v "${TASK_RESULT}" "${TASK_OUTPUT}" 1>&2
rm -rfv "${TASK_DIR}" 1>&2
