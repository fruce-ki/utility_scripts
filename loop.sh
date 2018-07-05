#!/bin/bash

function usage() {
    echo "Usage:"
    echo "      $0 [-f NUMBER_FROM -t NUMBER_TO] | [-v STRING_VALUE [-v STRING_VALUE] [...]]] -c 'COMMAND_STRING'"
    echo "Provide either the start and end of a numeric range, or multiple string literals, over which the command is to be iterated. Values will be substituted for '{task_id}' occurances in the command."
    exit 1
}

# Parse options.
while getopts 'f:t:v:c:h:' flag; do
  case "${flag}" in
    f) from="${OPTARG}" ;;
    t) to="${OPTARG}" ;;
    v) values+=("${OPTARG}") ;;
    c) command="${OPTARGS}";;
    h) usage ;;
  esac
done
shift "$((OPTIND-1))"

if [ ! -z $from ]; then
    if [ ! -z $to ]; then
	values=($(seq $from $to))
    fi
fi

command="$@"

for value in "${values[@]}"; do
	echo "${value}..."
	${command//'{val}'/$value}
done
