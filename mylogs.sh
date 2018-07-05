#!/bin/bash

# Log a command and then execute it.
if [ "$1" = '-e' ]; then
    comm=$(echo "$@" | sed 's/[$*\|"<>'"'"']/\\&/g')
    echo $(date "+%Y-%m-%d %H:%M:%S")$'\t\t'"mylogs.sh -e "$comm >> ./commands.log
    echo $(date "+%Y-%m-%d %H:%M:%S")$'\t\t'"${@:2}" >> ./subcommands.log
    ${@:2}
fi

# Log a message.
if [ "$1" = '-m' ]; then
comm=$(echo "${@:2}" | sed 's/[$*\|"<>'"'"']/\\&/g')
echo $(date "+%Y-%m-%d %H:%M:%S")$'\t\t'$comm >> ./messages.log
fi
