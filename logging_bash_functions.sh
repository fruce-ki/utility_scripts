# Log current comment command. Use like this:
# Log only if successful: <command> && logme
# Log regardless        : <command> ; logme
logme() {
    comm=$(fc -ln -1 | sed 's/^[[:space:]]*//')
    echo $(date "+%Y-%m-%d %H:%M:%S")$'\t\t'$(echo ${comm} | sed 's/&&logme/\&\& logme/') >> ./commands.log
    echo $(date "+%Y-%m-%d %H:%M:%S")$'\t\t'$(echo ${comm} | sed 's/[[:space:]]*[&;]*[[:space:]]*logme[[:space:]]*$//') >> ./subcommands.log
}

# Retroactively log the previous command.
loglast() {
    comm=$(fc -ln -2 | head -n 1 | sed 's/^[[:space:]]*//')
    echo $(date "+%Y-%m-%d %H:%M:%S")$'\t\t'${comm}' && logme' >> ./commands.log
    echo $(date "+%Y-%m-%d %H:%M:%S")$'\t\t'${comm} >> ./subcommands.log
}
