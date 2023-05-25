#!/bin/sh

set -e

waitfile="$1"
shift
cmd="$@"

until test -e $waitfile
do
    >$2 echo "Waiting for file [$waitfile]."
    sleep 1
done

>&2 echo "Found file [$waitfile]."
exec $cmd
