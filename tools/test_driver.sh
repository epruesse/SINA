#!/bin/bash
# get last argument
#env

log_base="${@:(-1):1}"
export LD_LIBRARY_PATH=$ARBHOME/lib
#exec "$@" --log_level=all --output_format=xml --report_level=detailed \
#          --report_sink="${log_base}.report.xml" \
#          --log_sink="${log_base}.log.xml"
#exec "$@" --log_level=all --report_level=no --output_format=xml > ${log_base}.log.xml
exec "$@" --log_level=all --report_level=detailed > ${log_base}.log


