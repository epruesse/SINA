#!/bin/bash
# get last argument
log_base="${@:(-1):1}"

exec "$@" \
     --logger=JUNIT,all,"${JUNIT_REPORT_PATH}${log_base}.xml" \
     --report_level=detailed > ${log_base}.log


