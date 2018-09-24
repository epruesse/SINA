#!/bin/bash
# get last argument
log_base="${@:(-1):1}"
log_base_noslash=test$(echo $log_base | tr / _)
exec "$@" \
     --logger=JUNIT,all,"${JUNIT_REPORT_PATH}${log_base_noslash}.xml" \
     --report_level=detailed > ${log_base}.log -- $TEST_ARGS


