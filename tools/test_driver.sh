#!/bin/bash
# get last argument
log_base="${@:(-1):1}"
log_base_noslash=test$(echo $log_base | tr / _)
exec "$@" \
     --logger=JUNIT,all,"${JUNIT_REPORT_PATH}${log_base_noslash}.xml" \
     --logger=HRF \
     --report_level=detailed \
     --no_color_output \
     -- $TEST_ARGS \
     > ${log_base}.log \


