set -Eo pipefail

# copy FDs
exec 5<&0
exec 6>&1
exec 7>&2

test_num=0
test_count=$(grep -c "begin_test" "$0")
echo "1..$test_count"

test_helper_tmpdirs_=()
trap test_helper_cleanup_ EXIT
script_error_disable_=
trap script_error ERR


script_error() {
    if [ -z "$script_error_disable_" ]; then
	n="${#FUNCNAME[@]}"
	test_err="${test_err}Error in bash, callstack:
"
	let n--;
	for ((; n>0; --n)); do
	    let m=n-1
	    test_err="${test_err} $n:  ${BASH_SOURCE[$n]}:${BASH_LINENO[$m]} in ${FUNCNAME[$n]}()
"
	done
    fi
}

maketmpdir() {
    local tmp
    tmp=$(mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir')
    test_helper_tmpdirs_+=($tmp)
    echo $tmp
}

test_helper_cleanup_() {
    for dir in "${test_helper_tmpdirs_[@]}"; do
	echo "Removing $dir"
	rm -rf "$dir"
    done
}

capture_stdout() {
    echo "Running \"$*"\"
    script_error_disable_=1
    output=$(eval "$*" | tee /dev/fd/6)
    exit_code=$?
    script_error_disable_=
}
capture_stderr() {
    echo "Running \"$*"\"
    script_error_disable_=1
    output=$(eval "$*" 2>&1 >&6 | tee /dev/fd/7)
    exit_code=$?
    script_error_disable_=
}
capture_stdouterr() {
    echo "Running \"$*"\"
    script_error_disable_=1
    output=$(eval "$*" 2>&1 | tee /dev/fd/6)
    exit_code=$?
    script_error_disable_=
}

begin_test() {
    let ++test_num
    test_name=$1
    test_err=
}

end_test() {
    if [ -n "$test_err" ]; then
	echo -n "not "
    fi
    echo "ok $test_num - $test_name"
    echo "$test_err"
}

assert_exit_success() {
    if [ $exit_code -ne 0 ]; then
	test_err="${test_err}# command exited with $exit_code
"
    fi
}

assert_exit_failure() {
    if [ $exit_code -eq 0 ]; then
	test_err="${test_err}# command exited with $exit_code
"
    fi
}

assert_output_contains() {
    if echo "$output" | grep -q "$1"; then
	:
    else
	test_err="${test_err}# command output did not include '$1'
"
    fi
}

assert_output_not_contains() {
    if echo "$output" | grep -q "$1"; then
	test_err="${test_err}# command output did not include '$1'
"
    fi
}

assert_output_count() {
    if test $(echo "$output" | grep -c "$1") $2 $3; then
	:
    else
	test_err="${test_err}# command output did not match $@
"
    fi
}
