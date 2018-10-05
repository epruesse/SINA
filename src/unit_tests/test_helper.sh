set -o pipefail

# copy FDs
exec 5<&0
exec 6>&1
exec 7>&2

test_num=0
test_count=$(grep -c "begin_test" "$0")
echo "1..$test_count"

test_helper_tmpdirs_=()
trap test_helper_cleanup_ EXIT

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
    output=$(eval "$*" | tee /dev/fd/6)
    exit_code=$?
}
capture_stderr() {
    echo "Running \"$*"\"
    output=$(eval "$*" 2>&1 >&6 | tee /dev/fd/7)
    exit_code=$?
}
capture_stdouterr() {
    echo "Running \"$*"\"
    output=$(eval "$*" 2>&1 | tee /dev/fd/6)
    exit_code=$?
}

begin_test() {
    let test_num++
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
	test_err="$test_err
# command exited with $exit_code 
"
    fi
}

assert_exit_fail() {
    if [ $exit_code -eq 0 ]; then
	test_err="$test_err
# command exited with $exit_code
"
    fi
}

assert_output_contains() {
    if echo $output | grep -q "$1"; then
	:
    else
	test_err="$test_err
# command output did not include '$1'
"
    fi
}