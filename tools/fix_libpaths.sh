#!/bin/bash
set -e

UNAME=$(uname)

logrun() {
    echo "$@"
    "$@"
}

get_rpaths_Darwin() {
    otool -l "$1" | grep LC_RPATH -A2 |grep path | awk '{print $2}'
}
get_rpaths_Linux() {
    patchelf --print-rpath "$1" | tr : '\n'
}
get_rpaths() {
    get_rpaths_$UNAME "$@"
}

get_needed_Darwin() {
    otool -L "$1" | tail -n +2 | awk '{print $1}'
}
get_needed_Linux() {
    patchelf --print-needed "$1"
}
get_needed() {
    get_needed_$UNAME "$@"
}

remove_rpaths_Darwin() {
    for rpath in `get_rpaths_Darwin $1`; do
	logrun install_name_tool -delete_rpath $rpath "$1"
    done
}
remove_rpaths_Linux() {
    logrun patchelf --remove-rpath "$1"
}
remove_rpaths() {
    remove_rpaths_$UNAME "$@"
}

set_rpath_Darwin() {
    remove_rpath_Darwin "$1"
    logrun install_name_tool -add_rpath "$2" "$1"
}
set_rpath_Linux() {
    logrun patchelf --set-rpath "$2" "$1"
}
set_rpath() {
    set_rpath_$UNAME "$@"
}

copy_libs_needed() {
    rpaths=`get_rpaths $1`
    dtneed=`get_needed $1 | grep -v /usr/lib`
    dest=$(dirname $1)/../lib
    for lib in $dtneed; do
	name="${lib#@rpath/}"
	if test -e "$dest/$name"; then
	    continue;
	fi
	source=
	for path in "" $rpaths /usr/local/lib; do
	    if test -e "$path/$name"; then
		source="$path/$name"
		break;
	    fi
	done
	if test -z "$source"; then
	    echo Missing $name
	else
	    echo "$source -> $dest/$name"
	    cat $source > $dest/$name
	fi
    done
}
		    

set_rpath_Darwin() {
    remove_rpaths_Darwin "$1"
    add_rpath_Darwin "$1" "$2"
}

make_rpath_relative_Darwin() {
    set_rpath_Darwin "$1" @loader_path/../lib
}

make_absolute_path_relative_Darwin() {
    for lib in `get_needed_Darwin $1`; do
	if test ${lib:0:1} = "/"; then
	    name=$(basename $lib)
	    dir="$(dirname $lib)"
	    if test "$dir" != "/usr/lib"; then
		logrun install_name_tool -change $lib @rpath/$name -add_rpath $dir "$1"
	    fi
	fi
    done
}

if test "$0" != "-bash"; then
    if test -n "$1"; then
	make_absolute_path_relative_Darwin "$1"
	copy_libs_needed "$1"
	make_rpath_relative_Darwin "$1"
    fi
fi
