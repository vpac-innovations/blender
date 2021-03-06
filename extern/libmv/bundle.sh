#!/bin/sh

if [ "x$1" = "x--i-really-know-what-im-doing" ] ; then
  echo Proceeding as requested by command line ...
else
  echo "*** Please run again with --i-really-know-what-im-doing ..."
  exit 1
fi

BRANCH="master"

repo="git://git.blender.org/libmv.git"
tmp=`mktemp -d`

git clone -b $BRANCH $repo $tmp/libmv

git --git-dir $tmp/libmv/.git --work-tree $tmp/libmv log -n 50 > ChangeLog

find libmv -type f -not -iwholename '*.svn*' -exec rm -rf {} \;
find third_party -type f -not -iwholename '*.svn*' -not -iwholename '*third_party/ceres*' \
    -not -iwholename '*third_party/SConscript*' -not -iwholename '*third_party/CMakeLists.txt*' \
    -exec rm -rf {} \;

cat "files.txt" | while read f; do
  mkdir -p `dirname $f`
  cp $tmp/libmv/src/$f $f
done

rm -rf $tmp

chmod 664 ./third_party/glog/src/windows/*.cc ./third_party/glog/src/windows/*.h ./third_party/glog/src/windows/glog/*.h

sources=`find ./libmv -type f -iname '*.cc' -or -iname '*.cpp' -or -iname '*.c' | grep -v _test.cc | grep -v test_data_sets | sed -r 's/^\.\//\t\t/' | sort -d`
headers=`find ./libmv -type f -iname '*.h' | grep -v test_data_sets | sed -r 's/^\.\//\t\t/' | sort -d`

third_sources=`find ./third_party -type f -iname '*.cc' -or -iname '*.cpp' -or -iname '*.c' | grep -v glog | grep -v gflags | grep -v ceres | sed -r 's/^\.\//\t\t/' | sort -d`
third_headers=`find ./third_party -type f -iname '*.h' | grep -v glog | grep -v gflags | grep -v ceres | sed -r 's/^\.\//\t\t/' | sort -d`

third_glog_sources=`find ./third_party -type f -iname '*.cc' -or -iname '*.cpp' -or -iname '*.c' | grep glog | grep -v windows | sed -r 's/^\.\//\t\t\t/' | sort -d`
third_glog_headers=`find ./third_party -type f -iname '*.h' | grep glog | grep -v windows | sed -r 's/^\.\//\t\t\t/' | sort -d`

third_gflags_sources=`find ./third_party -type f -iname '*.cc' -or -iname '*.cpp' -or -iname '*.c' | grep gflags | grep -v windows | sed -r 's/^\.\//\t\t/' | sort -d`
third_gflags_headers=`find ./third_party -type f -iname '*.h' | grep gflags | grep -v windows | sed -r 's/^\.\//\t\t/' | sort -d`

tests=`find ./libmv -type f -iname '*_test.cc' | sort -d | awk ' { name=gensub(".*/([A-Za-z_]+)_test.cc", "\\\\1", $1); printf("\t\tBLENDER_SRC_GTEST(\"libmv_%s\" \"%s\" \"libmv_test_dataset;extern_libmv;extern_ceres\")\n", name, $1) } '`

src_dir=`find ./libmv -type f -iname '*.cc' -exec dirname {} \; -or -iname '*.cpp' -exec dirname {} \; -or -iname '*.c' -exec dirname {} \; | sed -r 's/^\.\//\t\t/' | sort -d | uniq`
src_third_dir=`find ./third_party -type f -iname '*.cc' -exec dirname {} \; -or -iname '*.cpp' -exec dirname {} \; -or -iname '*.c' -exec dirname {} \;  | grep -v ceres | sed -r 's/^\.\//\t\t/'  | sort -d | uniq`
src=""
win_src=""
for x in $src_dir $src_third_dir; do
  t=""

  if test  `echo "$x" | grep -c glog ` -eq 1; then
    continue;
  fi

  if stat $x/*.cpp > /dev/null 2>&1; then
    t="    src += env.Glob('`echo $x'/*.cpp'`')"
  fi

  if stat $x/*.c > /dev/null 2>&1; then
    if [ -z "$t" ]; then
      t="    src += env.Glob('`echo $x'/*.c'`')"
    else
      t="$t + env.Glob('`echo $x'/*.c'`')"
    fi
  fi

  if stat $x/*.cc > /dev/null 2>&1; then
    if [ -z "$t" ]; then
      t="    src += env.Glob('`echo $x'/*.cc'`')"
    else
      t="$t + env.Glob('`echo $x'/*.cc'`')"
    fi
  fi

  if test `echo $x | grep -c "windows\|gflags" ` -eq 0; then
    if [ -z "$src" ]; then
      src=$t
    else
      src=`echo "$src\n$t"`
    fi
  else
    if [ -z "$win_src" ]; then
      win_src=`echo "    $t"`
    else
      win_src=`echo "$win_src\n    $t"`
    fi
  fi
done

cat > CMakeLists.txt << EOF
# ***** BEGIN GPL LICENSE BLOCK *****
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# The Original Code is Copyright (C) 2011, Blender Foundation
# All rights reserved.
#
# Contributor(s): Blender Foundation,
#                 Sergey Sharybin
#
# ***** END GPL LICENSE BLOCK *****

# NOTE: This file is automatically generated by bundle.sh script
#       If you're doing changes in this file, please update template
#       in that script too

set(INC
	.
)

set(INC_SYS
)

set(SRC
	libmv-capi.h
)

if(WITH_LIBMV OR WITH_GTESTS OR (WITH_CYCLES AND WITH_CYCLES_LOGGING))
	list(APPEND INC
		third_party/gflags
		third_party/gflags/gflags
		third_party/glog/src
		third_party/ceres/include
		third_party/ceres/config
		../../intern/guardedalloc
	)

	list(APPEND INC_SYS
		../Eigen3
		\${PNG_INCLUDE_DIRS}
		\${ZLIB_INCLUDE_DIRS}
	)

	if(WIN32)
		list(APPEND INC
			third_party/glog/src/windows
		)

		if(NOT MINGW)
			list(APPEND INC
				third_party/msinttypes
			)
		endif()
	endif()

	add_definitions(
		-DWITH_LIBMV_GUARDED_ALLOC
		-DGOOGLE_GLOG_DLL_DECL=
		-DLIBMV_NO_FAST_DETECTOR=
	)
endif()

if(WITH_LIBMV)
	TEST_SHARED_PTR_SUPPORT()
	if(SHARED_PTR_FOUND)
		if(SHARED_PTR_TR1_MEMORY_HEADER)
			add_definitions(-DCERES_TR1_MEMORY_HEADER)
		endif()
		if(SHARED_PTR_TR1_NAMESPACE)
			add_definitions(-DCERES_TR1_SHARED_PTR)
		endif()
	else()
		message(FATAL_ERROR "Unable to find shared_ptr.")
	endif()

	list(APPEND SRC
		intern/autotrack.cc
		intern/camera_intrinsics.cc
		intern/detector.cc
		intern/frame_accessor.cc
		intern/homography.cc
		intern/image.cc
		intern/logging.cc
		intern/reconstruction.cc
		intern/track_region.cc
		intern/tracks.cc
		intern/tracksN.cc
${sources}
${third_sources}

		intern/autotrack.h
		intern/camera_intrinsics.h
		intern/detector.h
		intern/frame_accessor.h
		intern/homography.h
		intern/image.h
		intern/logging.h
		intern/reconstruction.h
		intern/track_region.h
		intern/tracks.h
		intern/tracksN.h
${headers}

${third_headers}
	)


	if(WITH_GTESTS)
		blender_add_lib(libmv_test_dataset "./libmv/multiview/test_data_sets.cc" "${INC}" "${INC_SYS}")

${tests}
	endif()
else()
	list(APPEND SRC
		intern/stub.cc
	)
endif()

blender_add_lib(extern_libmv "\${SRC}" "\${INC}" "\${INC_SYS}")

if(WITH_LIBMV)
	add_subdirectory(third_party)
endif()

# make GLog a separate target, so it can be used for gtest as well.
if(WITH_LIBMV OR WITH_GTESTS OR (WITH_CYCLES AND WITH_CYCLES_LOGGING))
	# We compile GLog together with GFlag so we don't worry about
	# adding extra lib to linker.
	set(GLOG_SRC
${third_gflags_sources}

${third_gflags_headers}
	)

	if(WIN32)
		list(APPEND GLOG_SRC
			third_party/glog/src/logging.cc
			third_party/glog/src/raw_logging.cc
			third_party/glog/src/utilities.cc
			third_party/glog/src/vlog_is_on.cc
			third_party/glog/src/windows/port.cc

			third_party/glog/src/utilities.h
			third_party/glog/src/stacktrace_generic-inl.h
			third_party/glog/src/stacktrace.h
			third_party/glog/src/stacktrace_x86_64-inl.h
			third_party/glog/src/base/googleinit.h
			third_party/glog/src/base/mutex.h
			third_party/glog/src/base/commandlineflags.h
			third_party/glog/src/stacktrace_powerpc-inl.h
			third_party/glog/src/stacktrace_x86-inl.h
			third_party/glog/src/config.h
			third_party/glog/src/stacktrace_libunwind-inl.h
			third_party/glog/src/windows/glog/raw_logging.h
			third_party/glog/src/windows/glog/vlog_is_on.h
			third_party/glog/src/windows/glog/logging.h
			third_party/glog/src/windows/glog/log_severity.h
			third_party/glog/src/windows/port.h
			third_party/glog/src/windows/config.h

			third_party/gflags/windows_port.cc
			third_party/gflags/windows_port.h
		)
	else()
		list(APPEND GLOG_SRC
${third_glog_sources}

${third_glog_headers}
		)
	endif()

	blender_add_lib(extern_glog "\${GLOG_SRC}" "\${INC}" "\${INC_SYS}")
endif()
EOF

cat > SConscript << EOF
#!/usr/bin/python

# NOTE: This file is automatically generated by bundle.sh script
#       If you're doing changes in this file, please update template
#       in that script too

import sys
import os

Import('env')

defs = []
incs = '.'

if env['WITH_BF_LIBMV'] or (env['WITH_BF_CYCLES'] and env['WITH_BF_CYCLES_LOGGING']):
    defs.append('GOOGLE_GLOG_DLL_DECL=')
    defs.append('WITH_LIBMV_GUARDED_ALLOC')
    defs.append('LIBMV_NO_FAST_DETECTOR')

    incs += ' ../Eigen3 third_party/gflags third_party/gflags/gflags third_party/glog/src third_party/ceres/include third_party/ceres/config ../../intern/guardedalloc'
    incs += ' ' + env['BF_PNG_INC']
    incs += ' ' + env['BF_ZLIB_INC']

    if env['OURPLATFORM'] in ('win32-vc', 'win32-mingw', 'linuxcross', 'win64-vc', 'win64-mingw'):
        incs += ' ./third_party/glog/src/windows ./third_party/glog/src/windows/glog'
        if env['OURPLATFORM'] in ('win32-vc', 'win64-vc'):
            incs += ' ./third_party/msinttypes'
    else:
        incs += ' ./third_party/glog/src'

if env['WITH_BF_LIBMV']:
    if not env['WITH_SHARED_PTR_SUPPORT']:
        print("-- Unable to find shared_ptr which is required for compilation.")
        exit(1)

    if env['SHARED_PTR_HEADER'] == 'tr1/memory':
        defs.append('CERES_TR1_MEMORY_HEADER')
    if env['SHARED_PTR_NAMESPACE'] == 'std::tr1':
        defs.append('CERES_TR1_SHARED_PTR')

    src = env.Glob('intern/*.cc')
    src.remove('intern' + os.sep + 'stub.cc')
$src
else:
    src = env.Glob("intern/stub.cc")

src = [src for src in src if src.find('_test.cc') == -1]

env.BlenderLib(libname = 'extern_libmv', sources=src, includes=Split(incs), defines=defs, libtype=['extern', 'player'], priority=[20,137])

if env['WITH_BF_LIBMV'] or (env['WITH_BF_CYCLES'] and env['WITH_BF_CYCLES_LOGGING']):
    glog_src = []
    glog_src += env.Glob("third_party/gflags/*.cc")
    if env['OURPLATFORM'] in ('win32-vc', 'win32-mingw', 'linuxcross', 'win64-vc', 'win64-mingw'):
        glog_src += ['./third_party/glog/src/logging.cc', './third_party/glog/src/raw_logging.cc', './third_party/glog/src/utilities.cc', './third_party/glog/src/vlog_is_on.cc']
        glog_src += ['./third_party/glog/src/windows/port.cc']
    else:
        glog_src.remove('third_party/gflags/windows_port.cc')
        glog_src += env.Glob("third_party/glog/src/*.cc")

    env.BlenderLib(libname = 'extern_glog', sources=glog_src, includes=Split(incs), defines=defs, libtype=['extern', 'player'], priority=[20,137])

if env['WITH_BF_LIBMV']:
    SConscript(['third_party/SConscript'])
EOF
