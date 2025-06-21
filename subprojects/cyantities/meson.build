# The blueprint Cyantities meson.build file.
#
# Author: Malte J. Ziebarth (mjz.science@fmvkb.de)
#
# Copyright (C) 2024 Malte J. Ziebarth
#
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by
# the European Commission - subsequent versions of the EUPL (the "Licence");
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#
# https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the Licence is distributed on an "AS IS" basis,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and
# limitations under the Licence.
#
#
# Usage
# =====
#
# Include this file in your subprojects as in this example,
# and use the dependency
#    cyantities_dep = dependency(
#        'cyantities',
#        fallback : ['cyantities', 'libcyantities_dep']
#    )
# in your main meson.build file to include and link to Cyantities.
#

project('cyantities_dependency', 'cpp', 'cython',
        default_options : ['optimization=3'])


# Get C++ compiler and Python environment:

cc = meson.get_compiler('cpp')
python = import('python').find_installation()
dep_py = python.dependency()


#
# Find the Cyantities headers here:
#
incpath_cyantities = run_command(
  'cyantitiesinclude',
  check : false
)
if incpath_cyantities.returncode() != 0
    message(incpath_cyantities.stdout().strip())
    message(incpath_cyantities.stderr().strip())
    error('Could not determine Cyantities include path.')
else
    incpath_cyantities = incpath_cyantities.stdout().strip()
endif

#
# Find the cyantities library:
#
lib_path_cyantities = run_command(
  'cyantitieslibdir',
  check : false
)
if lib_path_cyantities.returncode() != 0
    message(lib_path_cyantities.stdout().strip())
    message(lib_path_cyantities.stderr().strip())
    error('Could not determine Cyantities include path.')
else
    lib_path_cyantities = lib_path_cyantities.stdout().strip()
endif
libcyantities = cc.find_library(
    'cyantities',
    dirs : [lib_path_cyantities],
    has_headers : [
        'cyantities/boost.hpp',
        'cyantities/quantitywrap.hpp',
        'cyantities/unit.hpp'
    ],
    header_include_directories : include_directories(incpath_cyantities),
    required : true,
    static : true
)
message('lib_path_cyantities:', lib_path_cyantities)

libcyantities_dep = declare_dependency(
    dependencies : libcyantities,
    include_directories : include_directories(incpath_cyantities)
)