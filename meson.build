# Adapted from Cyantities/exaples/gravity/meson.build
project(
    'etascatgen',
    'cpp',
    'cython',
    default_options : ['optimization=3']
)


#
# Dependencies:
#
boost_dep = dependency('boost')

cyantities_dep = dependency(
    'cyantities',
    fallback : ['cyantities', 'libcyantities_dep']
)



#
# Find the NumPy headers here:
#
python = import('python').find_installation()
dep_py = python.dependency()

incpath_np = run_command(
  'gravelspoon',
  check : false
)
if incpath_np.returncode() != 0
    message(incpath_np.stdout().strip())
    message(incpath_np.stderr().strip())
    error('Could not determine NumPy include path.')
else
    incpath_np = incpath_np.stdout().strip()
endif


#
# Include directory:
#
incdir = include_directories(['cpp/include', incpath_np])

add_project_arguments(['-std=c++23','-Wall'], language : 'cpp')


#
# The C++ sources:
#
libetascatgen = static_library(
    'etascatgen',
    [
        'cpp/src/catgen_M_t.cpp'
    ],
    include_directories: incdir,
    dependencies: [boost_dep, cyantities_dep]
)

#
# Finally compile the extension module:
#
python.extension_module(
    'backend',
    ['etascatgen/backend.pyx'],
    dependencies : [dep_py, cyantities_dep],
    include_directories : [incdir],
    link_with: libetascatgen,
    override_options : ['cython_language=cpp']
)