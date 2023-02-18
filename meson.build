project('sim_geom', 'cpp',
  default_options : [ 'warning_level=3', 'buildtype=debugoptimized', 'cpp_std=c++11' ]
)

#find dependencies
sim_deps = [dependency('meep'), dependency('hdf5', language: 'cpp')]
inc_dirs = [include_directories('/usr/include/doctest')]

#set executable names
app_name = 'sim_geom'
tst_name = 'tst'

#source files
app_srcs = files('src/cgs_read.cpp', 'src/cgs.cpp', 'src/data_utils.cpp', 'src/disp.cpp', 'src/main.cpp')
tst_srcs = files('src/cgs_read.cpp', 'src/cgs.cpp', 'src/data_utils.cpp', 'src/disp.cpp', 'src/main_test.cpp')

#check whether 32 bit or 64 bit precision for saved fields is used
app_args = []
if get_option('float_bits') == 64
  app_args += '-DSTO_PREC_32=0'
else
  app_args += '-DSTO_PREC_32=1'
endif

executable(app_name, app_srcs, dependencies: sim_deps, include_directories : inc_dirs, c_args: app_args)

#copy the junction files
src = join_paths(meson.source_root(), 'junctions', get_option('junction_type'), 'junc.geom')
dst = join_paths(meson.build_root(), 'junc.geom')
run_command('rm', dst, check: true)
run_command('cp', '-s', src, dst, check: true)
src = join_paths(meson.source_root(), 'junctions', get_option('junction_type'), 'params.conf')
dst = join_paths(meson.build_root(), 'params.conf')
run_command('rm', dst, check: true)
run_command('cp', '-s', src, dst, check: true)

#only build tests if debug is set
if get_option('buildtype') == 'debug' or get_option('buildtype') == 'debugoptimized'
  executable(tst_name, tst_srcs, dependencies: sim_deps, include_directories : inc_dirs)

  #copy tests
  src = join_paths(meson.source_root(), 'tests')
  dst = join_paths(meson.build_root(), 'tests')
  message('copying @0@ to @1@ ...'.format(src,dst))
  run_command('cp', '-rs', src, dst, check: false)
endif