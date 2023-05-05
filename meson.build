project(
    'ntCard',
    'cpp',
    version: '1.3.0',
    license: 'MIT',
    default_options: [
        'cpp_std=c++17',
        'werror=true',
        'warning_level=3',
        'optimization=3',
    ],
)

conf = configuration_data()
conf.set('VERSION', meson.project_version())
configure_file(
    input: 'include/ntcard/version.hpp',
    output: 'version.hpp',
    configuration: conf,
)

compiler = meson.get_compiler('cpp')

dependencies = [
    dependency('openmp'),
    compiler.find_library('btllib', static: true),
    compiler.find_library('nthash', static: true),
]

include_dirs = [
    include_directories('include'),
    include_directories('vendor/argparse/include'),
]

executable(
    'ntcard',
    'src/main.cpp',
    dependencies: dependencies,
    include_directories: include_dirs,
    install: true,
)
