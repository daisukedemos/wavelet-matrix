from waflib import Options
import sys

top = '.'
out = 'build'
subdirs = ['lib']

def options(opt):
    opt.load('compiler_cxx unittest_gtest')
    opt.recurse(subdirs)

def configure(conf):
    conf.env.CXXFLAGS += ['-W', '-Wall', '-Wextra', '-Wold-style-cast', '-g', '-O2', '-fopenmp']
#   conf.env.CXXFLAGS += ['-W', '-Wall', '-Wextra', '-Wold-style-cast', '-g', '-fopenmp']
    conf.load('compiler_cxx unittest_gtest')
    conf.check_cfg(package= 'rsdic', args='--cflags --libs', uselib_store='RSDIC')
    conf.recurse(subdirs)


def build(bld):
    bld.recurse(subdirs)
