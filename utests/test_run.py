#!/usr/bin/env python2.7

import sys
import re
import numpy as np
import math

# functions for testing test field contractions
# TODO maybe it deserves a separate file/module?
def pw_prod_matr(latsize, *p_x0_list):
    """ return a normalized matrix of products of plane waves:
        M[{{i^(a)}}] = \sum_x[ exp(I * \sum_a[ p^(a)_{i^(a)} . (x - x^(a)_i) ]) ] / \sum_x[ 1 ]
            latsize     lattice size
            p_x0_list = [ ( [ p^(1), ... ], x^(a) ), ... ]
    """
    p  = [ np.array(t[0]) for t in p_x0_list ]
    xx = np.array([ t[1] for t in p_x0_list ])
    ls = np.array(latsize, float)
    res_sh = tuple(len(p_i) for p_i in p)

    def ft_val(*idx):
        pp = np.array([ p[i][idx[i]] for i in range(len(idx)) ])
        if np.allclose(0, pp.sum(axis=0) % ls):
            phase = (pp[1:] * (xx[1] - xx[1:]) / ls).sum()
            return np.exp(2j * math.pi * phase)
        else: 
            return 0.

    return np.fromfunction(np.vectorize(ft_val, otypes=[np.complex128]),
                           res_sh, dtype=int)

def space2full(c, axis, c_t=0):
    res = [ c_i for c_i in c ]
    res.insert(axis, c_t)
    return res

def full2space(c, axis):
    idx = range(len(c))
    del(idx[axis])
    return c[idx]


gamma_dgr = np.array([ 
    [[  0,  0,  0, 1j], [  0,  0, 1j,  0], [  0,-1j,  0,  0], [-1j,  0,  0,  0]],
    [[  0,  0,  0, -1], [  0,  0,  1,  0], [  0,  1,  0,  0], [ -1,  0,  0,  0]],
    [[  0,  0, 1j,  0], [  0,  0,  0,-1j], [-1j,  0,  0,  0], [  0, 1j,  0,  0]],
    [[  0,  0,  1,  0], [  0,  0,  0,  1], [  1,  0,  0,  0], [  0,  1,  0,  0]] ])

def make_Gamma16(g_list):
    assert 4 == len(g_list)
    g1, g2, g3, g4 = g_list
    from numpy import dot, identity
    return \
            [ identity(4),                      # 0 = 0000
              g1,                               # 1 = 0001
              g2,                               # 2 = 0010
              dot(g1, g2),                      # 3 = 0011
              g3,                               # 4 = 0100
              dot(g1, g3),                      # 5 = 0101
              dot(g2, g3),                      # 6 = 0110
              dot(dot(g1, g2), g3),             # 7 = 0111

              g4,                               # 8 = 1000
              dot(g1, g4),                      # 9 = 1001
              dot(g2, g4),                      #10 = 1010
              dot(dot(g1, g2), g4),             #11 = 1011
              dot(g3, g4),                      #12 = 1100
              dot(dot(g1, g3), g4),             #13 = 1101
              dot(dot(g2, g3), g4),             #14 = 1110
              dot(dot(g1, g2), dot(g3, g4)) ]   #15 = 1111

Gamma16_dgr = np.array(make_Gamma16(gamma_dgr))

def tensor_matrix_dot(t, m, axis_t):
    """ right-multiply tensor . matrix, using 'axis_t'th index of the tensor 
        preserving the meaning of all the other tensor indices
    """
    assert 2 == len(m.shape)
    return np.rollaxis(np.tensordot(t, m, axes=(axis_t, 0)), -1, axis_t)
def matrix_tensor_dot(m, t, axis_t):
    """ left-multiply tensor . matrix, using 'axis_t'th index of the tensor 
        preserving the meaning of all the other tensor indices
    """
    assert 2 == len(m.shape)
    return np.rollaxis(np.tensordot(t, m, axes=(axis_t, 1)), -1, axis_t)

def cplx_norm2(a): 
    a   = np.asarray(a)
    return (a.real*a.real).sum() + (a.imag*a.imag).sum()

# reporting
def rdiff(a, b):
    if 0 == np.abs(a).sum():
        if 0 == np.abs(b).sum(): return 0
        else: return 1.
    else:
        if 0 == np.abs(a).sum(): return 1.
        else: return float(np.abs(a - b).sum()) / (np.abs(a).sum() + np.abs(b).sum())

# qlua code generation
def qlua_string(v):
    if list == type(v):
        if (len(v) <= 0): return '{}'
        s = '{ %s' % qlua_string(v[0])
        for x in v[1:]:
            s += ', %s' % qlua_string(x)
        return s + ' }'
    elif int == type(v) or float == type(v): return str(v)
    elif complex == type(v): return "complex(%s,%s)" % (str(v.real), str(v.imag))
    elif str == type(v): return "'%s'" % v
    else:
        raise ValueError, "cannot convert to qlua string: %s: %s" % (str(type(v)), str(v))

def make_qlua_test_param(p, qlua_test_param, test_dir, qlua_qlib=None):
    f   = open(qlua_test_param, 'w')
    if None != qlua_qlib:
        f.write("package.path = package.path .. ';%s/?.qlua'\n" % qlua_qlib)
    f.write('test_dir = %s\n' % qlua_string(test_dir))
    for k,v in p.iteritems():
        f.write('%s = %s\n' % (k, qlua_string(v)))
    f.close()

def run_qlua_sh(qlua_bin, qlua_test_param, qlua_test_src, 
            qlua_stdout, qlua_stderr,
            **run_param):
    import os
    cmd = ('time %s %s %s >%s 2>%s </dev/null' % (
                        qlua_bin, qlua_test_param, qlua_test_src,
                        qlua_stdout, qlua_stderr))
    print cmd
    res = os.system(cmd)

    return (0 == res)

def run_qlua_mpi(qlua_bin, qlua_test_param, qlua_test_src, 
            qlua_stdout, qlua_stderr,
            **run_param):
    import os
    if (run_param.has_key('mpi')): mpiopts = run_param['mpi']
    else: mpiopts = ''
    if (run_param.has_key('qmp')): qmpopts = run_param['qmp']
    else: qmpopts = ''
    if (run_param.has_key('mpirun')): mpirun = run_param['mpirun']
    else: mpirun = 'mpirun'
    cmd = '%s %s %s %s %s %s >%s 2>%s </dev/null ' % (
                    mpirun,
                    mpiopts,
                    qlua_bin, 
                    qmpopts,
                    qlua_test_param,
                    qlua_test_src,
                    qlua_stdout, 
                    qlua_stderr)
    print cmd
    res = os.system(cmd)
    return (0 == res)

def run_qlua_cobalt(qlua_bin, qlua_test_param, qlua_test_src, 
            qlua_stdout, qlua_stderr,
            **run_param):
    # TODO implement actual call
    raise NotImplementedError
    return False 

def timer(): 
    time0 = time.time() 
    def dt():  
        return time.time() - time0 
    return dt 

def parse_comma_list(s):
    """ parse a list of A=B,C=D,E,F... parameters """
    l = []
    d = {}
    for p in s.split(','):
        rp = re.match('^([^=]+)=(.*)$', p)
        if None == rp: l.append(p)
        else: d[rp.groups()[0]] = rp.groups()[1]
    return l, d



# main loop
def call_main(what, run_name, case, input, run_param_str=None):
    """
            what            run|ch_log|ch_out
            run_name        run name
            case            test case
            input           test input
            run_param_str   run parameters
    """
    import os

    def get_env(env_name, default=''):
        if os.environ.has_key(env_name): return os.environ[env_name]
        else: return default

    qlua_top        = get_env('QLUA_TOP', '.')
    qlua_bin        = get_env('QLUA_BIN', qlua_top + '/qlua')
    qlua_utests_dir = get_env('QLUA_UTESTS_DIR', qlua_top + '/utests') 
    qlua_test_src   = get_env('QLUA_TEST_SRC', qlua_utests_dir + '/' + case['qlua_src'])
    qlua_qlib       = get_env('QLUA_QLIB', qlua_top + '/qlib')
    scratch_dir     = get_env('SCRATCH_DIR', '.') + '/qlua_tests'

    def debug_print_param():
        print "test_name\t=",       test_name
        print "test_src\t=",        test_src
        print "qlua_top\t=",        qlua_top
        print "qlua_bin\t=",        qlua_bin
        print "qlua_utests_dir\t=", qlua_utests_dir
        print "qlua_test_src\t=",   qlua_test_src
        print "qlua_qlib\t=",       qlua_qlib
        print "scratch_dir\t=",     scratch_dir

    #debug_print_param()


    test_dir        = scratch_dir + '/' + run_name
    qlua_test_param = test_dir + '/param.qlua'
    qlua_stdout     = test_dir + '/stdout'
    qlua_stderr     = test_dir + '/stderr'
    
    def print_res(run_name, res):
        if res: print run_name, "\tOK"
        else:   print run_name, "\tFail"

    # TODO parse statements like "run,ch_out" and do the actions (in order?)
    if ('run' == what):
        if (None == run_param_str):
            raise ValueError, "run: <params> are not set"
        run_l, run_d = parse_comma_list(run_param_str)

        if (os.system('mkdir -p %s' % test_dir)):
            raise OSError, "cannot create dir '%s'" % test_dir

        make_qlua_test_param(input, qlua_test_param, test_dir, qlua_qlib)

        if (not run_d.has_key('run') or 'sh' == run_d['run']):
            res = run_qlua_sh(qlua_bin, qlua_test_param, qlua_test_src, 
                qlua_stdout, qlua_stderr, **run_d)
        elif ('mpi' == run_d['run']):
            res = run_qlua_mpi(qlua_bin, qlua_test_param, qlua_test_src, 
                qlua_stdout, qlua_stderr, **run_d)
        elif ('cobalt' == run_d['run']):
            res = run_qlua_cobalt(qlua_bin, qlua_test_param, qlua_test_src, 
                qlua_stdout, qlua_stderr, **run_d)
        else:
            raise ValueError, "unknown run mode='%s'" % run_d['run']

        print_res('[run] ' + run_name, res)

    elif ('ch_log' == what): 
        res = case['check_log'](input, test_dir)
        print_res('[log] ' + run_name, res)

    elif ('ch_out' == what): 
        res = case['check_out'](input, test_dir)
        print_res('[out] ' + run_name, res)

    else:
        raise ValueError, "unknown function '%s'" % what


if (__name__ == '__main__'):
    import sys
    if len(sys.argv) < 5:
        raise ValueError, """usage:
$ test_run.py  <test>.py  run|ch_log|ch_out
        <case>,...
        <input#>,...
        <run_param>,...
"""
        sys.exit(1)
    test_name, what_list_str, case_list_str, in_list_str = sys.argv[1:5]
    if 6 <= len(sys.argv): run_param_str = sys.argv[5]
    else: run_param_str = None

    # FIXME specify locals/globals to keep test variables separate from the current scope
    execfile(test_name)
    
    if 'all' == what_list_str: what_list = ['run', 'ch_log', 'ch_out' ]
    else: what_list = parse_comma_list(what_list_str)[0]

    # TODO find case ## by names
    if 'all' == case_list_str: i_case_list = range(len(case_list))  # replace with test's 'local' var
    else: i_case_list = map(int, parse_comma_list(case_list_str)[0])

    if 'all' == in_list_str: i_in_list = range(len(input_list))     # replace with test's 'local' var
    else: i_in_list = map(int, parse_comma_list(in_list_str)[0])

    for i_case in i_case_list:
        for i_in in i_in_list:
            for what in what_list:
                call_main(what, 
                          ('%s_%d' %(case_list[i_case]['case'], i_in)),
                          case_list[i_case], input_list[i_in], 
                          run_param_str)
