import math
import numpy as np
import tables as tb
import sys
import os

func_name   = 'save-q2pt'
test_name   = 'test1'

list_test = [
   {'name'            : func_name + '.' + test_name + '-0',
    'path'            : '/test1',
    'geom'            : [6,6,6,12],
    't_axis'          : 3,
    'src0'            : [0,0,0,0],
    'prop_mom_sh'     : [[0,0,0], [0,0,0], [0,0,0], [0,0,0]],
    'prop_mom_c0'     : [0,0,0,0],
    'prop_exp_sh'     : [.0,.0,.0,.0],
    'list_mom'        : [[0,0,0]],
    'list_col'        : [[1,0,0],[0,1,0],[0,0,1]],
    'list_exp'        : [.0,-.0, -.0] },

   {'name'            : func_name + '.' + test_name + '-1',
    'path'            : '/test1',
    'geom'            : [6,6,6,12],
    't_axis'          : 3,
    'src0'            : [1,2,3,4],
    'prop_mom_sh'     : [[1,0,0],[0,1,0],[0,0,1],[1,1,0]],
    'prop_mom_c0'     : [3,2,1,5],
    'prop_exp_sh'     : [.15,.16,.18,.20],
    'list_mom'        : [[0,0,0], [1,0,0], [0,1,0], [0,0,1]],
    'list_col'        : [[1,0,0],[0,1,0],[0,0,1]],
    'list_exp'        : [.11,-.12, -.13] }
    #'list_exp'        : [.0,-.0, -.0] }
    ]

def make_test_res(geom, src0, t_axis,
                  prop_mom_sh, prop_mom_c0, prop_exp_sh, 
                  list_mom, list_col, list_exp):
    """ create table for test1 
        geom        lattice geometry
        src0        starting point for color vector plane waves
        t_axis      time axis index
        prop_mom_sh prop pw shift
        prop_mom_c0 prop pw start point
        prop_exp_sh prop time-dependent exp
        list_mom    list of ev momenta
        list_col[col]   list of color vectors
        list_exp[col]   t-exponents for color vectors
    """

    geom        = np.asarray(geom)
    src0        = np.asarray(src0)
    lt          = geom[t_axis]
    vol3        = geom.prod() / lt
    t0      = src0[t_axis]
    prop_mom_c0 = np.asarray(prop_mom_c0)
    prop_mom_sh = np.asarray(prop_mom_sh)
    prop_exp_sh = np.asarray(prop_exp_sh)

    n_c     = len(list_col)
    n_mom   = len(list_mom)
    n_v     = n_mom * n_c


    # color part [c1, c2]
    list_col    = np.asarray(list_col)
    col12   = np.array([ [ np.dot(c1.conj(), c2) 
                           for c2 in list_col ]
                         for c1 in list_col ], 
                       dtype=np.complex128) 
    # complex [t1, t2, p1, c1, p2, c2, s1, s2]
    col12   = col12[None, None, None, :, None, :, None, None]

    # spin part [s1, s2]
    s12     = np.identity(4)
    # complex [t1, t2, p1, c1, p2, c2, s1, s2]
    s12     = s12[None, None, None, None, None, None, :, :] 

    # spatial mom part [p1, p2, s2] 
    # = delta(list_mom[p1]==list_mom[p2]+prop_mom_sh[s2])
    # ---  *exp(2pi*I*list_mom[p1].(src0-prop_mom_c0))
    #   *exp(2pi*I*prop_mom_sh[p1].(src0-prop_mom_c0))
    prop_mom_sh = np.asarray(prop_mom_sh)
    list_mom    = np.asarray(list_mom)
    mom12   = np.fromfunction(
                (lambda p1,p2, s2:
                    np.where(np.all(
                        list_mom[p1] == list_mom[p2] + prop_mom_sh[s2], axis=-1),
                        1, 0)),
                (n_mom, n_mom, 4), 
                dtype=int)
    s_axis  = range(len(geom)); del s_axis[t_axis]
    #mom12   = mom12 * np.exp(2j * math.pi * np.dot(list_mom, 
    #                   ((src0 - prop_mom_c0) / np.array(geom, float))[s_axis]))[:,None,None]
    mom12   = mom12 * np.exp(2j * math.pi * np.dot(prop_mom_sh, 
                       ((src0 - prop_mom_c0) / np.array(geom, float))[s_axis]))
    # complex [t1, t2, p1, c1, p2, c2, s1, s2]
    mom12   = mom12[None, None, :, None, :, None, None, :]

    # time exp part [t1, t2, c1, c2, s2]
    t1      = np.arange(lt)[:, None, None, None, None]
    le1     = np.asarray(list_exp)[None, None, :, None, None]
    t2      = np.arange(lt)[None, :, None, None, None]
    le2     = np.asarray(list_exp).conj()[None, None, None, :, None]
    dt      = (lt + t1 - t2) % lt
    le_p    = np.asarray(prop_exp_sh)[None, None, None, None, :]
    exp12   = np.exp(t1 * le1 + t1 * le2 + dt * le_p)  # sic!
    # complex [t1, t2, p1, c1, p2, c2, s1, s2]
    exp12   = exp12[:, :, None, :, None, :, None, :]

    # complex [t1, t2, p1, c1, p2, c2, s1, s2]
    #print col12.shape
    #print s12.shape
    #print mom12.shape
    #print exp12.shape
    # complex [t1, t2, n1, n2, s1, s2]
    return (vol3 * col12 * s12 * mom12 * exp12).reshape((lt, lt, n_v, n_v, 4, 4))

def make_test_res_p(p):
    return make_test_res(p['geom'], p['src0'], p['t_axis'],
                         p['prop_mom_sh'], p['prop_mom_c0'], p['prop_exp_sh'], 
                         p['list_mom'], p['list_col'], p['list_exp'])

def qlua_string(v):
    if list == type(v):
        if (len(v) <= 0): return '{}'
        s = '{ %s' % qlua_string(v[0])
        for x in v[1:]:
            s += ', %s' % qlua_string(x)
        return s + ' }'
    elif int == type(v) or float == type(v): return str(v)
    elif str == type(v): return "'%s'" % v
    elif complex == type(v): return "complex(%s,%s)" % (str(v.real), str(v.imag))
    else:
        raise ValueError, "cannot convert to qlua string: %s: %s" % (str(type(v)), str(v))


qlua_bin    = './qlua'
qlua_test   = 'utests/save-q2pt/test1a.qlua'

def call_qlua(p, path_prefix=''):
    qlua_case   = path_prefix + p['name'] + '.qlua'
    hdf5_out    = path_prefix + p['name'] + '.h5'
    f   = open(qlua_case, 'w')
    for k,v in p.iteritems():
        f.write('%s = %s\n' % (k, qlua_string(v)))
    f.write('hdf5_out = %s\n' % qlua_string(hdf5_out))
    f.close()

    os.system('rm -f %s' % hdf5_out)
    return os.system(qlua_bin + ' ' + qlua_case + ' ' + qlua_test)


def check_test(p, path_prefix=''):
    hdf5_out    = path_prefix + p['name'] + '.h5'
    h5  = tb.openFile(hdf5_out)
    y   = h5.getNode('/' + p['path'])
    y   = y[...,0] + 1j * y[...,1]
    y1  = make_test_res(p['geom'], p['src0'], p['t_axis'],
                        p['prop_mom_sh'], p['prop_mom_c0'], p['prop_exp_sh'], 
                        p['list_mom'], p['list_col'], p['list_exp'])
    res = np.allclose(y, y1)
    
    h5.close()

    return res
    

if (__name__ == '__main__'):
    for i,p in enumerate(test_list):
        qlua_stat = call_qlua(p)
        if 0 != qlua_stat:
            qlua_stat_str   = 'Fail(%d)' % qlua_stat
            check_stat      = False
            check_stat_str  = 'None'
            continue
        else:
            qlua_stat_str   = 'OK'
            check_stat = check_test(p)
            if check_stat:  check_stat_str = 'OK'
            else: check_stat_str = 'Fail'
        print "%d\t%s\t%s\t%s\n" % (i, p['name'], qlua_stat_str, check_stat_str)
