import math
import numpy as np
import tables as tb
import sys
import os

"""
h5file  = 'squark-wf.test1.hdf5'
h5path  = 'test1'

h5f     = tb.openFile(h5file, mode='r')
sqwf    = h5f[h5path][:]
"""



def make_test_res(geom, src0, t_axis, 
               list_psrc, list_mom, list_col, list_exp, c0):
    """ create table for test1 
        geom        lattice geometry
        src0        starting point for color vector plane waves
        t_axis      time axis index
        list_psrc   list of src projection momenta
        list_mom    list of ev momenta
        list_col    list of color vectors
        list_exp    t-exponents for color vectors
        c0          starting point for src Fourier transform
    """
    geom    = np.asarray(geom)
    src0    = np.asarray(src0)
    c0      = np.asarray(c0)
    lt      = geom[t_axis]
    vol3    = geom.prod() / lt
    t0      = src0[t_axis]
    n_p     = len(list_psrc)
    n_c     = len(list_col)
    n_mom   = len(list_mom)
    n_v     = n_mom * n_c

    # color part
    list_col    = np.asarray(list_col)
    col012  = np.array([ [ [ np.linalg.det([c1, c2, c3])
                             for c3 in list_col ]
                           for c2 in list_col ]
                         for c1 in list_col ], 
                       dtype=np.complex128)
    # spatial mom part
    list_psrc   = np.asarray(list_psrc)
    list_mom    = np.asarray(list_mom)
    def equal_sum_mom(l1, l2, l3, lsum):
        l1  = np.asarray(l1)
        l2  = np.asarray(l2)
        l3  = np.asarray(l3)
        def func(i,j,k, i_p):
            return np.where(l1[i] + l2[j] + l3[k] == lsum[i_p], 1, 0)
        return func

    mom012  = np.fromfunction(
                (lambda i,j,k, i_p: 
                    np.where(np.all(
                        list_mom[i] + list_mom[j] + list_mom[k] == list_psrc[i_p], axis=-1),
                        1, 0)),
                (n_mom, n_mom, n_mom, n_p), 
                dtype=int)
    s_axis  = range(len(geom)); del s_axis[t_axis]
    mom012  = mom012 * np.exp(2j * math.pi 
                        * np.dot(list_psrc, ((src0 - c0)/np.array(geom, float))[s_axis]))
    #print np.dot(list_psrc, ((src0 - c0)/np.array(geom, float))[s_axis])

    # time dep
    list_exp    = np.asarray(list_exp)
    exp012      = np.exp(list_exp.sum() * ((np.arange(lt) + lt - t0) % lt))

    #print col012[None, :, None, :, None, :, None, None].shape
    #print mom012[:, None, :, None, :, None, None, :].shape
    #print exp012[..., :, None].shape
    return (vol3* col012[None, :, None, :, None, :, None, None] \
                * mom012[:, None, :, None, :, None, None, :] \
                * exp012[..., :, None]).reshape((n_v, n_v, n_v, lt, n_p))





test_list = [
  { 'name'      : 'squark-wf.test1-1.hdf5',
    'path'      : 'test1',
    'geom'      : [ 6, 6, 6, 12 ],
    't_axis'    : 3,
    'src0'      : [ 0, 0, 0, 0 ],
    'c0'        : [ 0, 0, 0, 0 ],
    'list_mom'  : [ [0,0,0], [1,0,0] ],
    'list_psrc' : [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ],
    'list_col'  : [ [1,0,0], [0,1,0], [0,0,1] ],
    'list_exp'  : [ 0., 0., 0.] },
  { 'name'      : 'squark-wf.test1-2.hdf5',
    'path'      : 'test1',
    'geom'      : [ 6, 6, 6, 12 ],
    't_axis'    : 3,
    'src0'      : [ 0, 0, 0, 0 ],
    'c0'        : [ 0, 0, 0, 0 ],
    'list_mom'  : [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ],
    'list_psrc' : [ [0,0,0], 
                    [1,0,0], [0,1,0], [0,0,1],
                    [1,1,0], [1,0,1], [0,1,1],
                    [1,1,1], 
                    [2,0,0], [3,0,0], [4,0,0] ],
    'list_col'  : [ [1,0,0], [0,1,0], [0,0,1] ],
    'list_exp'  : [ 0.1, 0.2, -0.3] },
  { 'name'      : 'squark-wf.test1-3.hdf5',
    'path'      : 'test1',
    'geom'      : [ 6, 6, 6, 12 ],
    't_axis'    : 3,
    'src0'      : [ 0, 0, 0, 0 ],
    'c0'        : [ 0, 0, 0, 0 ],
    'list_mom'  : [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ],
    'list_psrc' : [ [0,0,0], 
                    [1,0,0], [0,1,0], [0,0,1],
                    [1,1,0], [1,0,1], [0,1,1],
                    [1,1,1], 
                    [2,0,0], [3,0,0], [4,0,0] ],
    'list_col'  : [ [-0.38668154,  0.91713577,  0.59997194],
                    [ 0.20150464, -0.19782945,  0.35816456],
                    [-0.9190735 ,  1.04076121,  1.02879648] ],
    'list_exp'  : [ 0.1, -0.2, 0.4 ] }
    ]


qlua_bin    = './qlua'
qlua_test   = 'utests/squark-wf/test1a.qlua'

def qlua_string(v):
    if list == type(v):
        if (len(v) <= 0): return '{}'
        s = '{ %s' % qlua_string(v[0])
        for x in v[1:]:
            s += ', %s' % qlua_string(x)
        return s + ' }'
    elif int == type(v) or float == type(v): return str(v)
    elif str == type(v): return "'%s'" % v
    else:
        raise ValueError, "cannot convert to qlua string: %s: %s" % (str(type(v)), str(v))

def call_qlua(p):
    qlua_case = p['name'] + '.' + p['path'] + '.qlua'
    f   = open(qlua_case, 'w')
    for k,v in p.iteritems():
        f.write('%s = %s\n' % (k, qlua_string(v)))
    f.close()

    os.system('rm -f %s' % p['name'])
    return os.system(qlua_bin + ' ' + qlua_case + ' ' + qlua_test)



def check_test(p):
    h5  = tb.openFile(p['name'])
    y   = h5.getNode('/' + p['path'])
    y   = y[...,0] + 1j * y[...,1]
    y1  = make_test_res(p['geom'], p['src0'], p['t_axis'],
                        p['list_psrc'], p['list_mom'], p['list_col'], p['list_exp'], 
                        p['c0'])
    res = np.allclose(y, y1)
    
    h5.close()

    return res
    

if (__name__ == '__main__'):
    for i,p in enumerate(test_list):
        stat = call_qlua(p)
        if 0 != stat:
            print "%s exited with stat=%d\n" % (qlua_bin, stat)
            continue

        res = check_test(p)
        print i, res
