import math
import numpy as np
import tables as tb


def make_test_res(geom, src0, t_axis, list_psrc, list_mom, list_col, list_exp, c0):
    """ create table for test1: v123_ft = det(v^1,v^2,v^3)_x exp(-i*psrc*x) 
        with v^a_i = exp(list_exp[a]*t) * list_col[j] * exp(i*list_mom[k]*x), i=index(j,k)
        geom                lattice geometry
        src0                starting point for color vector plane waves
        t_axis              time axis index
        list_psrc[i_psrc]   list of src projection momenta
        list_mom[i_mom]     list of ev momenta
        list_col[i_col]     list of color vectors
        list_exp[3]         t-exponents for color vectors
        c0                  starting point for src Fourier transform
        result[n_v, n_v, n_v, lt, n_p], n_v = n_mom * n_col
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

def make_test_res_p(p):
    return make_test_res(
        p['latsize'], p['laph_pw_x0'], p['t_axis'], p['ft_p3_list'], 
        p['laph_pw3_list'], p['laph_col_list'], p['laph_texp'], p['ft_x0'])


        
def check_output_data(p, test_dir='.'):
    h5  = tb.openFile(test_dir + '/' + p['h5_file'])
    y   = h5.getNode('/' + p['h5_path'])
    y   = y[...,0] + 1j * y[...,1]
    h5.close()
    
    y1  = make_test_res_p(p)
    
    res = np.allclose(y, y1)

    # TODO check metainfo

    return res
    
def check_qlua_run(p, test_dir='.'):
    f = open(test_dir + '/stdout', 'r')
    for l in f:
        if (l == 'QLUA_RUN_SUCCESS\n'):
            return True
    return False

case_list = [
    {   'case'      : 'squark-wf-test2a',
        'qlua_src'  : 'squark-wf/test2a.qlua',
        'make_data' : make_test_res_p,
        'check_log' : check_qlua_run,
        'check_out' : check_output_data }
]

input_list = [
  { 'h5_file'       : 'test0.hdf5',
    'h5_path'       : 'test1',
    'latsize'       : [ 6, 6, 6, 12 ],
    't_axis'        : 3,
    'laph_pw_x0'    : [ 0, 0, 0, 0 ],
    'ft_x0'         : [ 0, 0, 0, 0 ],
    'laph_pw3_list' : [ [0,0,0], [1,0,0] ],
    'ft_p3_list'    : [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ],
    'laph_col_list' : [ [1,0,0], [0,1,0], [0,0,1] ],
    'laph_texp'     : [ 0., 0., 0.] },
  { 'h5_file'       : 'test1.hdf5',
    'h5_path'       : 'test1',
    'latsize'       : [ 6, 6, 6, 12 ],
    't_axis'        : 3,
    'laph_pw_x0'    : [ 0, 0, 0, 0 ],
    'ft_x0'         : [ 0, 0, 0, 0 ],
    'laph_pw3_list' : [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ],
    'ft_p3_list'    : [ [0,0,0], 
                        [1,0,0], [0,1,0], [0,0,1],
                        [1,1,0], [1,0,1], [0,1,1],
                        [1,1,1], 
                        [2,0,0], [3,0,0], [4,0,0] ],
    'laph_col_list' : [ [1,0,0], [0,1,0], [0,0,1] ],
    'laph_texp'     : [ 0.1, 0.2, -0.3] },
  { 'h5_file'       : 'test2.hdf5',
    'h5_path'       : 'test1',
    'latsize'       : [ 6, 6, 6, 12 ],
    't_axis'        : 3,
    'laph_pw_x0'    : [ 0, 0, 0, 0 ],
    'ft_x0'         : [ 0, 0, 0, 0 ],
    'laph_pw3_list' : [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ],
    'ft_p3_list'    : [ [0,0,0], 
                        [1,0,0], [0,1,0], [0,0,1],
                        [1,1,0], [1,0,1], [0,1,1],
                        [1,1,1], 
                        [2,0,0], [3,0,0], [4,0,0] ],
    'laph_col_list' : [ [-0.38668154,  0.91713577,  0.59997194],
                        [ 0.20150464, -0.19782945,  0.35816456],
                        [-0.9190735 ,  1.04076121,  1.02879648] ],
    'laph_texp'     : [ 0.1, -0.2, 0.4 ] }
    ]

