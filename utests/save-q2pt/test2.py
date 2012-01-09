import math
import numpy as np
import tables as tb

def make_metainfo_p(p):
    """
    result: dictionary {attr : np.array}
    """
    return {
        'index_order'   : [ 't_snk', 't_src', 'vec_snk', 'vec_src', 
                            'spin_snk', 'spin_src', 're_im' ],
        'latsize'       : np.array(p['latsize']),
        't_axis'        : p['t_axis'],
        'nvec'          : len(p['laph_pw3_list']) * len(p['laph_col_list'])
        }

def check_metainfo(p, attr_f):
    attr = make_metainfo_p(p)
    for k in attr.keys():
        if (not attr_f.has_key(k)):
            print "%s: no key\n" % k
            return False
        if (not (attr_f[k] == attr[k]).all()):
            print "%s: data mismatch: '%s' != '%s'" % (k, attr_f[k], attr[k])
            return False
    return True

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
    t0          = src0[t_axis]
    prop_mom_c0 = np.asarray(prop_mom_c0)
    prop_mom_sh = np.asarray(prop_mom_sh)
    prop_exp_sh = np.asarray(prop_exp_sh)
    list_mom    = np.asarray(list_mom)

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
    x0 = full2space(src0, t_axis)
    xp = full2space(prop_mom_c0, t_axis)
    ls = full2space(geom, t_axis)
    mom12   = pw_prod_matr(ls, (-list_mom, x0), (list_mom, x0), (prop_mom_sh, xp))
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
    return (vol3 * col12 * s12 * mom12 * exp12).reshape((lt, lt, n_v, n_v, 4, 4))


def make_test_res_p(p):
    return make_test_res(p['latsize'], p['laph_pw_x0'], p['t_axis'],
                         p['prop_pw3_sh'], p['prop_pw_x0'], p['prop_texp_sh'], 
                         p['laph_pw3_list'], p['laph_col_list'], p['laph_texp'])


def check_output_data(p, test_dir='.'):
    h5  = tb.openFile(test_dir + '/' + p['h5_file'])
    y   = h5.getNode('/' + p['h5_path'])
    attr_f  = y.attrs
    y   = y[...,0] + 1j * y[...,1]
    h5.close()

    y1  = make_test_res_p(p)
    
    cmp_data    = np.allclose(y, y1)
    cmp_meta    = check_metainfo(p, attr_f.__dict__)    

    print y.shape, y1.shape
    dy  = y - y1    
    print math.sqrt((y*y.conj()).sum()), \
          math.sqrt((y1*y1.conj()).sum()), \
          math.sqrt((dy*dy.conj()).sum())
    
    return (cmp_data and cmp_meta)

def check_qlua_run(p, test_dir='.'):
    f = open(test_dir + '/stdout', 'r')
    for l in f:
        if (l == 'QLUA_RUN_SUCCESS\n'):
            return True
    return False
    

case_list = [
    {   'case'      : 'save-q2pt-test2a',
        'qlua_src'  : 'save-q2pt/test2a.qlua',
        'make_data' : make_test_res_p,
        'make_meta' : make_metainfo_p,
        'check_log' : check_qlua_run,
        'check_out' : check_output_data },
    {   'case'      : 'save-q2pt-test2b',
        'qlua_src'  : 'save-q2pt/test2b.qlua',
        'make_data' : make_test_res_p,
        'make_meta' : make_metainfo_p,
        'check_log' : check_qlua_run,
        'check_out' : check_output_data }
]

input_list = [
   {'h5_file'           : 'test.h5',
    'h5_path'           : '/test1',
    'latsize'           : [6,6,6,12],
    't_axis'            : 3,
    'laph_pw_x0'        : [0,0,0,0],
    'prop_pw3_sh'       : [[0,0,0], [0,0,0], [0,0,0], [0,0,0]],
    'prop_pw_x0'        : [0,0,0,0],
    'prop_texp_sh'      : [.0,.0,.0,.0],
    'laph_pw3_list'     : [[0,0,0]],
    'laph_col_list'     : [[1,0,0],[0,1,0],[0,0,1]],
    'laph_texp'         : [.0,-.0, -.0] },

   {'h5_file'           : 'test.h5',
    'h5_path'           : '/test1',
    'latsize'           : [6,6,6,12],
    't_axis'            : 3,
    'laph_pw_x0'        : [1,2,3,4],
    'prop_pw3_sh'       : [[1,0,0],[0,1,0],[0,0,1],[1,1,0]],
    'prop_pw_x0'        : [3,2,1,5],
    'prop_texp_sh'      : [ .15, .16, .18-.12j, .20],
    'laph_pw3_list'     : [[0,0,0], [1,0,0], [0,1,0], [0,0,1]],
    'laph_col_list'     : [[1,0,0],[0,1,0],[0,0,1]],
    'laph_texp'         : [ .11, -.12, -.13+.08j ] },

   {'h5_file'           : 'test.h5',
    'h5_path'           : '/test1',
    'latsize'           : [6,6,6,12],
    't_axis'            : 3,
    'laph_pw_x0'        : [1,2,3,4],
    'prop_pw3_sh'       : [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ],
    'prop_pw_x0'        : [3,2,1,5],
    'prop_texp_sh'      : [-.015, .16, .18, .20],
    'laph_pw3_list'     : [ [0,0,0],
                            [1,0,0], [0,1,0], [0,0,1],
                            [1,1,0], [1,0,1], [0,1,1],
                            [1,1,1],
                            [2,0,0], [3,0,0], [4,0,0] ],
    'laph_col_list'     : [ [-0.38668154,  0.91713577,  0.59997194],
                            [ 0.20150464, -0.19782945,  0.35816456],
                            [-0.9190735 ,  1.04076121,  1.02879648] ],
    'laph_texp'         : [ .07, -.14, .21] }
    ]


