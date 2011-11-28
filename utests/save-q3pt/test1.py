import math
import numpy as np
import tables as tb


def make_metainfo_p(p):
    """
    result: dictionary {attr : np.array}
    """
    # TODO
    raise NotImplementedError
    return None

def make_data(latsize, t_axis, 
            laph_pw3_list, laph_pw_x0, laph_col_list, laph_texp,
            prop_pw3_sh, prop_pw_x0, prop_texp_sh, 
            ft_p3_list, ft_x0, 
            tstep, src_snk_dt_max):
    """ create table for test1 
        latsize             lattice geometry
        t_axis              time axis index
        laph_pw3_list       list of ev momenta
        laph_pw_x0          starting point for color vector plane waves
        laph_col_list       list of color vectors
        laph_texp[col]      t-exponents for color vectors
        prop_pw3_sh[spin]   prop pw shift
        prop_pw_x0          prop pw start point
        prop_texp_sh[spin]  prop time-dependent exp
        ft_p3_list          FT momenta list
        ft_x0               FT starting point
        tstep               time step
        src_snk_dt_max      maximal src-snk dist
    result: complex [i_t12op, nsnk, nsrc, spin_snk, spin_src, n_op, n_qmom]
    """
    N = None

    latsize         = np.asarray(latsize)
    lt              = latsize[t_axis]
    vol3            = latsize.prod() / lt
    laph_pw3_list   = np.asarray(laph_pw3_list)
    laph_pw_x0      = np.asarray(laph_pw_x0)
    laph_col_list   = np.asarray(laph_col_list)
    laph_texp       = np.asarray(laph_texp)
    prop_pw_x0      = np.asarray(prop_pw_x0)
    prop_pw3_sh     = np.asarray(prop_pw3_sh)
    prop_texp_sh    = np.asarray(prop_texp_sh)
    ft_p3_list      = np.asarray(ft_p3_list)
    ft_x0           = np.asarray(ft_x0)

    # list_t12op[i_t12op,3]
    list_t12op = []
    for tsrc in range(0, lt, tstep):
        for dt in range(2*tstep, src_snk_dt_max + 1, tstep):
            tsnk = (lt + tsrc + dt) % lt
            for dtau in range(tstep, dt, tstep):
                t_op = (lt + tsrc + dtau) % lt
                list_t12op.append((tsrc, tsnk, t_op))
    list_t12op  = np.asarray(list_t12op)
    n_t12op     = len(list_t12op)

    def dist(a, b, lt):    
        dpos = (lt + a - b) % lt
        dneg = (lt + b - a) % lt
        return np.where(dpos < dneg, dpos, dneg)

    # [i_t12op]
    dist_tsrc   = dist(list_t12op[:,2], list_t12op[:,0], lt)
    dist_tsnk   = dist(list_t12op[:,2], list_t12op[:,1], lt)
    lin_laph    = (lt + list_t12op[:,2] - laph_pw_x0[t_axis]) % lt

    # color part [i_t12op, c1, c2]: laph_texp, laph_col
    col12   = np.array([ [ np.dot(c1.conj(), c2) for c2 in laph_col_list ]
                         for c1 in laph_col_list ], 
                       dtype=np.complex128)[N,:,:] \
              * np.exp(lin_laph[:,N,N] * (laph_texp[N,:,N].conj() + 
                                          laph_texp[N,N,:]))
    # reshape [i_t12op, p1, c1, p2, c2, s1, s2, i_op, i_qmom]
    col12   = col12[:, N, :, N, :, N, N, N, N]

    # spin part [i_t12op, s1, s2, iGamma] : prop_texp_sh
    spin12  = Gamma16_dgr.transpose(1,2,0)[N,:,:,:] \
              * np.exp(dist_tsnk[:,N,N,N] * prop_texp_sh[N,:,N,N].conj() +
                       dist_tsrc[:,N,N,N] * prop_texp_sh[N,N,:,N])
    # reshape [i_t12op, p1, c1, p2, c2, s1, s2, i_op, i_qmom]
    spin12  = spin12[:, N, N, N, N, :, :, :, N]

    # spatial mom part [p1, p2, s1, s2, i_qmom] : laph_pw3, prop_pw3, ft_pw3
    x0 = full2space(laph_pw_x0, t_axis)
    xp = full2space(prop_pw_x0, t_axis)
    xq = full2space(ft_x0, t_axis)
    ls = full2space(latsize, t_axis)
    mom12   = pw_prod_matr(ls, (-laph_pw3_list, x0), (laph_pw3_list, x0), 
                               (-prop_pw3_sh, xp), (prop_pw3_sh, xp),
                               (ft_p3_list, xq))
    # reshape [i_t12op, p1, c1, p2, c2, s1, s2, i_op, i_qmom]
    mom12   = mom12[N, :, N, :, N, :, :, N, :]

    # complex [i_t12op, n1, n2, s1, s2, i_op, i_qmom]
    n_v     = len(laph_col_list) * len(laph_pw3_list)
    n_qmom  = len(ft_p3_list)
    res = (col12 * spin12 * mom12).reshape((n_t12op, n_v, n_v, 4, 4, 16, n_qmom))
    return res

def make_data_p(p):
    return make_data(p['latsize'], p['t_axis'], 
            p['laph_pw3_list'], p['laph_pw_x0'], p['laph_col_list'], p['laph_texp'],
            p['prop_pw3_sh'], p['prop_pw_x0'], p['prop_texp_sh'], 
            p['ft_p3_list'], p['ft_x0'], 
            p['tstep'], p['src_snk_dt_max'])

def check_output_data(p, test_dir='.'):
    h5  = tb.openFile(test_dir + '/' + p['h5_file'])
    y   = h5.getNode('/' + p['h5_path'])
    y   = y[...,0] + 1j * y[...,1]
    h5.close()

    y1  = make_data_p(p)

    res = np.allclose(y, y1)
    dy  = y - y1
    print y.shape, y1.shape
    print math.sqrt((y*y.conj()).sum()), math.sqrt((y1*y1.conj()).sum()), math.sqrt((dy*dy.conj()).sum())

    # TODO check metainfo

    return res

def check_qlua_run(p, test_dir='.'):
    f = open(test_dir + '/stdout', 'r')
    for l in f:
        if (l == 'QLUA_RUN_SUCCESS\n'):
            return True
    return False

case_list = [
    {   'case'      : 'save-q3pt-test1a',
        'qlua_src'  : 'save-q3pt/test1a.qlua',
        'make_data' : make_data_p,
        'make_meta' : make_metainfo_p,
        'check_log' : check_qlua_run,
        'check_out' : check_output_data }
]

input_list = [
   {'h5_file'       : 'test.h5',
    'h5_path'       : '/test1',
    'latsize'       : [6,6,6,12],
    't_axis'        : 3,
    'laph_pw3_list' : [[0,0,0]],
    'laph_pw_x0'    : [0,0,0,0],
    'laph_col_list' : [[1,0,0],[0,1,0],[0,0,1]],
    'laph_texp'     : [.0,-.0, -.0],
    'prop_pw3_sh'   : [[0,0,0], [0,0,0], [0,0,0], [0,0,0]],
    'prop_pw_x0'    : [0,0,0,0],
    'prop_texp_sh'  : [.0,.0,.0,.0],
    'ft_p3_list'    : [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ],
    'ft_x0'         : [ 0, 0, 0, 0 ],
    'tstep'         : 2,
    'src_snk_dt_max': 6 },
   {'h5_file'       : 'test.h5',
    'h5_path'       : '/test1',
    'latsize'       : [6,6,6,12],
    't_axis'        : 3,
    'laph_pw3_list' : [[1,1,1]],
    'laph_pw_x0'    : [1,2,3,4],
    'laph_col_list' : [[1,0,0],[0,1,0],[0,0,1]],
    'laph_texp'     : [ 0.1+0.12j, -0.2+0.17j, 0.3-0.19j],
    'prop_pw3_sh'   : [[1,0,0], [0,1,0], [0,0,1], [-1,0,0]],
    'prop_pw_x0'    : [4,5,0,8],
    'prop_texp_sh'  : [ -0.11+0.13j, 0.12+0.14j, -0.13-0.08j, 0.14-0.06j],
    'ft_p3_list'    : [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ],
    'ft_x0'         : [3,4,5,0],
    'tstep'         : 2,
    'src_snk_dt_max': 6 },
   {'h5_file'       : 'test.h5',
    'h5_path'       : '/test1',
    'latsize'       : [6,6,6,12],
    't_axis'        : 3,
    'laph_pw3_list' : [[0,0,0], [1,0,0], [0,1,0], [0,0,1]],
    'laph_pw_x0'    : [1,2,3,4],
    'laph_col_list' : [ [-0.38668154,  0.91713577,  0.59997194],
                        [ 0.20150464, -0.19782945,  0.35816456],
                        [-0.9190735 ,  1.04076121,  1.02879648] ],
    'laph_texp'     : [ 0.1+0.12j, -0.2+0.23j, 0.3-0.15j ],
    'prop_pw3_sh'   : [[0,0,0], [1,0,0], [0,1,0], [0,0,1]],
    'prop_pw_x0'    : [4,5,0,8],
    'prop_texp_sh'  : [ -0.21-0.15j, 0.22-0.17j, -0.23+0.18j, 0.24+0.16j],
    'ft_p3_list'    : [ [0,0,0], 
                        [1,0,0], [0,1,0], [0,0,1],
                        [1,1,0], [1,0,1], [0,1,1],
                        [1,1,1], 
                        #[2,0,0], [3,0,0] ],
                        [2,0,0], [3,0,0], [4,0,0] ],
    #'ft_p3_list'    : [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ],
    'ft_x0'         : [3,1,4,11],
    'tstep'         : 2,
    'src_snk_dt_max': 6 }
]
