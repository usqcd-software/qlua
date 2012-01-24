import math
import numpy as np
import tables as tb
import aff

def make_dirprop_mat(c0, c_ac, c_as, c_bc, c_bs):
    r_c = np.arange(3)
    r_s = np.arange(4)
    N = None
    return (c0 + c_ac *r_c[:,N,N,N] + c_as *r_s[N,:,N,N] 
            + c_bc *r_c[N,N,:,N] + c_bs *r_s[N,N,N,:])

def make_test_res(geom, 
    frw_lin_coeff, frw_ft_mom, 
    bkw_lin_coeff, bkw_ft_mom, ft_x0):
    
    frw_mat = make_dirprop_mat(*frw_lin_coeff)
    bkw_mat = make_dirprop_mat(*bkw_lin_coeff)
    bkw_mat_h = matrix_tensor_dot(
                    Gamma16_dgr[-1],
                    tensor_matrix_dot(
                        bkw_mat.conj().transpose(2,3,0,1), # herm.conj
                        Gamma16_dgr[-1], 3),
                    1)
    print frw_mat.shape, bkw_mat_h.shape
    res = np.empty((16, 3,4,3,4), np.complex128)
    for iG in range(16):
        res[iG] = np.tensordot(
                    bkw_mat_h,
                    matrix_tensor_dot(Gamma16_dgr[iG], frw_mat, 1),
                    axes=([2,3],[0,1]))
    return frw_mat, bkw_mat, res

def make_test_res_p(p):
    return make_test_res(p['latsize'], p['frw_lin_coeff'], p['frw_ft_mom'],
                         p['bkw_lin_coeff'], p['bkw_ft_mom'], p['ft_x0'])

def check_output_data(p, test_dir='.'):
    aff_r   = aff.Reader(test_dir + '/' + p['aff_file'])

    frw_p   = np.array(aff_r.read(p['frw_kpath'])).reshape((3,4,3,4))
    bkw_p   = np.array(aff_r.read(p['bkw_kpath'])).reshape((3,4,3,4))
    qv      = np.empty((16, 3,4,3,4), np.complex128)
    for iG in range(16):
        d       = aff_r.read(p['qv_kpath'] + ('/g%d' % iG))
        qv[iG]  = np.array(d).reshape((3,4,3,4))

    frw_p1, bkw_p1, qv1 = make_test_res_p(p)
    cmp_data    = (np.allclose(frw_p, frw_p1) 
                    and np.allclose(bkw_p, bkw_p1) 
                    and np.allclose(qv, qv1))

    #print y.shape, y1.shape

    dqv = qv - qv1
    print ("%e\t%e\t%e" % (cplx_norm2(frw_p), cplx_norm2(frw_p1), rdiff(frw_p, frw_p1)))
    print ("%e\t%e\t%e" % (cplx_norm2(bkw_p), cplx_norm2(bkw_p1), rdiff(bkw_p, bkw_p1)))
    print ("%e\t%e\t%e" % (cplx_norm2(qv), cplx_norm2(qv1), rdiff(qv, qv1)))
    
    return cmp_data


def check_qlua_run(p, test_dir='.'):
    f = open(test_dir + '/stdout', 'r')
    for l in f:
        if (l == 'QLUA_RUN_SUCCESS\n'):
            return True
    return False

case_list = [
    {   'case'      : 'save-npr-test1',
        'qlua_src'  : 'save-npr/test1.qlua',
        'make_data' : make_test_res_p,
        'check_log' : check_qlua_run,
        'check_out' : check_output_data }
]
input_list = [
   {'aff_file'          : 'test.aff',
    'qv_kpath'      : '/test1/qvertex',
    'latsize'           : [6,6,6,12],
    'frw_kpath'         : '/test1/frw_prop',
    'frw_ft_mom'        : [1,2,3,4],
    'frw_lin_coeff'     : [1., .3+.2j, .2+.3j, -.4+.5j, .25-.17j],
    'bkw_kpath'         : '/test1/bkw_prop',
    'bkw_ft_mom'        : [0,1,-1,1],
    'bkw_lin_coeff'     : [2., .13+.12j, .42+.23j, -.04+.05j, .5-.7j],
    'ft_x0'             : [0,1,2,3] }
]
