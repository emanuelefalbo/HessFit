#!/usr/bin/env python3

import sys
import numpy as np
from numpy.linalg import inv
from scipy import optimize
import matplotlib.pyplot as plt

def load_data(fname):
    data = []
    with open(fname,'r') as f:
        for line in f:
            data.append(line.split())
        data = np.array(data, dtype=float)    
    return data

def build_A_matrix(xval):
    deg2rad = np.pi/180
    xval = xval*deg2rad
    cos_1 =  (1 + np.cos(xval)  )
    cos_2 =  (1 - np.cos(2*xval))
    cos_3 =  (1 + np.cos(3*xval))
    cos_4 =  (1 - np.cos(4*xval))
    A_matrix = np.vstack((cos_1, cos_2, cos_3, cos_4)).T
    return A_matrix

def get_rmsd(y_pred, y_data):
    MSE = np.square(np.subtract(y_pred, y_data)).mean()
    RMSD = np.sqrt(MSE)
    print('   MSE   = {:.5f} '.format(MSE) )
    print(' % RMSD  = {:.5f} '.format(RMSD*100) )

def print_output(fname, coeffs, A_matrix, qm_rel, mm_rel):
    yfit = np.array(
                     coeffs[0]*A_matrix[:,0] +
                     coeffs[1]*A_matrix[:,1] +
                     coeffs[2]*A_matrix[:,2] +
                     coeffs[3]*A_matrix[:,3] +
                     mm_rel)

    with open(fname, 'w') as f:
        for i, item in enumerate(yfit):
            print(i, qm_rel[i], mm_rel[i], item, file = f, sep='  ')
    f.close()

    return yfit


#####################
# NON LINEAR FITTING
#####################

def non_linear_f(x, a1, a2, a3, a4):
    y = a1*(1 + np.cos(x)) + \
        a2*(1 - np.cos(2*x)) + \
        a3*(1 - np.cos(3*x)) + \
        a4*(1 - np.cos(4*x)) 
    return y

def residual(p, x, y):
    return y - non_linear_f(x, *p)

def main():
    fname = sys.argv[1]
    data = load_data(fname)
    xval = data[:,0]
    qm_rel = data[:,1] - min(data[:,1])
    mm_rel = data[:,2] - min(data[:,2])
    yval = qm_rel - mm_rel
    A_mat = build_A_matrix(xval)
    
    # LLSQ 
    coeffs = np.dot(np.linalg.inv(np.dot(A_mat.T, A_mat)), np.dot(A_mat.T, yval))
    # cov_matrix = np.diag(residual(coeffs, xval, yval)**2)
    print(coeffs)
    print('', 60*'=')
    print(' Torsional Coefficients after LLSQ (kcal/mol): ')
    print('', 60*'=')
    [ print(' Coefficient # {} = {:.4f}'.format(i+1, x*6.275030E02)) for i, x in enumerate(coeffs)]
    fname = 'linear_fitted.txt'
    yfit_linear = print_output(fname, coeffs, A_mat, qm_rel, mm_rel)
    get_rmsd(yfit_linear, yval)

     # WLLSQ 
    kT = 3.167*1E0-6* 1298.15 # Boltzman constant in au at Room Temperature
    w = np.exp(-yval/kT)
    # print(w)
    W_matrix = np.diag(w)
    AW = np.dot(A_mat.T, W_matrix)
    prefac = np.dot(AW, A_mat)
    fac = np.dot(AW, yval)
    w_coeffs = np.dot(np.linalg.inv(prefac), fac)
    print(w_coeffs)
    print('', 60*'=')
    print(' Torsional Coefficients after WLLSQ (kcal/mol): ')
    print('', 60*'=')
    [ print(' Coefficient # {} = {:.4f}'.format(i+1, x*6.275030E02)) for i, x in enumerate(w_coeffs)]
    fname = 'w_linear_fitted.txt'
    yfit_w_linear = print_output(fname, w_coeffs, A_mat, qm_rel, mm_rel)
    get_rmsd(yfit_w_linear, yval)

    # # NON-LLSQ
    # # print('', 60*'=')
    # # print(' Torsional Coefficients after Non-linear LLSQ (kcal/mol): ')
    # # print('', 60*'=')
    # # # p0 = [0.1, 0.1, 0.1, 0.1]
    # # p0 = [1., 1., 1., 1.]
    # # popt, pcov = optimize.leastsq(residual, p0, args=(xval, yval))
    # # [ print(' Coefficient # {} = {:.4f}'.format(i+1, x*6.275030E02)) for i, x in enumerate(popt)]
    # # fname = 'non-linear_fitted.txt'
    # # yfit_non_linear = print_output(fname, popt, A_mat, qm_rel, mm_rel)
    # # get_rmsd(yfit_non_linear, yval)


    # # PLOT DATA
    plt.plot(xval, qm_rel, label='QM', marker='x')
    plt.plot(xval, mm_rel, label='MM', marker='d')
    plt.plot(xval, yfit_linear, label='LLSQ', marker='o')
    plt.xlabel('Angle (°)')
    plt.ylabel('$\Delta$E (kcal/mol)')
    # plt.plot(xval, yfit_w_linear, label='WLLSQ', marker='v')
    # plt.plot(xval, yfit_non_linear, label='NON-LLSQ', marker='v')
    plt.legend()
    plt.show()



if __name__ == '__main__':
    main()
