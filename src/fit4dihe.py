#!/usr/bin/env python3

import argparse
import numpy as np
from numpy.linalg import inv
from scipy import optimize
import matplotlib.pyplot as plt
import pandas as pd


def build_parser():
    """
    Build options for parser
    """
    parser = argparse.ArgumentParser()
    txt = "file containg phi angle, QM, and MM energies in 1st, 2nd, and 3rd columns"
    parser.add_argument('file', type=str, help=txt)
    txt = 'Plot the PES with matplotlib'
    parser.add_argument('--plot', action="store_true", help=txt)
    txt = 'Name of output png file'
    parser.add_argument('--name', default='scan', help=txt)
    return parser


def load_data(fname):
    data = []
    with open(fname,'r') as f:
        for line in f:
            data.append(line.split())
        data = np.array(data, dtype=float)    
    return data

def oplsa_matrix(xval, phase):
    deg2rad = np.pi/180
    xval = xval*deg2rad
    if phase is None:
       phase = np.zeros(4)

    cos_0 = np.ones(len(xval))
    cos_1 =  (1 + np.cos(xval + phase[0] ) )
    cos_2 =  (1 - np.cos(2*xval + phase[1]))
    cos_3 =  (1 + np.cos(3*xval + phase[2]))
    cos_4 =  (1 - np.cos(4*xval + phase[3]))
    A_matrix = np.vstack((cos_0, cos_1, cos_2, cos_3, cos_4)).T
    # A_matrix = np.vstack((cos_1, cos_2, cos_3, cos_4)).T

    return A_matrix

def RyBe_matrix(xval, phase):
    deg2rad = np.pi/180
    xval = xval*deg2rad
    if phase is None:
       phase = np.zeros(6)

    cos_0 = np.ones(len(xval))
    cos_1 =  -np.cos(xval) 
    cos_2 = (np.cos(xval))**2
    cos_3 = -(np.cos(xval))**3 
    cos_4 = (np.cos(xval))**4
    cos_5 = -(np.cos(xval))**5
    A_matrix = np.vstack((cos_0, cos_1, cos_2, cos_3, cos_4, cos_5)).T

    return A_matrix

def get_rmsd(y_pred, y_data):
    MSE = np.square(np.subtract(y_pred, y_data)).mean()
    RMSD = np.sqrt(MSE)
    print('   MSE   = {:.5f} '.format(MSE) )
    print(' % RMSD  = {:.5f} '.format(RMSD*100) )

def get_fit(model, fname, coeffs, A_matrix, xval, qm_rel, mm_rel):
    if model == 0:                                   # 0 triggers on the OPLSA type 
       yfit = np.array(
                     coeffs[0]* A_matrix[:,0] +
                     coeffs[1]* A_matrix[:,1] +
                     coeffs[2]* A_matrix[:,2] +
                     coeffs[3]* A_matrix[:,3] +
                     coeffs[4]* A_matrix[:,4] +
                     mm_rel)
                    #  0.5* coeffs[0]* A_matrix[:,0] +
                    #  0.5* coeffs[1]* A_matrix[:,1] +
                    #  0.5* coeffs[2]* A_matrix[:,2] +
                    #  0.5* coeffs[3]* A_matrix[:,3] +
                    #  mm_rel)
    else:                                            # Otherwise Ryckart-Bellemans potential
        yfit = np.array( coeffs[0]* A_matrix[:,0] +
                         coeffs[1]* A_matrix[:,1] +
                         coeffs[2]* A_matrix[:,2] +
                         coeffs[3]* A_matrix[:,3] +
                         coeffs[4]* A_matrix[:,4] +
                         coeffs[5]* A_matrix[:,5] +
                         mm_rel)

    with open(fname, 'w') as f:
        for i, item in enumerate(yfit):
            f.write(f' {xval[i]}  {qm_rel[i]:.6f}  {mm_rel[i]:.6f}  {item:.6f} \n')
    f.close()

    return yfit
    
#####################
# NON LINEAR FITTING
#####################

def non_linear_f(x, a1, a2, a3, a4, p1, p2, p3, p4):
    y = a1*(1 + np.cos(x + p1)) + \
        a2*(1 - np.cos(2*x + p2)) + \
        a3*(1 - np.cos(3*x + p3)) + \
        a4*(1 - np.cos(4*x + p4)) 
    return y

def solve_LLSQ(A_matrix, yval):
    # coeffs = np.dot(np.linalg.inv(np.dot(A_matrix.T, A_matrix)), np.dot(A_matrix.T, yval))
    coeffs = np.linalg.lstsq(A_matrix, yval)[0]
    print(coeffs)
    # cov_matrix = np.diag(residual(coeffs, xval, yval)**2)
    print('', 60*'=')
    print(' Torsional Coefficients after LLSQ (kcal/mol): ')
    print('', 60*'=')
    [ print(' Coefficient # {} = {:.4f}'.format(i+1, x)) for i, x in enumerate(coeffs)]

    return coeffs


def rb2opls(coeffs):
    f4 = - coeffs[4] / 4
    f3 = - coeffs[3] / 2
    f2 = - coeffs[2] - coeffs[4] / 4
    f1 = - 2* coeffs[1] - (3/2)*coeffs[3]
    x = np.array([f1, f2, f3, f4])

    # A = np.array([[1/2, 1, 1/2, 0],
    #               [-1/2, 0, 3/2, 0],
    #               [0, -1, 0, 4],
    #               [0, 0, -2, 0]])
    # b = coeffs[:4]
    # x = np.linalg.solve(A, b)
    return x              
    # return c_opls

def main():
    parser = build_parser()
    args = parser.parse_args()
    fname = args.file
    data = load_data(fname)
    xval = data[:,0] # - 180.0
    qm_rel = ( data[:,1] - min(data[:,1]) ) *6.275030E02  # Eh to kcal/mol
    mm_rel = ( data[:,2] - min(data[:,2]) ) *6.275030E02
    yval =  ( qm_rel - mm_rel ) 
    pd.Series(qm_rel).to_csv('qm_rel.csv', index=False)
    pd.Series(mm_rel).to_csv('mm_rel.csv', index=False)
    
    A_oplsa = oplsa_matrix(xval, None)
    A_rybe = RyBe_matrix(xval, None)
    # LLSQ
    coeffs_oplsa = solve_LLSQ(A_oplsa, yval)
    coeffs_rybe = solve_LLSQ(A_rybe, yval)
    # c_oplsa_new = rb2opls(coeffs_rybe)
    # print(c_oplsa_new)
    fname_1 = 'oplsa_fitted.txt'
    fname_2 = 'rybe_fitted.txt'
    yfit_oplsa = get_fit(0, fname_1, coeffs_oplsa, A_oplsa, xval, qm_rel, mm_rel)
    yfit_rybe = get_fit(1, fname_2, coeffs_rybe, A_rybe, xval, qm_rel, mm_rel)
    # get_rmsd(yfit_oplsa, yval)
    # get_rmsd(yfit_rybe, yval)

    #  # WLLSQ 
    # kB = 8.617333262e-05
    # temp =  1000
    # kBT = kB * temp  # Boltzman constant in au at Room Temperature
    # # w = np.exp(-yval/kBT)
    # w = np.exp(-0.2 *np.sqrt(qm_rel))
    # print("")
    # print(w)
    # print("")

    # res = yval - np.dot(A_mat,coeffs)
    # # W_matrix = np.diag(w)
    # W_matrix = np.diag(res**2)
    # ATWI = np.dot(A_mat.T, np.linalg.inv(W_matrix))
    # prefac = np.dot(ATWI, A_mat)
    # fac = np.dot(ATWI, yval)
    # w_coeffs = np.dot(np.linalg.inv(prefac), fac)
    # print(w_coeffs)
    # print('', 60*'=')
    # print(' Torsional Coefficients after WLLSQ (kcal/mol): ')
    # print('', 60*'=')
    # [ print(' Coefficient # {} = {:.4f}'.format(i+1, x)) for i, x in enumerate(w_coeffs)]
    # fname = 'w_linear_fitted.txt'
    # A_mat_w = np.copy(A_mat)
    # yfit_w_linear = print_oplsTrs(fname, w_coeffs, A_mat_w, xval, qm_rel, mm_rel)
    # get_rmsd(yfit_w_linear, yval)


    # # NON-LLSQ
    # print('', 60*'=')
    # print(' Torsional Coefficients after Non-linear LLSQ (kcal/mol): ')
    # print('', 60*'=')
    # # p0 = [0.1, 0.1, 0.1, 0.1]
    # # p0 = np.ones(8)
    # p_0 = np.zeros(8)
    # # # popt, pcov = optimize.leastsq(residual, p0, args=(xval, yval))
    # popt, pcov = optimize.curve_fit(non_linear_f, xval, yval, p0=p_0)
    # # popt, pcov = optimize.curve_fit(residual, xval, yval)
    # [ print(' Coefficient # {} = {:.4f}'.format(i+1, x*6.275030E02)) for i, x in enumerate(popt)]
    # fname = 'non-linear_fitted.txt'
    # phase = popt[4:] * 180/np.pi
    # print(phase)
    # vs = popt[:4]
    # A_mat_shifted = oplstrs_matrix(xval, phase)
    # yfit_non_linear = print_oplsTrs(fname, vs, A_mat_shifted, xval, qm_rel, mm_rel)
    # get_rmsd(yfit_non_linear, yval)

    # # PLOT DATA
    if args.plot:
       fig = plt.subplots(1, figsize= (10,8), dpi=96 )
       plt.plot(xval, qm_rel, label='QM', marker='o', markersize=8, fillstyle='none', c='black', linestyle='dashed', dashes=(5, 10))
       plt.plot(xval, mm_rel, label='MM', marker='o', markersize=8, fillstyle='none', c='orangered', linestyle='dashed', dashes=(5, 10))
       plt.plot(xval, yfit_oplsa, label='Fitted$_{LLS}$', marker='o',  markersize=10, c='blue', alpha=0.6)
    #    plt.plot(xval, yfit_rybe, label='LLSQ$_{Ry-Be}$', marker='o',  markersize=10, c='green')
       # plt.plot(xval, yfit_w_linear, label='WLLSQ', marker='v')
       # plt.plot(xval, yfit_non_linear, label='NON-LLSQ', marker='v')
       plt.xlabel('Angle (°)', fontsize=18)
       plt.ylabel('$\Delta$E (kcal/mol)', fontsize=18)
       plt.xticks(fontsize=16)
       plt.yticks(fontsize=16)
       plt.legend(prop={'size':20})
       plt.tight_layout()
       plt.savefig(args.name +'.png')
       plt.show()



if __name__ == '__main__':
    main()
