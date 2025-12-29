#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def build_parser():
    parser = argparse.ArgumentParser(description="Fit torsional potential to QM data using LLSQ, WLLSQ, or Nonlinear methods.")
    parser.add_argument('file', type=str, help='Input file with angle, QM, and MM energies')
    parser.add_argument('--model', choices=['oplsa', 'rybe'], default='oplsa', help='Torsional potential model')
    parser.add_argument('--plot', action='store_true', help='Plot the energy surface')
    parser.add_argument('--weighted', action='store_true', help='Use weighted linear least squares')
    parser.add_argument('--nonlinear', action='store_true', help='Use nonlinear least squares')
    parser.add_argument('--name', default='scan', help='Output plot filename (without extension)')
    return parser


def load_data(fname):
    data = pd.read_csv(fname, delim_whitespace=True, header=None)
    x = data.iloc[:, 0].to_numpy()
    qm_abs = data.iloc[:, 1].to_numpy()
    mm_abs = data.iloc[:, 2].to_numpy()
    qm_rel = (qm_abs - np.min(qm_abs)) * 627.5030
    mm_rel = (mm_abs - np.min(mm_abs)) #* 627.5030
    return x, qm_rel, mm_rel


def oplsa_matrix(xval, phase=None):
    xrad = np.deg2rad(xval)
    if phase is None:
        phase = np.zeros(4)
    terms = [np.ones_like(xrad)]
    terms += [(1 + (-1)**(i+1) * np.cos((i+1) * xrad + phase[i])) for i in range(4)]
    return np.column_stack(terms)

def oplsa_matrix_4term(xval, phase):
    xrad = np.deg2rad(xval)
    return np.column_stack([
        (1 + np.cos(xrad + phase[0])),
        (1 - np.cos(2*xrad + phase[1])),
        (1 - np.cos(3*xrad + phase[2])),
        (1 - np.cos(4*xrad + phase[3]))
    ])

def rybe_matrix(xval):
    xrad = np.deg2rad(xval)
    cos_phi = np.cos(xrad)
    terms = [np.ones_like(cos_phi)]
    for i in range(1, 6):
        terms.append((-1)**i * cos_phi**i)
    return np.column_stack(terms)


def solve_llsq(A, yval, weights=None):
    if weights is not None:
        W = np.diag(weights)
        coeffs = np.linalg.lstsq(W @ A, W @ yval, rcond=None)[0]
    else:
        coeffs = np.linalg.lstsq(A, yval, rcond=None)[0]
    return coeffs


def get_fit(A, coeffs, mm_rel):
    return A @ coeffs + mm_rel


def get_rmsd(yfit, yref, label=''):
    mse = np.mean((yfit - yref)**2)
    rmsd = np.sqrt(mse)
    print(f'{label} RMSD = {rmsd:.5f} kcal/mol')
    return rmsd


def print_coefficients(coeffs, label=''):
    print(f'\n{label} Torsional Coefficients (kcal/mol):')
    for i, c in enumerate(coeffs[1:], 1):
        print(f'  Coefficient #{i}: {c:.4f}')


def non_linear_function(x, a1, a2, a3, a4, p1, p2, p3, p4):
    xrad = np.deg2rad(x)
    return (
        a1 * (1 + np.cos(xrad + p1)) +
        a2 * (1 - np.cos(2 * xrad + p2)) +
        a3 * (1 - np.cos(3 * xrad + p3)) +
        a4 * (1 - np.cos(4 * xrad + p4))
    )


def plot_results(xval, qm, mm, yfit, title, outname):
    plt.figure(figsize=(10, 7))
    plt.plot(xval, qm, 'ko--', label='QM')
    plt.plot(xval, mm, 'ro--', label='MM')
    plt.plot(xval, yfit, 'bo-', label=title)
    plt.xlabel('Angle (°)', fontsize=16)
    plt.ylabel('Energy (kcal/mol)', fontsize=16)
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{outname}.png')
    plt.show()


def main():
    args = build_parser().parse_args()
    xval, qm, mm = load_data(args.file)
    yval = qm - mm

    if args.model == 'oplsa':
        A = oplsa_matrix(xval)
    else:
        A = rybe_matrix(xval)

    # Linear Least Squares (with optional weighting)
    weights = None
    if args.weighted:
        weights = np.exp(-0.2 * np.sqrt(qm))

    coeffs = solve_llsq(A, yval, weights)
    fit_label = 'WLLSQ' if args.weighted else 'LLSQ'
    print_coefficients(coeffs, label=fit_label)
    yfit = get_fit(A, coeffs, mm)
    get_rmsd(yfit, qm, label=fit_label)

    # Save result
    pd.DataFrame({'Angle': xval, 'QM': qm, 'MM': mm, 'Fitted': yfit}).to_csv(f'{args.name}_fitted.csv', index=False)

    # Nonlinear fit (optional)
    if args.nonlinear and args.model == 'oplsa':
        p0 = np.zeros(8)
        popt, _ = curve_fit(non_linear_function, xval, yval, p0=p0)
        amps, phases = popt[:4], popt[4:]
        # A_shifted = oplsa_matrix(xval, phase=phases)
        A_shifted = oplsa_matrix_4term(xval, phases)
        yfit_nonlin = get_fit(A_shifted, amps, mm)
        print_coefficients(np.concatenate([[0], amps]), label='Nonlinear')
        get_rmsd(yfit_nonlin, qm, label='Nonlinear')
        pd.DataFrame({'Angle': xval, 'QM': qm, 'MM': mm, 'Fitted': yfit_nonlin}).to_csv(f'{args.name}_nonlinear.csv', index=False)

    # Plot
    if args.plot:
        title = 'Nonlinear Fit' if args.nonlinear else fit_label
        yplot = yfit_nonlin if args.nonlinear and args.model == 'oplsa' else yfit
        plot_results(xval, qm, mm, yplot, title, args.name)


if __name__ == '__main__':
    main()
