#!/usr/bin/env python3


def print_GauHarm(*arg):
    ele_ls = arg[0]
    tp_ls= arg[1]
    coord = arg[2]
    bond_tp_ls = arg[3]
    k_bond_ls = arg[4]
    bond_ln_ls = arg[5]
    angle_tp_ls = arg[6]
    k_angle_ls = arg[7]
    angle_ln_ls = arg[8]
    chg = arg[9]

    header_gjf ="""%mem=1GB
%nprocshared=1
%chk=GauHarm.chk
#p Amber=(SoftOnly,Print) nosymm  Freq=intmodes
 
Title

0 1
"""

    master_func = """
! Master function
NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.2
"""
    fname = 'GauHarm.gjf'
    with open(fname, 'w') as fout:
        fout.write(header_gjf)
        for m, p, l in zip(ele_ls, tp_ls, coord):
                    #  s1 = '  '.join(str(x) for x in p)
                     s2 = '  '.join((f'{x:.6f}') for x in l)
                     fout.write(f'{m}-{p}-0.00  {s2}\n')
        # fout.write(f'\n')
        fout.write(master_func)
        fout.write(f'! SMARTFIELD FF\n')
        fout.write(f'!Bonds\n')
        for k in range(len(bond_tp_ls)):
            msg = (
                  f'HrmStr1 {bond_tp_ls[k]}  1.0 ' 
                  f' {bond_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'!Angles\n')
        for k in range(len(angle_tp_ls)):
            msg = (
                  f'HrmBnd1 {angle_tp_ls[k]}  1.0 ' 
                  f' {angle_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'\n')


def print_GauNonBon(*arg):
    ele_ls = arg[0]
    tp_ls= arg[1]
    coord = arg[2]
    bond_tp_ls = arg[3]
    bond_ln_ls = arg[4]
    angle_tp_ls = arg[5]
    angle_ln_ls = arg[6]
    chg = arg[7]
    vdw = arg[8]

    header_gjf ="""%mem=1GB
%nprocshared=1
%chk=GauNonBon.chk
#p Amber=(SoftOnly,Print) nosymm Freq=intmodes
 
Title

0 1
"""

    master_func = """
! Master function
NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.2
"""
    fname = 'GauNonBon.gjf'
    with open(fname, 'w') as fout:
        fout.write(header_gjf)
        for m, p, l, q in zip(ele_ls, tp_ls, coord, chg):
                #   s1 = '  '.join(str(x) for x in p)
                s2 = '  '.join((f'{x:.6f}') for x in l)
                fout.write(f'{m}-{p}-{q}  {s2}\n')
        # fout.write(f'\n')
        fout.write(master_func)
        fout.write(f'! SMARTFIELD FF\n')
        fout.write(f'!Bonds\n')
        for k in range(len(bond_tp_ls)):
            msg = (
                  f'HrmStr1 {bond_tp_ls[k]}  0.0 ' 
                  f' {bond_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'!Angles\n')
        for k in range(len(angle_tp_ls)):
            msg = (
                  f'HrmBnd1 {angle_tp_ls[k]}  0.0 ' 
                  f' {angle_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'!VDW\n')
        for k in vdw:
            s = ' '.join(k)
            fout.write(f'{s} \n')
        fout.write(f'\n')