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
        


def print_GauInp(*arg):
    ele_ls = arg[0]
    tp_ls= arg[1]
    coord = arg[2]
    bond_tp_ls = arg[3]
    k_bond_ls = arg[4]
    bond_ln_ls = arg[5]
    angle_tp_ls = arg[6]
    k_angle_ls = arg[7]
    angle_ln_ls = arg[8]
    tors_tp_ls = arg[9]
    v1 = arg[10]
    v2 = arg[11]
    v3 = arg[12]
    phase = arg[13]
    hybrid_list = arg[14]
    chg = arg[15]
    """
    Write a gaussian input file:
    ele_ls: element list 
    tp_ls: atom type list
    coord: coordinates XYZ
    bond_tp_ls:  bond type list
    bond_ln_ls:  bond length list 
    k_bond_ls:  bond constant list
    angle_tp_ls: angle type list 
    k_angle_ls: angle constant list 
    angle_ln_ls : deg angle list

    """
    header_gjf ="""%mem=1GB
%nprocshared=1
%chk=SmartField4gau.chk
#p Amber=(SoftFirst,Print) nosymm geom=nocrowd opt Freq
 
Title

0 1
"""

    master_func = """
! Master function
NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.2
"""
    fname = 'SmartField4gau.gjf'
    with open(fname, 'w') as fout:
        fout.write(header_gjf)
        for m, p, l, q in zip(ele_ls, tp_ls, coord, chg):
                    #  s1 = '  '.join(str(x) for x in p)
                     s2 = '  '.join((f'{x:.6f}') for x in l)
                     fout.write(f'{m}-{p}-{q}  {s2}\n')
        # fout.write(f'\n')
        fout.write(master_func)
        fout.write(f'! SMARTFIELD FF\n')
        fout.write(f'!Bonds\n')
        for k in range(len(bond_tp_ls)):
            msg = (
                  f'HrmStr1 {bond_tp_ls[k]}  {k_bond_ls[k]:.3f} ' 
                  f' {bond_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'!Angles\n')
        for k in range(len(angle_tp_ls)):
            msg = (
                  f'HrmBnd1 {angle_tp_ls[k]}  {k_angle_ls[k]:.3f} ' 
                  f' {angle_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'! Torsions\n')
        for k in range(len(tors_tp_ls)):
            s = '  '.join((f'{x}') for x in phase[k,:])
            msg = (
                  f'AmbTrs {tors_tp_ls[k]} {s} '
                  f' {v1[k]:.2f} {v2[k]:.2f} {v3[k]:.2f} 0. ' 
                  f' {float(hybrid_list[k])}\n'
                  )
            fout.write(msg)
        fout.write(f'\n')


def print_AmbFrcmod(*arg):
    tp_ls= arg[0]
    bond_tp_ls = arg[1]
    k_bond_ls = arg[2]
    bond_ln_ls = arg[3]
    angle_tp_ls = arg[4]
    k_angle_ls = arg[5]
    angle_ln_ls = arg[6]
    tors_tp_ls = arg[7]
    v1 = arg[8]
    v2 = arg[9]
    v3 = arg[10]
    phase=arg[11]
    hybrid_list = arg[12]
    
    # print(v1)

    """
    Write a AMBER-like frcmod file:
    tp_ls: atom type list
    bond_tp_ls:  bond type list
    bond_ln_ls:  bond length list 
    k_bond_ls:  bond constant list
    angle_tp_ls: angle type list 
    k_angle_ls: angle constant list 
    angle_ln_ls : deg angle list

    """
    tp_ls_unique = set(tp_ls)
    fname = 'SmartField_frcmod.txt'
    with open(fname, 'w') as fout:
        fout.write('MASS\n')
        for m in tp_ls_unique:
            fout.write(f'{m}\n')
        fout.write(f'BOND\n')
        for k in range(len(bond_tp_ls)):
            msg = (
                  f'{bond_tp_ls[k]} {k_bond_ls[k]:.3f} ' 
                  f' {bond_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'ANGLE\n')
        for k in range(len(angle_tp_ls)):
            msg = (
                  f'{angle_tp_ls[k]}  {k_angle_ls[k]:.3f} ' 
                  f' {angle_ln_ls[k]:.3f}\n'
                  )
            fout.write(msg)
        fout.write(f'DIHE\n')
        for k in range(len(tors_tp_ls)):
            s = '  '.join((f'{x}') for x in phase[k,:])
            msg = (
                  f'{tors_tp_ls[k]} {s}'
                  f' {v1[k]:.2f} {v2[k]:.2f} {v3[k]:.2f} 0.0 ' 
                  f' {hybrid_list[k]}\n'
                  )
            fout.write(msg)
        fout.write(f'\n')
    