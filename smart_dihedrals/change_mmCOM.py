#!/usr/bin/env python3


import sys

def change_mm_COM(filename, angle):
    with open(filename, 'r') as f:
        text = []
        for line in f:
            text.append(line.split())

    # Break it in parts
    first_bit = []
    second_bit = []
    third_bit = []
    for i in text:
        first_bit.append(i)
        if len(i) == 2:
            break
    for i in text[8:]:
        second_bit.append(i)
        if i[0] == 'Variables:':
            break

    second_bit[angle[-1]-1][1] = angle[2]  # Replace Bond
    second_bit[angle[-1]-1][3] = angle[1]  #  .... Angle
    second_bit[angle[-1]-1][5] = angle[0]  #  .. Dihedral

    r1 = second_bit[angle[-1]-1][2]        # Bond Paramters
    r2 = second_bit[angle[-1]-1][4]        # Angle ....
    r3 = second_bit[angle[-1]-1][6]        #  Dihedral ...

    for i in text[len(second_bit)+8:]:
        third_bit.append(i)

    third_bit = [x for x in third_bit if x]   # Take out empty spaces in list

    for i, x in enumerate(third_bit):
        if x[0] == r3:
            third_bit[i][1] = '0.000 S 15 24.0 '


    fout = filename[:-4] + '_mod.gjf'
    text_new = first_bit + second_bit + third_bit
    with open(fout, 'w') as f2:
        for i in text_new:
            print(*i, file=f2)
            # print(i, file = fout)


def main():
    filename = sys.argv[1]
    # angle = [2, 1, 5, 8]
    angle = [1, 5, 8, 11]
    change_mm_COM(filename, angle)

    

if __name__ == '__main__':
    main()
        
