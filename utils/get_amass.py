# Dictionary of all elements matched with their atomic masses.
elements_dict = {'H' : 1.008,'HE' : 4.003, 'LI' : 6.941, 'BE' : 9.012,\
                 'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
                 'F' : 18.998, 'NE' : 20.180, 'NA' : 22.990, 'MG' : 24.305,\
                 'AL' : 26.982, 'SI' : 28.086, 'P' : 30.974, 'S' : 32.066,\
                 'CL' : 35.453, 'AR' : 39.948, 'K' : 39.098, 'CA' : 40.078,\
                 'SC' : 44.956, 'TI' : 47.867, 'V' : 50.942, 'CR' : 51.996,\
                 'MN' : 54.938, 'FE' : 55.845, 'CO' : 58.933, 'NI' : 58.693,\
                 'CU' : 63.546, 'ZN' : 65.38, 'GA' : 69.723, 'GE' : 72.631,\
                 'AS' : 74.922, 'SE' : 78.971, 'BR' : 79.904, 'KR' : 84.798,\
                 'RB' : 84.468, 'SR' : 87.62, 'Y' : 88.906, 'ZR' : 91.224,\
                 'NB' : 92.906, 'MO' : 95.95, 'TC' : 98.907, 'RU' : 101.07,\
                 'RH' : 102.906, 'PD' : 106.42, 'AG' : 107.868, 'CD' : 112.414,\
                 'IN' : 114.818, 'SN' : 118.711, 'SB' : 121.760, 'TE' : 126.7,\
                 'I' : 126.904, 'XE' : 131.294, 'CS' : 132.905, 'BA' : 137.328,\
                 'LA' : 138.905, 'CE' : 140.116, 'PR' : 140.908, 'ND' : 144.243,\
                 'PM' : 144.913, 'SM' : 150.36, 'EU' : 151.964, 'GD' : 157.25,\
                 'TB' : 158.925, 'DY': 162.500, 'HO' : 164.930, 'ER' : 167.259,\
                 'TM' : 168.934, 'YB' : 173.055, 'LU' : 174.967, 'HF' : 178.49,\
                 'TA' : 180.948, 'W' : 183.84, 'RE' : 186.207, 'OS' : 190.23,\
                 'IR' : 192.217, 'PT' : 195.085, 'AU' : 196.967, 'HG' : 200.592,\
                 'TL' : 204.383, 'PB' : 207.2, 'BI' : 208.980, 'PO' : 208.982,\
                 'AT' : 209.987, 'RN' : 222.081, 'FR' : 223.020, 'RA' : 226.025,\
                 'AC' : 227.028, 'TH' : 232.038, 'PA' : 231.036, 'U' : 238.029,\
                 'NP' : 237, 'PU' : 244, 'AM' : 243, 'CM' : 247, 'BK' : 247,\
                 'CT' : 251, 'ES' : 252, 'FM' : 257, 'MD' : 258, 'NO' : 259,\
                 'LR' : 262, 'RF' : 261, 'DB' : 262, 'SG' : 266, 'BH' : 264,\
                 'HS' : 269, 'MT' : 268, 'DS' : 271, 'RG' : 272, 'CN' : 285,\
                 'NH' : 284, 'FL' : 289, 'MC' : 288, 'LV' : 292, 'TS' : 294,\
                 'OG' : 294}

# List of all elements to allow for easy inclusion testing.
elements_list = ['H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA',\
                 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA', 'SC', 'TI',\
                 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE',\
                 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO',\
                 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE',\
                 'I', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM',\
                 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF',\
                 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB',\
                 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U',\
                 'NP', 'PU', 'AM', 'CM', 'BK', 'CT', 'ES', 'FM', 'MD', 'NO',\
                 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG', 'CN',\
                 'NH', 'FL', 'MC', 'LV', 'TS', 'OG']


print('Atomic Mass Calculator')
print('Developed by: Lukas Richters')
print('2017')
print()
print("Enter the formula for the substance by seperating the atomic symbols\
 from their subscripts. ie. C 6 H 12 O 6.")

# Get the substance's formula from user.
formula = input("Enter a formula or type 'q' to quit: ").upper().split()
print()

# Continue to ask for a new formula until the user inputs a "q" for quit.
while formula != ["Q"]:
    
    # Set the conditons for a new loop.
    atomic_mass_float = 0.0
    invalid_input = False
    coefficient = 1
    
    # Loop through each index in the formula.
    for i, ch in enumerate(formula):
        # If the character is an element, check if the next character in the
        # formula is an integer.
        if ch in elements_list:
            element_mass = elements_dict.get(ch)
            # If the next character is an integer, multiply the element's mass 
            # by that integer and add it to the running sum.
            if formula[i + 1].isdigit() == True:
                atomic_mass_float += element_mass * int(formula[i + 1])
            # If not, just add that element's mass to the sum.
            else:
                atomic_mass_float += element_mass
        # If the charcter is an integer
        elif ch.isdigit() == True:
            # If the first index is an integer, assume that that is a
            # coefficient for the formula and multiply the formula by that
            # number
            if i == 0:
                coefficient = int(ch)
            # Make sure the previous index wasn't an integer as well
            # (the integer needs an element to multiply), if it was, let the
            # program know that there is an error.
            elif formula[i - 1].isalpha() == False:
                invalid_input = True
            else:
                pass
        # If the only character in the formula is an integer, let the program
        # know that there is an error.
        elif ch.isdigit() == True and len(formula) == 1:
            invalid_input = True
        # If a character is not an element or an integer, let the program know
        # that there is an error.
        else:
            invalid_input = True
    
    # Muliply the entire atomic mass by the coefficient.
    atomic_mass_float *= coefficient
    
    # If there was an error in the program, display the appropriate error.
    if invalid_input == True:
        # Error for if there is only an integer entered.
        if ch.isdigit() == True and len(formula) == 1:
            print("You must enter at least one element for every integer\
 subscript.")
            print()
        # Error for if there is a float or non-elemnt entered.
        else:
            print("Please enter only the atomic symbols of elements or an integer\
 subscript.")
            print()
    
    # If there was no error, then print the claculated atomic mass
    else:
        print("Atomic Mass: ", round(atomic_mass_float, 3))
        print()
    
    # Ask for another formula.
    formula = input("Enter a formula or type 'q' to quit: ").upper().split()
    print()        
