import numpy as np


# Define a function to process each row
def process_row(row):
    # Split the row into individual values
    values = row.split()
    result = set()
    # for item in values:
    #     result.add(frozenset(lang for lang in item.split(" ") if lang))
    # print(result)


    # Print the values or perform any other processing
    print(values)

  
    
def read_structure_data(text):
    data = {"bonds": [], "angles": [], "torsions": []}
    current_section = None

    lines = text.split('\n')
    for line in lines:
        # Skip empty lines
        if not line.strip():
            continue

        # Detect section change
        if line.startswith('!'):
            current_section = line[1:].strip().lower()
            # print(current_section)
            continue

        # Parse data based on the current section
        if current_section:
            values = line.split()
            if current_section == 'bonds' and len(values) == 5:
                data["bonds"].append({
                    "type": values[0],
                    "atom1": values[1],
                    "atom2": values[2],
                    "value": float(values[3]),
                    "force_constant": float(values[4])
                })
            elif current_section == 'angles' and len(values) == 6:
                data["angles"].append({
                    "type": values[0],
                    "atom1": values[1],
                    "atom2": values[2],
                    "atom3": values[3],
                    "value": float(values[4]),
                    "force_constant": float(values[5])
                })
            elif current_section == 'torsions' and len(values) == 14:
                data["torsions"].append({
                    "type": values[0],
                    "atom1": values[1],
                    "atom2": values[2],
                    "atom3": values[3],
                    "atom4": values[4],
                    "phi_0": float(values[5]),
                    "k": float(values[6]),
                    "multiplicity": int(values[7]),
                    "periodicity": int(values[8]),
                    "phase": float(values[9]),
                    "E0": float(values[10]),
                    "E1": float(values[11]),
                    "E2": float(values[12]),
                    "E3": float(values[13])
                })

    return data


with open("SmartField4gau.gjf", 'r') as file:
    # Split the input text into rows
    # rows = file.strip().split('\n')
    content = file.read()
    # rows = content.split('\n')
    
    result = read_structure_data(content)
    
    tors_list =[]
    keys =['atom1', 'atom2', 'atom3', 'atom4']
    for i in result['torsions']:
        tors_list.append([i.get(key) for key in keys])
    print(tors_list)
    # print("")
    
    result_list = [' '.join(sublist) for sublist in tors_list]
    # [print(i,"; ", j.split()) for i,j in zip(result_list,sorted(result_list))]
    
    # def unique_list(mylist):
    #     x = set([tuple(sorted(s.split())) for s in mylist])
    #     return [" ".join(s) for s in x]
    
    # for i in unique_list(result_list):
    #     print(i)
    
    
    # tmp = set()
    # for item in result_list:
    #     tmp.add(frozenset(lang for lang in item.split(" ") if lang))

    # print(tmp)
    

         
    