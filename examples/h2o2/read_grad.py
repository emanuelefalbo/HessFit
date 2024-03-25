import sys

def store_any_file(fname):
    """ 
    Reading Any file 
    """
    with open(fname, 'r') as f:
        all_lines = []
        for line in f:
            all_lines.append(line.strip())

    return all_lines

def range2(start, end):
     return range(start, end+1)

def read_gradient_file(all_lines):
    gradient=[]
    with open(file_path, 'r') as file:
        for s in range(len(all_lines)): 
            if "Cartesian Gradient" in all_lines[s]:# Get no of Atoms
                N_force = int(all_lines[s].split()[4])
                nlines = int(N_force/5.0) + 1
                start = s
                for e in range2(start + 1, start + nlines):
                    gradient.append(all_lines[e].split() )

    return gradient

# Example usage:
file_path = sys.argv[1]  # Provide the path to your gradient file
text = store_any_file(file_path)
gradient = read_gradient_file(text)
print("Gradient:")
for row in gradient:
    print(row)

