import sys, os
import shutil
import math

def calculate_distance(str1, str2):
    x1, y1, z1 = str1.strip().split()
    x1 = float(x1)
    y1 = float(y1)
    z1 = float(z1)

    x2, y2, z2 = str2.strip().split()
    x2 = float(x2)
    y2 = float(y2)
    z2 = float(z2)

    dx = x1-x2
    dy = y1-y2
    dz = z1-z2

    distance = math.sqrt(dx**2 + dy**2 + dz**2)
    return distance

 
def exchange_coords(nearest_water_index, lines, line_number, mol_line_number):
    o1 = lines[line_number][21:45] 
    h1 = lines[line_number + 1][21:45]
    h2 = lines[line_number + 2][21:45]

    mol_o = lines[mol_line_number][21:45]
    mol_h1 = lines[mol_line_number + 1][21:45]
    mol_h2 = lines[mol_line_number + 2][21:45]

    lines[mol_line_number] = lines[mol_line_number][:21] + o1
    lines[mol_line_number + 1] = lines[mol_line_number + 1][:21] + h1
    lines[mol_line_number + 2] = lines[mol_line_number + 2][:21] + h2

    lines[line_number] = lines[line_number][:21] + mol_o
    lines[line_number + 1] = lines[line_number + 1][:21] + mol_h1
    lines[line_number + 2] = lines[line_number + 2][:21] + mol_h2
    return lines


def read_file(in_path):
    out_path = os.path.split(in_path)[0] + "1/" + os.path.split(in_path)[1]
    att = {}
    nearest_water = {}
    trapped_water = {}
    distances ={}
    trapped_water_distance = {}
    nearest_water_distance = 0.0
    nearest_water_index = 0
    line_number = 0
    mol_line_number = 0
    f = open(in_path, 'r')
    lines = f.readlines()
    for idx, line in enumerate(lines):
        if line[5:8] == "ATT":
            att['coord_string'] = line[21:45].strip()
        if line[5:8] == "MOL" and line[14:15] == "O":
            trapped_water['coord_string'] = line[21:45].strip()
            trapped_water_distance = calculate_distance(att['coord_string'], trapped_water['coord_string'])
            nearest_water_distance = trapped_water_distance
            mol_line_number = idx
        if (line[5:8] == "HOH" or line[5:8] == "SOL") and line[14:15] == "O":
            water_dist = calculate_distance(line[21:45].strip(), att['coord_string'])
            if water_dist < nearest_water_distance:
                nearest_water_distance = water_dist
                nearest_water_index = int(line[:5].strip())
                line_number = idx
    f.close()
    if nearest_water_index != 0:
        lines = exchange_coords(nearest_water_index, lines, line_number, mol_line_number)
        with open(out_path, 'w') as f:
            for line in lines:
                f.write(line)
        f.close()
        print(f"exchanged for {out_path}")
    else:
        shutil.copyfile(in_path, out_path)
    return

def read_folder(folder_path):
    for i in range(1, 101):
        read_file(f"{folder_path}/frame{i}.gro")
    return


def main():
    folder_path = sys.argv[1]
    read_folder(folder_path)


if __name__ == "__main__":
    main()
