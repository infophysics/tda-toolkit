from tda.tda import CubicalRipser2D, CubicalRipser3D, Perseus, Ripser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import struct

'''
    TDA package for Topological Data Analysis tools
'''

'''
    DIPHA Format file converter                             M. Tallon
    and converter to point cloud for RIPSER         (small edits by N. Carrara)

'''


def convert_csv_to_dipha(input_file, output_file):

    df = pd.read_csv(input_file, header=None)

    array_size = df.shape[1]
    iterations = df.shape[0]

    with open(output_file, 'wb') as f:
        symbol = struct.pack("<Q", 8067171840) # magic number
        f.write(symbol)
        symbol = struct.pack("<Q", 1) # 1 = Image file
        f.write(symbol)
        symbol = struct.pack("<Q", array_size * iterations) # image size
        f.write(symbol)
        symbol = struct.pack("<Q", 2) # 2 dimensions
        f.write(symbol)
        symbol = struct.pack("<Q", array_size) # image width
        f.write(symbol)
        symbol = struct.pack("<Q", iterations) # image height
        f.write(symbol)

        # add data for all points
        for row in range(df.shape[0]):
            for col in range(df.shape[1]):
                value = df.iloc[row, col]
                symbol = struct.pack("<d", value)
                f.write(symbol)
    print("Converted %s to DIPHA format file %s\n" % (input_file, output_file))


def save_array_to_dipha(input_array, output_file):

    df = pd.DataFrame(input_array)

    array_size = df.shape[1]
    iterations = df.shape[0]

    with open(output_file, 'wb') as f:
        symbol = struct.pack("<Q", 8067171840) # magic number
        f.write(symbol)
        symbol = struct.pack("<Q", 1) # 1 = Image file
        f.write(symbol)
        symbol = struct.pack("<Q", array_size * iterations) # image size
        f.write(symbol)
        symbol = struct.pack("<Q", 2) # 2 dimensions
        f.write(symbol)
        symbol = struct.pack("<Q", array_size) # image width
        f.write(symbol)
        symbol = struct.pack("<Q", iterations) # image height
        f.write(symbol)

        # add data for all points
        for row in range(df.shape[0]):
            for col in range(df.shape[1]):
                value = df.iloc[row, col]
                symbol = struct.pack("<d", value)
                f.write(symbol)
    print("Saved array to DIPHA format file %s\n" % output_file)


def convert_binary_cells_to_point_cloud(input_file, output_file):

    df = pd.read_csv(input_file, header=None)
    df = np.where(df==1)
    coords = np.asarray(list(zip(df[0],df[1])))
    dfOut = pd.DataFrame(coords)
    dfOut.to_csv(output_file, index=False, header=False)


def save_binary_cells_to_point_cloud(input_array, output_file):

    df = pd.DataFrame(input_array)
    df = np.where(df==1)
    coords = np.asarray(list(zip(df[0],df[1])))
    dfOut = pd.DataFrame(coords)
    dfOut.to_csv(output_file, index=False, header=False)


'''
    Persistent Homology 

'''


def plot_persistence_diagram(life_death_array):

    birth_times = [life_death_array[i][0] for i in range(len(life_death_array))]
    death_times = [life_death_array[i][1] for i in range(len(life_death_array))]

    max_birth = np.max(birth_times)
    max_death = np.max(death_times)
    max_both = max(max_birth, max_death)
    y = np.linspace(0, max_both, 2)

    plt.figure()
    plt.plot(y, y, color='g')
    plt.scatter(birth_times, death_times, color='b')
    plt.xlabel('Birth Time')
    plt.ylabel('Death Time')
    plt.title('Persistence Diagram')
    plt.grid(True)
    plt.show()


def plot_persistence_diagram_from_file(life_death_file):

    df = pd.read_csv(life_death_file, header=None)
    persist = [[df.values[i][1],df.values[i][2]] for i in range(len(df.values)) if df.values[i][0] != -1]
    plot_persistence_diagram(persist)