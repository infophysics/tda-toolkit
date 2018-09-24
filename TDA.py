from tda.tda import CubicalRipser2D, CubicalRipser3D, Perseus, Ripser

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