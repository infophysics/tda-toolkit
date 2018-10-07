from tda.tda import CubicalRipser2D, CubicalRipser3D, Perseus, Ripser, Filter2D, BottleneckDistance, Generator

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import struct
import csv
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
import scipy.interpolate as interp
import ipywidgets as widgets
from IPython.display import display
import warnings
warnings.filterwarnings('ignore')


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
        symbol = struct.pack("<Q", 8067171840)  # magic number
        f.write(symbol)
        symbol = struct.pack("<Q", 1)  # 1 = Image file
        f.write(symbol)
        symbol = struct.pack("<Q", array_size * iterations)  # image size
        f.write(symbol)
        symbol = struct.pack("<Q", 2)  # 2 dimensions
        f.write(symbol)
        symbol = struct.pack("<Q", array_size)  # image width
        f.write(symbol)
        symbol = struct.pack("<Q", iterations)  # image height
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
        symbol = struct.pack("<Q", 8067171840)  # magic number
        f.write(symbol)
        symbol = struct.pack("<Q", 1)  # 1 = Image file
        f.write(symbol)
        symbol = struct.pack("<Q", array_size * iterations)  # image size
        f.write(symbol)
        symbol = struct.pack("<Q", 2)  # 2 dimensions
        f.write(symbol)
        symbol = struct.pack("<Q", array_size)  # image width
        f.write(symbol)
        symbol = struct.pack("<Q", iterations)  # image height
        f.write(symbol)

        # add data for all points
        for row in range(df.shape[0]):
            for col in range(df.shape[1]):
                value = df.iloc[row, col]
                symbol = struct.pack("<d", value)
                f.write(symbol)
    print("Saved array to DIPHA format file %s\n" % output_file)


#   Converting Binary Cells to point cloud (M. Tallon)

def convert_binary_cells_to_point_cloud(input_file, output_file):
    df = pd.read_csv(input_file, header=None)
    df = np.where(df == 1)
    coords = np.asarray(list(zip(df[0], df[1])))
    dfOut = pd.DataFrame(coords)
    dfOut.to_csv(output_file, index=False, header=False)


def save_binary_cells_to_point_cloud(input_array, output_file):
    df = pd.DataFrame(input_array)
    df = np.where(df == 1)
    coords = np.asarray(list(zip(df[0], df[1])))
    dfOut = pd.DataFrame(coords)
    dfOut.to_csv(output_file, index=False, header=False)


'''
    Persistent Homology 

'''


#   Persistence and Barcode plotting (N. Carrara)

def plot_persistence_diagram(barcode, split=True):
    dims = [barcode[i][0] for i in range(len(barcode))]
    #  find unique dimensions    
    unique_dims = []
    for dim in dims:
        if dim not in unique_dims:
            unique_dims.append(dim)
    if len(unique_dims) == 1:
        fig, axs = plt.subplots(1)
        max_val = np.max(barcode[:][:])
        y = np.linspace(0, max_val, 2)

        axs.plot(y, y, color='g', linestyle='--')
        for j in range(len(unique_dims)):
            birth_times = [barcode[i][1] for i in range(len(barcode)) if barcode[i][0] == unique_dims[j]]
            death_times = [barcode[i][2] for i in range(len(barcode)) if barcode[i][0] == unique_dims[j]]
        axs.scatter(birth_times, death_times)
        axs.set_xlabel('Birth Time')
        axs.set_ylabel('Death Time')
        axs.set_title('Persistence Diagram for degree $H_%s$' % 0)
        axs.grid(True)
    else:
        #   Plot persistence diagrams for each degree separately
        if split:
            fig, axs = plt.subplots(1, len(unique_dims), figsize=(15, 5))
            for j in range(len(unique_dims)):
                birth_times = [barcode[i][1] for i in range(len(barcode)) if barcode[i][0] == unique_dims[j]]
                death_times = [barcode[i][2] for i in range(len(barcode)) if barcode[i][0] == unique_dims[j]]

                max_birth = np.max(birth_times)
                max_death = np.max(death_times)
                max_both = max(max_birth, max_death)
                y = np.linspace(0, max_both, 2)

                axs[j].plot(y, y, color='g', linestyle='--')
                axs[j].scatter(birth_times, death_times, color='b')
                axs[j].set_xlabel('Birth Time')
                axs[j].set_ylabel('Death Time')
                axs[j].set_title('Persistence Diagram for degree $H_%s$' % j)
                axs[j].grid(True)
        #   Or together
        else:
            fig, axs = plt.subplots(1, figsize=(15,10))
            max_val = np.max(barcode[:][:])
            y = np.linspace(0, max_val, 2)

            axs.plot(y, y, color='g', linestyle='--')
            for j in range(len(unique_dims)):
                birth_times = [barcode[i][1] for i in range(len(barcode)) if barcode[i][0] == unique_dims[j]]
                death_times = [barcode[i][2] for i in range(len(barcode)) if barcode[i][0] == unique_dims[j]]
                axs.scatter(birth_times, death_times, label='$H_%s$' % j)
            axs.set_xlabel('Birth Time')
            axs.set_ylabel('Death Time')
            axs.legend()
            axs.set_title('Persistence Diagram')
            axs.grid(True)

    plt.show()


def plot_persistence_diagram_from_file(life_death_file, split=True):
    df = pd.read_csv(life_death_file, header=None)
    persist = [[df.values[i][0], df.values[i][1], df.values[i][2]] for i in range(len(df.values))]
    plot_persistence_diagram(persist, split=split)


def plot_barcode_diagram(barcode):
    dims = [barcode[i][0] for i in range(len(barcode))]
    #  find unique dimensions
    unique_dims = []
    for dim in dims:
        if dim not in unique_dims:
            unique_dims.append(dim)
    fig, axs = plt.subplots(1, len(unique_dims), figsize=(15, 5))
    if len(unique_dims) == 1:
        for j in range(len(unique_dims)):
            birth_times = [barcode[i][1] for i in range(len(barcode)) if barcode[i][0] == unique_dims[j]]
            death_times = [barcode[i][2] for i in range(len(barcode)) if barcode[i][0] == unique_dims[j]]

            max_birth = np.max(birth_times)
            max_death = np.max(death_times)
            min_birth = np.min(birth_times)
            max_both = max(max_birth, max_death)
            scale = 1.0 / max_death
            for k in range(len(birth_times)):
                axs.plot([birth_times[k], death_times[k]], [(k + 1), (k + 1)], color='b', linewidth=10)
            axs.set_xlabel('Time')
            axs.set_ylabel('Component')
            axs.set_ylim(0, len(birth_times) + 1)
            axs.set_yticks([])
            axs.set_title('Barcode Diagram for degree $H_%s$' % 0)
            axs.grid(True)
    else:
        for j in range(len(unique_dims)):
            birth_times = [barcode[i][1] for i in range(len(barcode)) if barcode[i][0] == unique_dims[j]]
            death_times = [barcode[i][2] for i in range(len(barcode)) if barcode[i][0] == unique_dims[j]]

            max_birth = np.max(birth_times)
            max_death = np.max(death_times)
            min_birth = np.min(birth_times)
            max_both = max(max_birth, max_death)
            scale = 1.0 / max_death
            for k in range(len(birth_times)):
                axs[j].plot([birth_times[k], death_times[k]], [(k + 1), (k + 1)], color='b', linewidth=10)
            axs[j].set_xlabel('Time')
            axs[j].set_ylabel('Component')
            axs[j].set_ylim(0, len(birth_times) + 1)
            axs[j].set_yticks([])
            axs[j].set_title('Barcode Diagram for degree $H_%s$' % j)
            axs[j].grid(True)

    plt.show()


def plot_barcode_diagram_from_file(life_death_file):
    df = pd.read_csv(life_death_file, header=None)
    persist = [[df.values[i][0], df.values[i][1], df.values[i][2]] for i in range(len(df.values))]
    plot_barcode_diagram(persist)


'''
    Persistent Homology Dimension                   (M. Tallon)
'''


def compute_2DPHD(barcode, show_plot=True, output_file=''):
    dfOut = pd.DataFrame([], columns=['PH Dim'])
    df = pd.DataFrame(barcode)
    df = df.dropna()
    #df.D.astype('float')
    df["X"] = (df.iloc[:, 2].values + df.iloc[:, 1].values) / 2
    df["Y"] = np.arccos(df.iloc[:, 1].values / df.iloc[:, 2].values)
    dfDim1 = df.loc[df.iloc[:, 0] == 1]
    dfDim0 = df.loc[df.iloc[:, 0] == 0]
    if dfDim1.shape[0] == 0:
        dfOut.append({'PH Dim': 'NA'}, ignore_index=True)
        return
    xSorted = np.sort(dfDim1.iloc[:, 3].values)
    logPH = np.empty((dfDim1.shape[0], 2))
    logPH[:, 0] = xSorted
    logPH[0, 1] = logPH.shape[0]
    eqCount = 1
    for row in range(1, logPH.shape[0]):
        if logPH[row, 0] > logPH[row - 1, 0]:
            logPH[row, 1] = logPH[row - 1, 1] - eqCount
            eqCount = 1
        else:
            logPH[row, 1] = logPH[row - 1, 1]
            eqCount += 1
    logPH = np.log10(logPH)
    slope = np.polyfit(logPH[:, 0], logPH[:, 1], 1)
    print(slope)
    dfOut = dfOut.append({'PH Dim': -slope[0]}, ignore_index=True)

    if output_file:
        dfOut.to_csv(output_file + ".csv", index=False)

    # Plots
    if show_plot:
        plt.figure()
        plt.tight_layout(pad=4.4, w_pad=4.5, h_pad=4.0)
        fitline = np.empty(logPH.shape)
        fitline[:, 0] = logPH[:, 0]
        fitline[:, 1] = fitline[:, 0] * slope[0] + slope[1]
        plt.plot(fitline[:, 0], fitline[:, 1])
        plt.scatter(logPH[:, 0], logPH[:, 1])
        plt.title("PH Dimension:  {}".format(str(round(-slope[0], 2))))
        plt.xlabel("Log(X)")
        plt.ylabel("Log(F(X))")
        plt.savefig((output_file + ".png"))
        #plt.show()


def compute_2DPHD_from_file(input_file, show_plot=True, output_file=''):
    df = pd.read_csv(input_file, header=None)
    persist = [[df.values[i][0], df.values[i][1], df.values[i][2]] for i in range(len(df.values))]
    compute_2DPHD(persist, show_plot, output_file)





'''
    Sliding Window Embedding                       (Chris Tralie)
                                                Modified by N. Carrara
'''



def get_sliding_window_1D(x, dim, Tau, dT):
    N = len(x)
    NWindows = int(np.floor((N-dim*Tau)/dT)) # The number of windows
    if NWindows <= 0:
        print("Error: Tau too large for signal extent")
        return np.zeros((3, dim))
    X = np.zeros((NWindows, dim)) # Create a 2D array which will store all windows
    idx = np.arange(N)
    for i in range(NWindows):
        # Figure out the indices of the samples in this window
        idxx = dT*i + Tau*np.arange(dim) 
        start = int(np.floor(idxx[0]))
        end = int(np.ceil(idxx[-1]))+2
        if end >= len(x):
            X = X[0:i, :]
            break
        # Do spline interpolation to fill in this window, and place
        # it in the resulting array
        X[i, :] = interp.spline(idx[start:end+1], x[start:end+1], idxx)
    return X

def embedding_transform_1D(X):
    pca = PCA(n_components = 2)
    Y = pca.fit_transform(X)
    eigs = pca.explained_variance_
    return Y, eigs

def plot_sliding_window_1D(x, X, Y, Tau, dT):
    NWindows = len(X)
    dim = len(X[0])
    N = len(x)
    extent = Tau*dim
    plt.figure(figsize=(9.5, 3))
    ax = plt.subplot(121)
    ax.plot(x)
    ax.set_ylim((-2*max(x), 2*max(x)))
    ax.set_title("Original Signal")
    ax.set_xlabel("Sample Number")
    yr = np.max(x)-np.min(x)
    yr = [np.min(x)-0.1*yr, np.max(x)+0.1*yr]
    #ax.plot([extent, extent], yr, 'r')
    #ax.plot([0, 0], yr, 'r')     
    #ax.plot([0, extent], [yr[0]]*2, 'r')
    #ax.plot([0, extent], [yr[1]]*2, 'r')
    ax2 = plt.subplot(122)
    ax2.set_title("PCA of Sliding Window Embedding")
    ax2.scatter(Y[:, 0], Y[:, 1])
    ax2.set_aspect('equal', 'datalim')
    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    T = 40 # The period in number of samples
    NPeriods = 4 # How many periods to go through
    N = T*NPeriods # The total number of samples
    t = np.linspace(0, 2*np.pi*NPeriods, N+1)[0:N] # Sampling indices in time
    x = np.cos(t)  # The final signal
    Tau = 10
    dT = 1
    dim = 10
    X = get_sliding_window_1D(x, dim, Tau, dT)
    Y, eigs = embedding_transform_1D(X)
    plot_sliding_window_1D(x, X, Y, Tau, dT)

    with open("slide_test.csv", 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerows(X)
    #   run persistence using Ripser
    rips = Ripser()
    rips.ComputeBarcode("slide_test.csv", 2, 10, 1, "point-cloud",1)
    barcode = rips.getBarcode()
    plot_persistence_diagram(barcode)
    '''
    grid = [[1,1,1,1,1],[1,0,0,0,1],[1,0,0,0,1],[1,0,0,0,1],[1,1,1,1,1]]
    with open("square.csv", 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerows(grid)
    #   try 2D von neumann filter
    filt = Filter2D()
    filt.loadBinaryFromFile("square.csv")
    filt.filterBinaryVonNeumann(10)
    filt.saveBinaryFiltration("square2.csv")
    cube2D = CubicalRipser2D()
    convert_csv_to_dipha("square2.csv", "square_dipha.csv")
    cube2D.ComputeBarcode("square_dipha.csv", "test.csv", "DIPHA", "LINKFIND", 10, True)
    barcode = cube2D.getBarcode()
    plot_persistence_diagram(barcode)

    #   Now try with RIPSER directly
    rips = Ripser()
    # save as point cloud format
    save_binary_cells_to_point_cloud(grid, "square_cloud.csv")

    # Run ripser on this
    rips.ComputeBarcode("square_cloud.csv", 2, 10, 1, "point-cloud", 1)

    # Plot the barcode
    barcode2 = rips.getBarcode()
    plot_persistence_diagram(barcode2)

    code1 = [[barcode[i][1],barcode[i][2]] for i in range(len(barcode))]
    code2 = [[barcode2[i][1],barcode2[i][2]] for i in range(len(barcode))]
    #   prepare bottleneck distance files
    with open("code1.txt", 'w') as csvfile:
        writer = csv.writer(csvfile,delimiter="\t")
        writer.writerows(code1)

    with open("code2.txt", 'w') as csvfile:
        writer = csv.writer(csvfile,delimiter='\t')
        writer.writerows(code2)

    #   bottleneck distance test
    bottle = BottleneckDistance()
    distance = bottle.Distance("code1.txt", "code2.txt", 10)
    print(distance)
    '''