#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import KDTree

def plot_2d(sites,values):
    # plot 2d voronoi diagrams
    ptmin = np.amin(sites,axis=0)
    ptmax = np.amax(sites,axis=0)
    ptmin = ptmin - 0.1*(ptmax-ptmin)
    ptmax = ptmax + 0.1*(ptmax-ptmin)

    # create a KDTree
    tree = KDTree(sites)

    # create query points
    x = np.linspace(ptmin[0],ptmax[0],101)
    y = np.linspace(ptmin[1],ptmax[1],101)
    grid = np.meshgrid(x,y)
    positions = np.vstack(map(np.ravel,grid))

    # create voronoi structure using kdtree
    print(str(positions.shape))
    ind = tree.query(positions.T,return_distance=False)
    gvalues = np.reshape(values[ind], (len(y),len(x)))

    # plot figures
    plt.figure(figsize=(500,500/(x[-1]-x[0])*(y[-1]-y[0])))
    plt.imshow(gvalues, extent=(ptmin[0], ptmax[0], ptmin[1], ptmax[1]))
    plt.xlabel('X(km)')
    plt.ylabel('Y(km)')
    plt.colorbar()
    plt.show()

def plot_3d(sites, values, xslice=None, yslice=None, zslice=None):
    # plot 2d voronoi diagrams
    ptmin = np.amin(sites,axis=0)
    ptmax = np.amax(sites,axis=0)
    ptmin = ptmin - 0.05*(ptmax-ptmin)
    ptmax = ptmax + 0.05*(ptmax-ptmin)

    if(xslice is None):
        xslice = (ptmin[0]+ptmax[0])/2
    if(yslice is None):
        yslice = (ptmin[1]+ptmax[1])/2
    if(zslice is None):
        zslice = (ptmin[2]+ptmax[2])/2

    # create a KDTree
    tree = KDTree(sites)

    # create query points
    x = np.linspace(ptmin[0],ptmax[0],101)
    y = np.linspace(ptmin[1],ptmax[1],101)
    z = np.linspace(ptmin[2],ptmax[2],101)
    grid = np.meshgrid(x,y,z)
    positions = np.vstack(map(np.ravel,grid))

    # create voronoi structure using kdtree
    print(str(positions.shape))
    ind = tree.query(positions.T,return_distance=False)
    gvalues = np.reshape(values[ind], (len(y),len(x),len(z)))

    # plot figures
    xind = int((xslice-ptmin[0])/(x[1]-x[0]))
    yind = int((yslice-ptmin[1])/(y[1]-y[0]))
    zind = int((zslice-ptmin[2])/(z[1]-z[0]))
    # fig, ax = plt.subplots(3, 1, figsize=(200,600))
    # ax[0].imshow(gvalues[:,:,zind], extent=(ptmin[0], ptmax[0], ptmin[1], ptmax[1]))
    # ax[0].set_xlabel('X(km)')
    # ax[0].set_ylabel('Y(km)')
    # ax[1].imshow(gvalues[yind,:,:], extent=(ptmin[1], ptmax[1], ptmin[2], ptmax[2]))
    # ax[2].imshow(gvalues[:,xind,:], extent=(ptmin[0], ptmax[0], ptmin[2], ptmax[2]))
    plt.figure(figsize=(400,400/(x[-1]-x[0])*(y[-1]-y[0])))
    plt.imshow(gvalues[:,:,zind], extent=(ptmin[0], ptmax[0], ptmin[1], ptmax[1]), aspect='auto')
    plt.xlabel('X(km)')
    plt.ylabel('Y(km)')
    plt.colorbar()
    plt.show()
    plt.figure(figsize=(500,400/(x[-1]-x[0])*(z[-1]-z[0])))
    plt.imshow(gvalues[yind,:,:], extent=(ptmin[0], ptmax[0], ptmin[2], ptmax[2]), aspect='auto')
    plt.gca().invert_yaxis()
    plt.xlabel('X(km)')
    plt.ylabel('Z(km)')
    plt.colorbar()
    plt.show()
    plt.figure(figsize=(500,400/(y[-1]-y[0])*(z[-1]-z[0])))
    plt.imshow(gvalues[:,xind,:], extent=(ptmin[1], ptmax[1], ptmin[2], ptmax[2]), aspect='auto')
    plt.gca().invert_yaxis()
    plt.xlabel('Y(km)')
    plt.ylabel('Z(km)')
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Voronoi cells")
    parser.add_argument("fname", metavar='fname', nargs=1, help='Voronoi sites filename')
    parser.add_argument("-v", "--vtype", metavar='vtype', default='vp', help='plot velocity type')
    parser.add_argument("-x", "--xslice", metavar='xslice', default=None, type=float, help='slice location at x direction')
    parser.add_argument("-y", "--yslice", metavar='yslice', default=None, type=float, help='slice location at y direction')
    parser.add_argument("-z", "--zslice", metavar='zslice', default=None, type=float, help='slice location at z direction')

    args = parser.parse_args()
    fname = args.fname[0]
    vtype = args.vtype
    xslice = args.xslice
    yslice = args.yslice
    zslice = args.zslice

    # load Voronoi sites and values
    voronoi = np.loadtxt(fname, skiprows=1)
    nsites = voronoi.shape[0]
    nd = voronoi.shape[1]
    if(nd == 5):
        sites = voronoi[:, 0:2]
        if(vtype=='vp'):
            values = voronoi[:, 2]
        elif(vtype=='vs'):
            values = voronoi[:, 3]
        else:
            values = voronoi[:, 4]
        plot_2d(sites,values)
    elif(nd == 6):
        sites = voronoi[:, 0:3]
        if(vtype=='vp'):
            values = voronoi[:, 3]
        elif(vtype=='vs'):
            values = voronoi[:, 4]
        else:
            values = voronoi[:, 5]
        plot_3d(sites,values, xslice, yslice, zslice)
    else:
        print("Not supported dimensionality")
