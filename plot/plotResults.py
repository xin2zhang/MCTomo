##!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os.path

def readinp(filename):
    # read input file
    with open(filename,'r') as f:
        data = f.readlines()
    runInfo = {}
    for ln in data:
        if not '=' in ln:
            continue
        varname, val = ln.split('=')
        mname, key = varname.split('%')
        if("'" in val):
            runInfo[key.strip()] = val.strip().strip("'")
        elif('.' in val):
            runInfo[key.strip()] = float(val.strip().strip(','))
        else:
            runInfo[key.strip()] = int(val.strip().strip(','))
    return runInfo

def getNumberOfSamples(filename):
    # open txt file
    with open(filename,'r') as f:
        data = f.readlines()
    nsamples = int(data[2].split()[0])
    return nsamples

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot mean and stdev of the tomography results")
    parser.add_argument("base", metavar='base', nargs=1, help='base folder of results')
    parser.add_argument("-n", "--nchains", metavar='nchains', default=None, type=int, help='The number of chains')
    parser.add_argument("-v", "--vtype", metavar='vtype', default='vp', help='plot velocity type')
    parser.add_argument("-x", "--xslice", metavar='xslice', default=None, type=float, help='slice location at x direction')
    parser.add_argument("-y", "--yslice", metavar='yslice', default=None, type=float, help='slice location at y direction')
    parser.add_argument("-z", "--zslice", metavar='zslice', default=None, type=float, help='slice location at z direction')
    parser.add_argument("-c", "--cmap", metavar='cmap', default='seismic_r', help='colormap')

    args = parser.parse_args()
    base = args.base[0]
    nchains = args.nchains
    vtype = args.vtype
    xslice = args.xslice
    yslice = args.yslice
    zslice = args.zslice
    cmap = args.cmap

    # read input file to get some necessary information
    inp = os.path.join(base,'MCTomo.inp')
    if(os.path.isfile(inp)):
        runInfo = readinp(inp)
    else:
        print('input file does not exist in base folder')
        exit(1)
    nx = runInfo['NX']
    ny = runInfo['NY']
    nz = runInfo['NZ']
    if(nchains is None):
        nchains = runInfo['NUMBER_OF_TEMPERATURES']

    # recs
    recs = np.loadtxt(os.path.join(base,'receivers.dat'),skiprows=1)

    # read run time info file to get correct number of samples for each chain
    nsamples = np.zeros(nchains)
    mean = np.zeros((nchains,2*ny*nx,nz))
    var = np.zeros((nchains,2*ny*nx,nz))
    if(os.path.exists(os.path.join(base,'Results'))):
        for i in range(1,nchains+1):
            runInfoFile = 'Results/resume/run_time_info_'+str(i)+'.dat'
            nsamples[i-1] = getNumberOfSamples(os.path.join(base,runInfoFile))
            runMeanFile = 'Results/resume/run_time_average_'+str(i)+'.dat'
            runVarFile = 'Results/resume/run_time_var_'+str(i)+'.dat'
            cmean = np.loadtxt(os.path.join(base,runMeanFile), skiprows=1)
            cvar = np.loadtxt(os.path.join(base,runVarFile), skiprows=1)
            if(cmean.shape[0]==ny*nx):
                mean[i-1,0:ny*nx,:] = cmean
                var[i-1,0:ny*nx,:] = cvar
            elif(cmean.shape[0] == nx):
                mean[i-1,0:ny*nx,0] = np.reshape(cmean, ny*nx)
                var[i-1,0:ny*nx,0] = np.reshape(cvar, ny*nx)
                nz = 1
            elif(cmean.shape[0] == nx*2):
                mean[i-1,:,0] = np.reshape(cmean, 2*ny*nx)
                var[i-1,:,0] = np.reshape(cvar, 2*ny*nx)
                nz = 1
            else:
                mean[i-1,:,:] = cmean
                var[i-1,:,:] = cvar
    elif(os.path.exists(os.path.join(base,str(1)))):
        for i in range(1,nchains+1):
            runInfoFile = 'Results/resume/run_time_info_'+str(i)+'.dat'
            nsamples[i-1] = getNumberOfSamples(os.path.join(base,str(i),runInfoFile))
            runMeanFile = 'Results/resume/run_time_average_'+str(i)+'.dat'
            runVarFile = 'Results/resume/run_time_var_'+str(i)+'.dat'
            cmean = np.loadtxt(os.path.join(base,str(i),runMeanFile), skiprows=1)
            cvar = np.loadtxt(os.path.join(base,str(i),runVarFile), skiprows=1)
            if(cmean.shape[0]==ny*nx):
                mean[i-1,0:ny*nx,:] = cmean
                var[i-1,0:ny*nx,:] = cvar
            elif(cmean.shape[0] == nx):
                mean[i-1,0:ny*nx,0] = np.reshape(cmean, ny*nx)
                var[i-1,0:ny*nx,0] = np.reshape(cvar, ny*nx)
                nz = 1
            elif(cmean.shape[0] == nx*2):
                mean[i-1,:,0] = np.reshape(cmean, 2*ny*nx)
                var[i-1,:,0] = np.reshape(cvar, 2*ny*nx)
                nz = 1
            else:
                mean[i-1,:,:] = cmean
                var[i-1,:,:] = cvar
    else:
        print('Cannot find results files')
        exit(2)

    # get correct mean and var for all chains
    mean = np.sum(mean*nsamples[:,None,None], axis=0)/np.sum(nsamples)
    var = np.sum(var*nsamples[:,None,None], axis=0)/np.sum(nsamples)
    std = np.sqrt(var-mean**2)

    # plot slice section
    if(xslice is None):
        xslice = (runInfo['XMIN']+runInfo['XMAX'])/2
    if(yslice is None):
        yslice = (runInfo['YMIN']+runInfo['YMAX'])/2
    if(zslice is None):
        zslice = (runInfo['ZMIN']+runInfo['ZMAX'])/2

    # horizontal slice
    if(nz>1):
        dz = (runInfo['ZMAX'] - runInfo['ZMIN'])/(runInfo['NZ']-1)
    else:
        dz = runInfo['ZMAX'] - runInfo['ZMIN']
    indz = int( (zslice-runInfo['ZMIN'])/dz )
    fig, ax = plt.subplots(2,1)
    if(vtype == 'vp'):
        im = ax[0].imshow(np.reshape(mean[0:ny*nx,:],(nx,ny,nz))[:,:,indz].T, extent=(runInfo['XMIN'], runInfo['XMAX'], runInfo['YMIN'], runInfo['YMAX']),cmap=cmap)
        fig.colorbar(im, ax=ax[0])
        ax[0].scatter(recs[:,0],recs[:,1],marker='^')
        im = ax[1].imshow(np.reshape(std[0:ny*nx,:],(nx,ny,nz))[:,:,indz].T, extent=(runInfo['XMIN'], runInfo['XMAX'], runInfo['YMIN'], runInfo['YMAX']),cmap=cmap)
        fig.colorbar(im, ax=ax[1])
        ax[1].scatter(recs[:,0],recs[:,1],marker='^')
    else:
        im = ax[0].imshow(np.reshape(mean[ny*nx:2*ny*nx,:],(nx,ny,nz))[:,:,indz].T, extent=(runInfo['XMIN'], runInfo['XMAX'], runInfo['YMIN'], runInfo['YMAX']),cmap=cmap)
        fig.colorbar(im, ax=ax[0])
        im = ax[1].imshow(np.reshape(std[ny*nx:2*ny*nx,:],(nx,ny,nz))[:,:,indz].T, extent=(runInfo['XMIN'], runInfo['XMAX'], runInfo['YMIN'], runInfo['YMAX']),cmap=cmap)
        fig.colorbar(im, ax=ax[1])
    plt.show()

    if(nz>1):
        # xslice
        dx = (runInfo['XMAX'] - runInfo['XMIN'])/(runInfo['NX']-1)
        indx = int( (xslice-runInfo['XMIN'])/dx )
        fig, ax = plt.subplots(2,1)
        if(vtype == 'vp'):
            ax[0].imshow(np.reshape(mean[0:ny*nx,:],(nx,ny,nz))[indx,:,:].T, extent=(runInfo['YMIN'], runInfo['YMAX'], runInfo['ZMIN'], runInfo['ZMAX']),cmap=cmap)
            ax[1].imshow(np.reshape(std[0:ny*nx,:],(nx,ny,nz))[indx,:,:].T, extent=(runInfo['YMIN'], runInfo['YMAX'], runInfo['ZMIN'], runInfo['ZMAX']),cmap=cmap)
        else:
            ax[0].imshow(np.reshape(mean[ny*nx:2*ny*nx,:],(nx,ny,nz))[indx,:,:].T, extent=(runInfo['YMIN'], runInfo['YMAX'], runInfo['ZMIN'], runInfo['ZMAX']),cmap=cmap)
            ax[1].imshow(np.reshape(std[ny*nx:2*ny*nx,:],(nx,ny,nz))[indx,:,:].T, extent=(runInfo['YMIN'], runInfo['YMAX'], runInfo['ZMIN'], runInfo['ZMAX']),cmap=cmap)

        # yslice
        dy = (runInfo['YMAX'] - runInfo['YMIN'])/(runInfo['NY']-1)
        indy = int( (yslice-runInfo['YMIN'])/dy )
        fig, ax = plt.subplots(2,1)
        if(vtype == 'vp'):
            ax[0].imshow(np.reshape(mean[0:ny*nx,:],(nx,ny,nz))[:,indy,:].T, extent=(runInfo['XMIN'], runInfo['XMAX'], runInfo['ZMIN'], runInfo['ZMAX']))
            ax[1].imshow(np.reshape(std[0:ny*nx,:],(nx,ny,nz))[:,indy,:].T, extent=(runInfo['XMIN'], runInfo['XMAX'], runInfo['ZMIN'], runInfo['ZMAX']))
        else:
            ax[0].imshow(np.reshape(mean[ny*nx:2*ny*nx,:],(nx,ny,nz))[:,indy,:].T, extent=(runInfo['XMIN'], runInfo['XMAX'], runInfo['ZMIN'], runInfo['ZMAX']))
            ax[1].imshow(np.reshape(std[ny*nx:2*ny*nx,:],(nx,ny,nz))[:,indy,:].T, extent=(runInfo['XMIN'], runInfo['XMAX'], runInfo['ZMIN'], runInfo['ZMAX']))
