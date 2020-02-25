import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.append("/home/czw/.local/lib/python3.6/site-packages/")
from astrowidgets import ImageWidget

#
#
# Helper to make reading saved data format agnostic.

def pq2df(filename, sep=None):
    try:
        with open(filename, 'rb') as f:
            out = pd.read_parquet(f)
    except:
        if sep is None:
            out = pd.read_csv(filename)
        else:
            out = pd.read_csv(filename, sep=sep)
    return out

#
#
# Begin plotting section

import matplotlib.colors as mcolors

def strCheck(val):
    if val is None or isinstance(val, str):
        val = (val, )
    return val

def makePlot(dataframe, xLabels, yLabels,
             loglike=None, diff=False, yf=1.0,
             xlim=None, ylim=None,
             fig=None, ax=None, colorCol=None
            ):
    doShow = True
    if fig is None:
        fig = plt.figure()
        doShow = False
    if ax is None:
        ax = plt.gca()
        doShow = False
    ax.set_alpha(0.9)
    xLabels = strCheck(xLabels)
    yLabels = strCheck(yLabels)

    xAxisLabel = []
    yAxisLabel = []
    for xL, yL in zip(xLabels, yLabels):
        XX = dataframe[xL]
        xAxisLabel.append(xL)

        if diff:
            YY = dataframe[xL] - yf * dataframe[yL]
            yAxisLabel.append(f"{xL} - {yL}")
        else:
            YY = yf * dataframe[yL]
            yAxisLabel.append(yL)
        if colorCol is None:
            ax.scatter(XX[~np.isnan(XX)], YY[~np.isnan(YY)], marker='.')
        else:
            colorData = dataframe[colorCol][~np.isnan(XX)].to_numpy(dtype=int) # (dtype=int)/112.0
            print("ScaleRange:", np.min(colorData), np.max(colorData))
            colorData = (colorData - np.min(colorData)) / (np.max(colorData) - np.min(colorData))

            imm = ax.scatter(XX[~np.isnan(XX)], YY[~np.isnan(YY)], marker='.',
                             c=colorData)
            fig.colorbar(imm, ax=ax)

    if loglike is not None:
        if 'x' in loglike:
            ax.set_xscale('log')
        if 'y' in loglike:
            ax.set_yscale('log')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_xlabel(", ".join(xAxisLabel))
    ax.set_ylabel(", ".join(yAxisLabel))
    ax.grid(True)
    if doShow:
        fig.show()

def makeHistPlot(dataframe, xLabels, yLabels=None,
                 loglike=None, yf=1.0,
                 Nbins=100,
                 xlim=None, ylim=None,
                 fig=None, ax=None

                ):
    doShow = True
    if fig is None:
        fig = plt.figure()
        doShow = False
    if ax is None:
        ax = plt.gca()
        doShow = False
    xLabels = strCheck(xLabels)
    yLabels = strCheck(yLabels)

    xAxisLabel = []
    yAxisLabel = []
    for xL, yL in zip(xLabels, yLabels):
        if yL is not None:
            XX = dataframe[xL] - yf * dataframe[yL]
            xAxisLabel.append(f"{xL} - {yL}")
        else:
            XX = dataframe[xL]
            xAxisLabel.append(xL)
        ax.hist(XX[~np.isnan(XX)], bins=Nbins, range=xlim)

    if loglike is not None:
        if 'x' in loglike:
            ax.set_xscale('log')
        if 'y' in loglike:
            ax.set_yscale('log')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_xlabel(", ".join(xAxisLabel))
    ax.set_ylabel("Counts")
    ax.grid(True)

    if doShow:
        fig.show()

def makeHist2DPlot(dataframe, xLabels, yLabels, zLabels=None,
                   loglike=None, diff=False,
                   yf=1.0,
                   Nbins=100,
                   range=None,
                   fig=None, ax=None

                ):
    doShow = True
    if fig is None:
        fig = plt.figure()
        doShow = False
    if ax is None:
        ax = plt.gca()
        doShow is False
    xLabels = strCheck(xLabels)
    yLabels = strCheck(yLabels)
    zLabels = strCheck(zLabels)
    xAxisLabel = []
    yAxisLabel = []
    for xL, yL, zL in zip(xLabels, yLabels, zLabels):
        XX = dataframe[xL]
        xAxisLabel.append(xL)
        if diff is True:
            if zL is None:
                YY = dataframe[xL] - yf * dataframe[yL]
                yAxisLabel.append(f"{xL} - {yL}")
            else:
                YY = dataframe[zL] - yf * dataframe[yL]
                yAxisLabel.append(f"{zL} - {yL}")
        else:
            YY = dataframe[yL]
            yAxisLabel.append(yL)
        ax.hist2d(XX[np.isfinite(XX)], YY[np.isfinite(XX)], range=range, bins=Nbins, norm=mcolors.LogNorm())
    if loglike is not None:
        if 'x' in loglike:
            ax.set_xscale('log')
        if 'y' in loglike:
            ax.set_yscale('log')
    ax.set_xlabel(", ".join(xAxisLabel))
    ax.set_ylabel(", ".join(yAxisLabel))
    if doShow:
        fig.show()

#
#
# Begin range shortcuts because of lazy.

range0 = (-1, 1)
range1 = (-1e-1, 1e-1)
range3 = (-1e-3, 1e-3)
range4 = (-1e-4, 1e-4)
range5 = (-1e-5, 1e-5)



