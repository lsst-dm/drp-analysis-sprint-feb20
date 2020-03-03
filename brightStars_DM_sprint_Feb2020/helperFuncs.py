""" This script contains a few helper functions that are common to the
three notebooks for exploring detections around bright stars."""
import numpy as np
from lsst import geom
import time

def extractBrightStarInfo(tractInfo, mask_path, radMaxPix=0, radInPix=True, verbose=False):
    """ Extract mask information from the ds9 region files. This is pretty hacky, but is
    only used for these specific masks - ultimately, it would be better to select the bright
    objects we look at from an external catalog (Gaia)."""
    f = open(mask_path, 'r')

    # boundaries; shrink BBox so largest annuli can fit in entirely
    bBox = tractInfo.getBBox()
    bBox.grow(-int(radMaxPix)-1)
    wcs = tractInfo.getWcs()

    pix_centers = []
    mags = []
    radii = []

    for line in f:
        if line[:6] == 'circle':
            maskinfo, comm = line.split('#')
            objid, mag = comm.split(',')
            maskinfo = maskinfo.split(',')
            radius = maskinfo[2].split('d')[0] # keep everything before the d
            circle = [float(maskinfo[0][7:]), # remove "circle("
                        float(maskinfo[1])]
            sp = geom.SpherePoint(circle[0], circle[1], geom.degrees)
            cpix = wcs.skyToPixel(sp)
            if (cpix[0] >= bBox.beginX and
                cpix[0] < bBox.endX):
                if (cpix[1] >= bBox.beginY and
                    cpix[1] < bBox.endY):
                    # also save magnitude (as put down by Andy)...
                    mag = float(mag[5:-2]) # remove " mag:" and the trailing \n
                    mags += [mag]
                    # ... and mask radius
                    radii += [float(radius)]
                    # if the center of the mask is inside, keep it
                    pix_centers += [cpix]
        elif line[:3] == 'box': # ignore saturation spikes/bleed trails boxes
            pass
        else:
            if verbose:
                # check we didn't miss anything important
                print(line)
    f.close()
    mags = np.array(mags)
    # convert radii (in degrees) to arcsecpixels
    radii = np.array(radii) * 3600
    # and, if requested, to pixels
    if radInPix:
        radii /= geom.radToArcsec(wcs.getPixelScale())
    return pix_centers, mags, radii

def countInAnnulus(src, brightCen, radiiPix, x_name="base_SdssCentroid_x",
                   y_name="base_SdssCentroid_y"):
    radMaxPix = radiiPix[-1]
    # only consider detections within a square of width 2 outer radii
    withinSquare = ((src[x_name] >= brightCen[0]-radMaxPix) &
                    (src[x_name] < brightCen[0]+radMaxPix) &
                    (src[y_name] >= brightCen[1]-radMaxPix) &
                    (src[y_name] < brightCen[1]+radMaxPix))
    insideCents = np.array([obj.getCentroid() for obj in src[withinSquare]])
    # compute distance to bright object for all candidates
    dist = np.sqrt(np.sum((insideCents - brightCen)**2, axis=1))
    counts = []
    areas = []
    for radIn,radOut in zip(radiiPix,radiiPix[1:]):
        # count those within annulus
        inAnnulus = (dist >= radIn) & (dist < radOut)
        counts += [np.sum(inAnnulus)]
        areas += [np.pi *(radOut**2 - radIn**2)]
    return counts, areas

def countsBeyondMask(src, nbAnnuli, annSizePix, brightCenters, brightRadii, verbose=False):
    """Compute detected source counts around bright star masks."""
    allCounts, allAreas = [], []
    for j,brightCen in enumerate(brightCenters):
        start = time.time()
        # compute annuli radii for this bright object
        brad = brightRadii[j]
        #radiiPix = np.arange(brad, brad + annSizePix*(nbAnnuli + 1), annSizePix
        radiiPix = np.array([brad + k*annSizePix for k in range(nbAnnuli+1)])
        counts, areas = countInAnnulus(src, brightCen, radiiPix)
        end = time.time()
        if not j%500 and verbose:
            print(' > Bright object {} out of {}; time elapsed for this object: {}s'.format(
                j+1,len(brightCenters),end-start))
        allCounts += [counts]
        allAreas += [areas]
    allCounts = np.array(allCounts)
    allAreas = np.array(allAreas)
    return allCounts, allAreas

def boxSelector(centers, radii, x0, y0, xWidth, yWidth, extraBuffer=0,
                from_center=False):
    """ Box selector, to be used as a `furtherBrightObjectSelector` if we are
    not looking at whole tracts at a time.
    If `from_center`, check that centers fall within box; if `False`, that
    circles of radii `radii` do."""
    if from_center:
        buffers = np.zeros(radii.shape) + extraBuffer
    else:
        buffers = radii + extraBuffer
    inside = np.array([(cent[0] - buff > x0) &
                       (cent[0] + buff < x0 + xWidth) &
                       (cent[1] - buff > y0) &
                       (cent[1] + buff < y0 + yWidth)
                   for cent,buff in zip(centers, buffers)])
    return np.where(inside)[0]
