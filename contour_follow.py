#!/usr/bin/env python

"""A program to take a binmap and create polygon regions from it."""

import math
import sys
from rdp import rdp

import fitsgz
import numpy as N

def areaPolygon(pts):
    """Get area of polygon.

    Returns -ve number for clockwise polygon
    Returns +ve number for anticlockwise polygon
    """

    tot = 0.
    N = len(pts)
    for i in xrange(N):
        j = (i+1) % N
        tot += pts[i][0]*pts[j][1] - pts[j][0]*pts[i][1]
    return tot * 0.5

def getLineSegments(image):
    """Finds all the (pixel-sized) line segments which surround bins.

    Returns a dict of segments for each bin.
    """

    print "Finding all line segments"
    allout = {}

    yw, xw = image.shape

    for y in xrange(yw):
        for x in xrange(xw):
            val = image[y, x]
            if val not in allout:
                allout[val] = []
            out = allout[val]

            # add line representing edge of pixel
            # these go clockwise (so they join up)
            if y+1 == yw or image[y+1, x] != val:
                out.append( (x, y+1, x+1, y+1) )
            if y == 0 or image[y-1, x] != val:
                out.append( (x+1, y, x, y) )
            if x+1 == xw or image[y, x+1] != val:
                out.append( (x+1, y+1, x+1, y) )
            if x == 0 or image[y, x-1] != val:
                out.append( (x, y, x, y+1) )

    return allout

def mergeStraightSegments(lines):
    """Remove extra unneeded points from polygon."""

    for line in lines:

        i = 0
        while i < len(line)-2:
            # check pairs of lines to see whether they have the same
            # angle
            angle1 = math.atan2( line[i][0] - line[i+1][0],
                                 line[i][1] - line[i+1][1] )
            angle2 = math.atan2( line[i+1][0] - line[i+2][0],
                                 line[i+1][1] - line[i+2][1] )

            if abs(angle1-angle2) < 1e-4:
                # remove unneeded point
                del line[i+1]
            else:
                i += 1

def joinLineSegments(segments):
    """Join line segments into a line."""

    # copy input
    insegs = list(segments)

    outlines = []

    line = [ insegs.pop() ]

    while len(insegs) != 0:

        found = False
        for i, seg in enumerate(insegs):
            if seg[0] == line[-1][2] and seg[1] == line[-1][3]:
                insegs.pop(i)
                line.append(seg)
                found = True
                break
        if not found:
            outlines.append(line)
            line = [ insegs.pop() ]

    if line:
        outlines.append(line)
    return outlines

if __name__ == '__main__':

    image = fitsgz.open(sys.argv[1])[0].data

    allsegs = getLineSegments(image.astype(N.int32))

    for bin, segs in allsegs.iteritems():
        if bin < 0:
            continue

        print "Joining bin", bin
        lines = joinLineSegments(segs)
        mergeStraightSegments(lines)

        newlines = []
        for l in lines:
            newlines.append( rdp(l, 4.0) )
        lines = newlines
        #t.simplifyAreaPolygon(lines)

        include = []
        exclude = []
        for line in lines:
            area = areaPolygon(line)
            #pts = [("%g,%g" % (p[0]*25+0.5, p[1]*25+0.5)) for p in line]
            pts = [("%g,%g" % (p[0]+0.5, p[1]+0.5)) for p in line]
            txt = "polygon(%s)" % ','.join(pts)

            if area < 0.:
                # clockwise, so include
                include.append(txt)
            else:
                exclude.append(txt)

        fout = open('xaf_%i.reg' % bin, 'w')
        #fout.write('image\n')
        for x in include:
            print >>fout, x
        for x in exclude:
            print >>fout, "-%s" % x
        fout.close()

