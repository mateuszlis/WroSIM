#!/usr/bin/python
import sys

if len( sys.argv ) < 3:
    print "Usage: calcClusterSize.py clusters.dat clusterSize.dat"

result = []
frameNums = []
lastFrame = ""
count = 0
size = 0
with open( sys.argv[1], 'r' ) as f:
    for line in f.readlines():
        lineList = line.split("\t")
        frameNum = lineList[0]
        if frameNum != lastFrame:
            if count != 0:
                result.append( float( size ) / count )
                frameNums.append( lastFrame )
            lastFrame = frameNum
            count = 0
            size = 0
        size += int( lineList[1] ) * int( lineList[2] )
        count += int( lineList[2] )
    if count != 0:
        result.append( float( size ) / count )
        frameNums.append( lastFrame )
with open( sys.argv[2], 'w' ) as fWrite:
    for num in range( len( result ) ):
        fWrite.write( frameNums[num] + "\t" + str( result[num] ) + "\n" )




