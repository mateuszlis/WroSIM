#!/usr/bin/python

from sys import argv
if len( argv ) > 5 or len( argv ) < 5:
    print "Not enough arguments usage: range_sum.py file.dat binSize startTime endTime"
    exit
else:
    multiplier = int(argv[2]) 
    d = [ 0 for i in range( 5000 ) ]
    startTime = int( argv[3] )
    endTime = int( argv[4] )
    
    with open( argv[1], 'r' ) as f:
        lines = f.readlines()
        for line in lines:
            values = line.split()
            frame, key, value = int( values[0] ), int( values[1] ), int( values[2] )
            if frame < endTime and frame > startTime: 
                d[ key / multiplier ] += value

    with open( "range_summed.dat", "w" ) as f:
        for i in range( len( d ) ):
            if d[i] > 0:
                f.write( str( i * multiplier ) + " " + str( d[i] * i * multiplier ) + "\n" )


        

        
