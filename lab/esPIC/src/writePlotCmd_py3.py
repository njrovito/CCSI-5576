#!/usr/bin/python3

import os

import sys
import getopt
import glob

# ==
# ||
# || Main Program
# ||
# ==

def writePlotCmd(argv):



    # List of all ptcl files
    
    ptclList = glob.glob("ptcl*.plt")

    # Parse the filenames counting the number of PEs for this run

    PElist = []
    timeList = []
    
    for i in range(0,len(ptclList)):
        name = ptclList[i];
        tmp = name.split('_')

        PE  = tmp[1]
        if PE not in PElist:
            PElist.append(PE)

        tmp2 = tmp[2].split(".")
        time = tmp2[0]
        if time not in timeList:
            timeList.append(int(time))


    # Write combined plot command
    
    g = open('pc_qval','w')
    print (  'set xrange [-.1:1.1]', file = g )
    print (  'set yrange [-.1:1.1]', file = g )
        
    for time in sorted(timeList):
        print ( 'splot \\', file = g)
        for PE in sorted(PElist):
            print ( '"qval_' + str(PE) + '_' + str(time) + '.plt" w l , "ptcl_'+str(PE)+'_' + str(time) + '.plt" w p lw 1, \\', file = g )
        print ( '', file = g )
        print ( 'pause .01', file = g )
    g.close()

    # Write particles only
    
    g = open('pc_voltage','w')
#    print (  , 'reset'
    print (  'set xrange [-.1:1.1]', file = g )
    print (  'set yrange [-.1:1.1]', file = g )
        
    for time in sorted(timeList):
        print (  'splot \\', file = g )
        for PE in sorted(PElist):
            print ( '"voltage_' + str(PE) + '_' + str(time) + '.plt" w l , "ptcl_'+str(PE)+'_' + str(time) + '.plt" w p lw 1, \\', file = g )
        print ( '', file = g )
        print ( 'pause .01', file = g )
    g.close()


if __name__ == "__main__":
    writePlotCmd(sys.argv[1:])
