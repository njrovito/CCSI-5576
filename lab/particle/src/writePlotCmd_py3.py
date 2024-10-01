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
    
    g = open('pc','w')
    g.write("reset\n")
    g.write("set xrange [-.1:1.1]\n")
    g.write("set yrange [-.1:1.1]\n")
        
    for time in sorted(timeList):
        g.write( 'plot \\\n')
        for PE in sorted(PElist):
            g.write( '"mesh_' + str(PE) + '_' + str(time) + '.plt" w l , "ptcl_'+str(PE)+'_' + str(time) + '.plt" w p lw 3, \\\n')
        g.write( '\n')
        g.write( 'pause .01\n')
    g.close()

    # Write mesh only command
    
    g = open('pc_mesh','w')
    g.write( 'reset\n')
    g.write( 'set xrange [-.1:1.1]\n')
    g.write( 'set yrange [-.1:1.1]\n')
        
    for time in sorted(timeList):
        g.write( 'plot \\\n')
        for PE in sorted(PElist):
            g.write( '"mesh_' + str(PE) + '_' + str(time) + '.plt" w l , \\\n')
        g.write( '\n')
        g.write( 'pause .1\n')
    g.close()


if __name__ == "__main__":
    writePlotCmd(sys.argv[1:])
