#!/usr/bin/env python3
"""
Created on Mon Feb  5 08:02:16 2018

@author: Esther Steiner
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 13:37:11 2018

@author: Esther Steiner
"""
import os
import sys
import re
import math

# collect flagged command line arguments in a dictionary
# and positional arguments in a list

def getopts(argv):
    flag_opts = {}
    positional_args = []
    while argv:
        if argv[0][0] == '-':
            flag_opts[argv[0]] = argv[1]
            argv = argv[2:]
        else:
            positional_args.append(argv[0])
            argv=argv[1:]
    return (flag_opts, positional_args)

# input the xyz coordinates of the atoms from the given file
# 
# the atoms are numbered sequentially, in the order in which they 
# are listed in the file
#
# the function returns a list of records, with each record being
# a list of [atom, x, y, z]
#
def fileInputXYZ( filename, debug ) :
    headings = ['atom','x','y','z']
    records = []
    recordCount = 0
    fileOne = open(filename,'r')
    while True:
        line = fileOne.readline()
        if not line: break
        theData = line.strip().split()
        theData[1] = float(theData[1])
        theData[2] = float(theData[2])
        theData[3] = float(theData[3])

        #print( theData )
        recordCount += 1
        theData[0] = theData[0] + str(recordCount)
        records.append(theData)
    return (records)

#
# determine the 3-D distance between two points
#
def distance(point1, point2):
    return math.sqrt((point2[0]-point1[0])**2 + (point2[1]-point1[1])**2
                 + (point2[2]-point1[2])**2)

#
# atomID is a string with the atom type, such as C or H, and the 
# atom number (assigned sequentially from the input file)
#
# this function returns the atom type from that string
#
def getAtomType(atomID):
    atom = ''
    for ch in atomID:
        if ch.isalpha():
            atom += ch
    return atom

#
# atomID is a string with the atom type, such as C or H, and the 
# atom number (assigned sequentially from the input file)
#
# this function returns the atom number from that string as an integer
#
def getAtomNumber( atomID ):
    atom = ''
    for ch in atomID:
        if ch.isnumeric():
            atom += ch
    return int(atom)

# 
# coordinates is a list of atoms with x, y, z coordinates
#
# threshhold specifies the maximum acceptable bond length
# distances greater than threshhold will not be considered
# as a bond
#
# the function returns a list of bonds.  Each item in the list
# is a list [atom1, atom2, distance]
#
def getBonds( coordinates, threshhold ):
    bondList = []
    atom1 = 0
    while atom1 < len(coordinates):
        atom2 = atom1 + 1
        while atom2 < len(coordinates):
            dist = distance(coordinates[atom1][1:], coordinates[atom2][1:])
            if dist <= threshhold and getAtomType(coordinates[atom1][0]) != 'H' \
                                  and getAtomType(coordinates[atom2][0]) != 'H':
                bondList.append([coordinates[atom1][0],coordinates[atom2][0],dist])
            atom2 += 1
        atom1 += 1
    return bondList

#
# bondLists is a two-dimensional list
# 
# the first subscript is the molecule number, with [0] containing the reference molecule
# the second subscript is the bond number
#
# this function returns a two-dimensional array with the following format
#
#   bond    bond lengths                             deltas
#   name    reference  other1  other2  other3 ...    other1-ref  other2-ref other3-ref ...
#

def calculate_deltas( bondLists ):
    bondNumber = 0
    delta_table = []
    while bondNumber < len(bondLists[0]):

       # the first set of bonds is used as the reference molecule

       referenceBond = bondLists[0][bondNumber][:]

       # put the reference bond information in the row

       row = []
       for item in referenceBond:
           row.append (item)

       # print('the reference bond is',row)
       

       # now we compare the bond in the reference to the same 
       # bond in all of the other molecules and append the bond lengths 
       # for each molecule to the rows in the table

       i = 1
       while i < len(bondLists):
           if (bondLists[i][bondNumber][0] != referenceBond[0] or \
                   bondLists[i][bondNumber][1] != referenceBond[1]):
               # the atoms in the bonds aren't the same
               # terminate the program with error
               print('The bonds in the reference molecule are not the same as ',
                      '\nthe bonds in file #', (i+1))
               print('Reference:', referenceBond, 'Other:',bondLists[i][bondNumber])
               sys.exit()

           row.append(bondLists[i][bondNumber][2])
           i += 1
       # end of bondLists loop


       #
       # last, we go through the row and calculate the columns of deltas
       # other1-ref, other2-ref, etc. and append the deltas to the row
       #
       # here is the format of each row before adding deltas
       #
       # item [0]   item [1]   item [2]          item [3] ...
       # atom1      atom2      ref_bond_length   other1_bond_length

       i = 0
       numOthers = len(row) - 3

       while i < numOthers:
           # print('row #', bondNumber, 'contains',row)
           row.append((row[3 + i] - row[2]))
           i += 1

       bondNumber += 1
       delta_table.append(row)
    
    # end of bondNumber loop

    # the delta table is complete

    print('\nDelta Table')
    print('Atom_1 Atom_2')
    print('Files processed: ', end = '')
    for name in fileList:
        print(name,end=' ')
    print( ) 
    for abond in delta_table:
        print(' '.join(map(str,abond)))
    return delta_table
#
# main
#os.chdir('G:/python programs')
from sys import argv
myargs = getopts(argv[1:])
opts = myargs[0]
fileList = myargs[1]
if (len(opts) > 0):
    if '-b' in opts:
        # print('base filename specified = ', opts['-b'])
        baseName = opts['-b']
        fileList = [baseName+'-N-opt.xyz', baseName+'-A2-opt.xyz',baseName+'-B2-opt.xyz']
elif (len(opts) == 0 and len(fileList) >= 2 ): 
    # print('the files to use are:', fileList)
    baseName = fileList[0][0:-4]
    # print('base filename will be ', baseName)
else:
    print('''Improper usage.  You have two choices: \
            calculateDeltas.py -b basename
            or
            calculateDeltas.py filename-neutral filename-cr0 filename-cr1
            ''')
    sys.exit()
bondLengths = []
fileData = []

#calculate the bond lengths for each file

for fileNumber in range( len(fileList) ):
    print('processing file:  ', fileList[fileNumber])
    fileData.append(fileInputXYZ(fileList[fileNumber], False))

    # distances above 1.9 angstroms will be considered non-bonding

    bondLengths.append( getBonds(fileData[fileNumber], 1.9) )
    # print('Data read from: ',fileList[fileNumber])
    # for atom in fileData[fileNumber]:
    #    print (atom)
    # print("")
    # print('Bond length table:')
    # for bond in bondLengths[fileNumber]:
    #     print(bond)

#create a summary table with delta values

deltaTable = calculate_deltas(bondLengths)


