# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 09:53:35 2015

@author: dylanbarth
"""

##############################
#####   Virus Model     ######
##############################

####################################
# Defining variables and functions #
####################################
import itertools
import numpy
import math

Nx = 10    # Range of cells on x axis
Ny = 1      # Range of cells on y axis
ts = 1e-01     # Time step for model
TauI = 12   # Avg time for infection
TauE = 6    # Avg time for eclipse
ne = 60     # 
ni = 60     #
probi = .5  # Threshold for probablility of cell to cell infection
nmean = 0.5 # Mean of probability of cell to cell infection

#############

rho = 224
beta = 3e-06
c = 0.056
deltx = 50e-09

def indx(condition, axis):                  # Translate array index to int index           
    if axis == 'yaxis':
        axis = 0
    if axis == 'xaxis':
        axis = 1
    if axis == 'zaxis':
        axis = 2                            #This function allows us to get an integer
    loc = numpy.where(condition)            #index number out of a numpy array (vector)
    lloc = list(str(loc[axis:axis +1]))     # where you input 'xaxis', 'yaxis', or 'zaxis'
    list1 = []
    list2 = []
    for i in lloc:
        f = i.isdigit()
        if f is True:
            list1.append(int(i))
    n = len(list1)
    for i in list1:
        m = list1.index(i)
        list2.append(i*(10**(n-m-1)))
        list1[m:m+1] = [0]
    locint = sum(list2)
    return loc
    return locint

def Te():                                   #Picks a random number from gamma
    return numpy.random.gamma(1/math.sqrt(ne), TauE/math.sqrt(ne))
    
def Ti():                                   #Picks a  random number from gamma
    return numpy.random.gamma(1/math.sqrt(ni), TauI/math.sqrt(ni))
    
def P1():
    return numpy.random.normal(0.5, 0.2)    #Picks a random number from gaussian

def P2():
    return numpy.random.normal(0.5, 0.2)    #Picks a random number from gaussian
    

###########################################################
#  Produces a 1x100 matrix of healthy cells (h)           #
#  Produces a time matrix for after eclipse phase (e)     #
#  Produces a time matrix for after infection phase (i)   #
#  Produces a time matrix hor healthy cells (t)           #
#  Produces a univeral time matrix (uv)                   #
###########################################################

clist = ['h']*Nx
cells = numpy.matrix(clist)
ecl = numpy.zeros((Ny,Nx))
inf = numpy.zeros((Ny,Nx))
th = numpy.zeros((Ny,Nx))
ut = numpy.zeros((Ny,Nx))

############################################################
#   Straight up not sure what I'm doing here
#   here's to hoping it works!
#                                       _______
#                              ________/_____
#                              \     /     |
#                               \~~~~~~~~~|
#                                \  `    |
#                                 \    `|  
#                                  \ ` |
#                                   \ |
#                                    T   
#                                  -----   
############################################################
for j in cells == 'i':
    Itn = 1
    Vtn = 1

Vt1n = rho*()


################################################################
#   Infects a random cell, now seen as (i)                     #
################################################################

rand = numpy.random.randint(1,Nx)
cells[0,rand] = 'i'
inf[0,rand] = Ti()


###################################################################
#                                                                 #
#                        Runs simulation                          #
#                                                                 #
###################################################################

step = itertools.cycle('ABCDE')
while (cells.all != 'd') is True:
#####################################
#       The Universal Time          #
#       is kept here (ut)           #
#####################################
    next(step)
    if next(step) == 'A':
        ut = ut + numpy.matrix([ts]*Nx)
        print ut
        print 'ut'
#####################################
#       The Healthy Cells' time     #
#       is kept here (th)           #
#####################################
    if next(step) == 'B':
        hindex = numpy.where(cells == 'h')[1]    
        for j in hindex:      
            th[0,j] = th[0,j] + ts
            print th
            print 'th'
#####################################
#    Eclipse phase -> Infection     #
#                                   #
#####################################
    if next(step) == 'C':
        eindex = numpy.where(cells == 'e')[1]
        print (ecl + th)    
        print 'ecl + th'
        print ut
        print 'ut2'
        for j in eindex:
            if (ecl[0,j] + th[0,j]) < ut[0,j]:
                cells[0,j] = 'i'
                inf[0,j] = inf[0,j] + Ti()
                print "----------------IT WORKED--------------"
            else:
                print 'nope'
#####################################
#       Infection spreads &         #
#       Infectious cells die        #
#####################################
    if next(step) == 'D':
        iindex = numpy.where(cells == 'i')[1]
        for j in iindex:
            if j == 'd':
                break
            if P2()>0.5:
                if j == 0:
                    break
                if cells[0,(j -1)] == 'h':
                    cells[0,(j -1)] = 'e'
                    ecl[0,(j -1)] = ecl[0, j -1] + Te()
            if P1()>0.5:
                if j == (Nx-1):
                    break
                if cells[0,(j +1)] == 'h':
                    cells[0,(j +1)] = 'e'
                    ecl[0,(j +1)] = ecl[0,(j +1)] + Te()
            if ut[0,j] > (inf[0,j] + ecl[0,j] + th[0,j]):
                cells[0,j] = 'd'       
        print inf
        print 'inf'
#####################################
#       Prints status of cells      #
#                                   #
#####################################
    if next(step) == 'E':
        print cells
########################################
#       If the infection can't         #
#       spread, this breaks the loop   #
########################################
    #if cells.all != ('h') or ('d'):           #ASK DR D!
    #    break
 

