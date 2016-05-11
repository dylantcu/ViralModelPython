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
import sys
import numpy
import math
from matplotlib import pyplot as p

Nx = 8      # Range of cells on x axis
Ny = 8      # Range of cells on y axis
Nt = 100    # Days simulation runs for
ts = 1e-2   # Time step for model
TauI = 12   # Avg time for infection
TauE = 6    # Avg time for eclipse
ne = 60     # 
ni = 60     #
probi = .5  # Threshold for probablility of cell to cell infection
nmean = 0.5 # Mean of probability of cell to cell infection
scale = 10000 # scaling factor between k and ts
r = 0
#############

rho = 2e07
D = 1.1e-11
c = 0.056
deltx = 50e-06
deltxprime = deltx #*10e6
dumb = 0



def Te():                                   #Picks a random number from gamma
    return numpy.random.gamma(TauE, TauE/math.sqrt(ne))
    
def Ti():                                   #Picks a  random number from gamma
    return numpy.random.gamma(TauI, TauI/math.sqrt(ni))
    
def PG1():
    return numpy.random.normal(0.5, 0.2)    #Picks a random number from gaussian

def PU1():
    return numpy.random.uniform(low=0.0, high=1.0, size=None)

        
def deleteContent(fName):                   #Deletes past content of a text file
    with open(fName, "w"):
        pass
    
def diffVModel(vindexy):
    g = vindexy[:,j]
    vtemp[g[1],g[0],1] = vtemp[g[1],g[0],2]
    if f % scale == 0:
        vir[g[1],g[0],(f/scale)] = vtemp[g[1],g[0],1]
        print "divisible"
    
    if cells[g[1],g[0]] == 'i':                 #Is the cell infected?#
        rho2 = rho
    else:
        rho2 = 0

    if g[0] == 0:
        vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(-2*vtemp[g[1],g[0],1] + vtemp[g[1],(g[0]+1),1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]                                    
    elif g[0] == Nx-1:
       vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(vtemp[g[1],(g[0] -1),1] - 2*vtemp[g[1],g[0],1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]
    elif g[1] == 0:
        vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(-2*vtemp[g[1],g[0],1] + vtemp[(g[1]+1),g[0],1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]                                    
    elif g[1] == Ny-1:
        vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(vtemp[(g[1] -1),g[0],1] - 2*vtemp[g[1],g[0],1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]                         
    else:
        vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(vtemp[(g[1]-1),g[0],1] - 2*vtemp[g[1],g[0],1] + vtemp[(g[1]+1),g[0],1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]
    
###########################################################
#  Produces a Nx x Ny matrix of healthy cells (h)         #
#  Produces a time matrix for after eclipse phase (e)     #
#  Produces a time matrix for after infection phase (i)   #
#  Produces a matrix that holds amount of virus (vir)     #
#  Produces a time matrix hor healthy cells (t)           #
#  Produces a univeral time matrix (uv)                   #
###########################################################

cells = numpy.empty([Ny,Nx],dtype=str )
cells.fill('h')
ecl = numpy.zeros((Ny,Nx))
test = numpy.zeros((Ny,Nx))
inf = numpy.zeros((Ny,Nx))
vir = numpy.zeros((Ny,Nx,Nt))
vtemp = numpy.zeros((Ny,Nx,3))
th = numpy.zeros((Ny,Nx))
ut = numpy.zeros((Ny,Nx))
tsmatrix = numpy.empty([Ny,Nx],dtype=float)
tsmatrix.fill(ts)


################################################################
#   Deletes past versions of cells_over_time                   #
################################################################
deleteContent('cells_over_time.txt')

################################################################
#   Infects a random cell, now seen as (i)                     #
################################################################

randx = numpy.random.randint(1,Nx)
randy = numpy.random.randint(1,Ny)
cells[randy,randx] = 'i'
inf[randy,randx] = Ti()
f = 0

###################################################################
#                                                                 #
#                        Runs simulation                          #
#                                                                 #
###################################################################

step = itertools.cycle('ABCDEFG')
while (cells.all != 'd') is True:
#####################################
#       The Universal Time          #
#       is kept here (ut)           #
#####################################
    next(step)
    r = r+1
    if next(step) == 'A':
        ut = ut + tsmatrix
        #print ut
        #print 'ut'
#####################################
#       The Healthy Cells' time     #
#       is kept here (th)           #
#####################################
    if next(step) == 'B':
        hindex = numpy.where(cells == 'h', ts, 0)    
        th = th + hindex
        #print th
        #print 'th'
#####################################
#    Eclipse phase -> Infection     #
#                                   #
#####################################
    if next(step) == 'C':
        eindex = numpy.where(cells == 'e')[1]
        eindey = numpy.where(cells == 'e')[0]
        eindexy = numpy.vstack((eindex, eindey))
        w = range(len(eindex))
        if eindex.size != 0:
            for j in w:
                g = eindexy[:,j]
                if (ecl[g[1],g[0]] + th[g[1],g[0]]) < ut[g[1],g[0]]:
                    cells[g[1],g[0]] = 'i'
                    inf[g[1],g[0]] = inf[g[1],g[0]] + Ti()
#####################################
#       Infection spreads &         #
#       Infectious cells die        #
#####################################
#       g[0] = x coordinate          #
#       g[1] = y coordinate          #
#####################################
    if next(step) == 'D':
        
        iindex = numpy.where(cells == 'i')[1]
        iindey = numpy.where(cells == 'i')[0]
        iindexy = numpy.vstack((iindex, iindey))
        w = range(len(iindex))
        if iindex.size != 0:
            for j in w:
                g = iindexy[:,j]
                if PG1()>0.5:
                    if g[0] == 0:
                        break
                    if cells[g[1],(g[0] -1)] == 'h':
                        cells[g[1],(g[0] -1)] = 'e'
                        ecl[g[1],(g[0] -1)] = ecl[g[1],(g[0] -1)] + Te()
                if PG1()>0.5:
                    if g[0] == (Nx-1):
                        break
                    if cells[g[1],(g[0] +1)] == 'h':
                        cells[g[1],(g[0] +1)] = 'e'
                        ecl[g[1],(g[0] +1)] = ecl[g[1],(g[0] +1)] + Te()
                if PG1()>0.5:
                    if g[1] == 0:
                        break
                    if cells[(g[1] -1),g[0]] == 'h':
                        cells[(g[1] -1),g[0]] = 'e'
                        ecl[(g[1] -1),g[0]] = ecl[(g[1] -1),g[0]] + Te()
                if PG1()>0.5:
                    if g[1] == (Ny-1):
                        break
                    if cells[(g[1] +1),g[0]] == 'h':
                        cells[(g[1] +1),g[0]] = 'e'
                        ecl[(g[1] +1),g[0]] = ecl[(g[1] +1),g[0]] + Te()
                if PG1()>0.5:                        
                    if g[0] == 0:
                        break
                    if g[1] == (Ny-1):
                        break
                    if cells[(g[1] +1),(g[0]-1)] == 'h':
                        cells[(g[1] +1),(g[0]-1)] = 'e'
                        ecl[(g[1] +1),(g[0]-1)] = ecl[(g[1] +1),(g[0]-1)] + Te()
                if PG1()>0.5:                        
                    if g[0] == (Nx-1):
                        break
                    if g[1] == (Ny-1):
                        break
                    if cells[(g[1] +1),(g[0]+1)] == 'h':
                        cells[(g[1] +1),(g[0]+1)] = 'e'
                        ecl[(g[1] +1),(g[0]+1)] = ecl[(g[1]+1),(g[0]+1)] + Te()                       
                if ut[g[1],g[0]] > (inf[g[1],g[0]] + ecl[g[1],g[0]] + th[g[1],g[0]]):
                    cells[g[1],g[0]] = 'd'
            
        #print (inf + ecl + th)
        #print 'total'
        #print ut
        #print 'ut'
#####################################
#       Prints status of cells      #
#                                   #
#####################################
    if next(step) == 'E':
        print cells
        
        eclipse = open('eclipse.txt', 'w')
        healthy = open('healthy.txt', 'w')
        infected = open('infected.txt', 'w')
        time = open('time.txt', 'w')        
        cells_over_time = open('cells_over_time.txt', 'a')
        
        for i in cells:
            cells_over_time.write(i)
            dumb = dumb + 1
            if dumb % Nx-1 == 0:
                cells_over_time.write('\n')
                dumb = 0
            
        
        numpy.savetxt('eclipse.txt', ecl)
        numpy.savetxt('healthy.txt', th)
        numpy.savetxt('infected.txt', inf)
        numpy.savetxt('time.txt', ut)
        
        cells_over_time.close()
        eclipse.close()
        healthy.close()
        infected.close()
        time.close()
        print r
        if r > 30000:
            sys.exit('DONE')

#####################################
#       differential equation       #
#       for viral diffusion         #
#####################################
    if next(step) == 'F':
        vindex = numpy.where(cells != 'x')[1]
        vindey = numpy.where(cells != 'x')[0]
        vindexy = numpy.vstack((vindex, vindey))
        w = range(len(vindex))
        if vindex.size != 0:
            for j in w:
                g = vindexy[:,j]
                vtemp[g[1],g[0],1] = vtemp[g[1],g[0],2]
                if f % scale == 0:
                    vir[g[1],g[0],(f/scale)] = vtemp[g[1],g[0],1]
                    print "divisible"
                
                if cells[g[1],g[0]] == 'i':                 #Is the cell infected?#
                    rho2 = rho
                else:
                    rho2 = 0

                if g[0] == 0:
                    vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(-2*vtemp[g[1],g[0],1] + vtemp[g[1],(g[0]+1),1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]                                    
                elif g[0] == Nx-1:
                   vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(vtemp[g[1],(g[0] -1),1] - 2*vtemp[g[1],g[0],1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]
                elif g[1] == 0:
                    vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(-2*vtemp[g[1],g[0],1] + vtemp[(g[1]+1),g[0],1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]                                    
                elif g[1] == Ny-1:
                    vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(vtemp[(g[1] -1),g[0],1] - 2*vtemp[g[1],g[0],1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]                         
                else:
                    vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(vtemp[(g[1]-1),g[0],1] - 2*vtemp[g[1],g[0],1] + vtemp[(g[1]+1),g[0],1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]

#####################################
#       differential equation       #
#       for viral diffusion         #
#####################################
    if next(step) == 'Z':
        vindex = numpy.where(cells != 'x')[1]
        vindey = numpy.where(cells != 'x')[0]
        vindexy = numpy.vstack((vindex, vindey))
        w = range(len(vindex))
        if vindex.size != 0:
            for j in w:
                g = vindexy[:,j]
                vtemp[g[1],g[0],1] = vtemp[g[1],g[0],2]
                if f % scale == 0:
                    vir[g[1],g[0],(f/scale)] = vtemp[g[1],g[0],1]
                    print "divisible"
                
                if cells[g[1],g[0]] == 'i':                 #Is the cell infected?#
                    rho2 = rho
                else:
                    rho2 = 0

                if g[0] == 0:
                    vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(-2*vtemp[g[1],g[0],1] + vtemp[g[1],(g[0]+1),1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]                                    
                elif g[0] == Nx-1:
                   vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(vtemp[g[1],(g[0] -1),1] - 2*vtemp[g[1],g[0],1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]
                elif g[1] == 0:
                    vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(-2*vtemp[g[1],g[0],1] + vtemp[(g[1]+1),g[0],1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]                                    
                elif g[1] == Ny-1:
                    vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(vtemp[(g[1] -1),g[0],1] - 2*vtemp[g[1],g[0],1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]                         
                else:
                    vtemp[g[1],g[0],2] = (rho2 - c*vtemp[g[1],g[0],1] + D*(vtemp[(g[1]-1),g[0],1] - 2*vtemp[g[1],g[0],1] + vtemp[(g[1]+1),g[0],1])*(1/deltxprime**2))*ts + vtemp[g[1],g[0],1]

#####################################
#       Probability for virus       #
#       to infect (diffusion)       #
#####################################
        if next(step) == 'G':
            vindex = numpy.where(cells != 'x')[1]
            vindey = numpy.where(cells != 'x')[0]
            vindexy = numpy.vstack((vindex, vindey))
            w = range(len(vindex))
            for j in w:
                g = vindexy[:,j]
                pinfect = vtemp[g[1],g[0],2]*D*ts
                if PU1() > pinfect:
                    if cells[g[1],g[0]] == 'h':
                        cells[g[1],g[0]] == 'i'
##########################################
#        Prints v/t graph and exits      #
##########################################                  
        if next(step) == 'H':
            f=f+1
            if f == Nt*scale - 1:
                for i in range(0,Nx):
                    t = numpy.arange(0,Nt,1)
                    v = vir[0,i,0:Nt]
                    print v
                    print 'v'
                    p.title("cell " + str(i))
                    p.plot(t,numpy.log(v))
                    p.show()
                print ecl[:,:]
                print 'ecl'
                print th[:,:]
                print "th"
                print ut[:,:]
                print "ut"
                sys.exit("PROGRAM DONE! TIME LIMIT REACHED")

#Notes:
#D(delta*t)/h^2<1/2
#Heisenburg for diffusion^

#FIX:
#numpy.unique to replace where function in step F
# look into using two for loops for ^^^
#find out how many cores you have and parallelize
#try really hard to exit program when cells die
#Look at cells_over_time for weirdly long lasting infectious phase
#In your photos^^^
#Email Parker virus matrix
#Put comments everywhere