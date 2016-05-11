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
from multiprocessing import Pool

Nx = 3      # Range of cells on x axis
Ny = 3      # Range of cells on y axis
Nt = 100    # Days simulation runs for
ts = 1e-2   # Time step for model
TauI = 12   # Avg time for infection
TauE = 6    # Avg time for eclipse
ne = 60     # 
ni = 60     #
probi = 0.5  # Threshold for probablility of cell to cell infection
nmean = 0.5 # Mean of probability of cell to cell infection
scale = 10000 # scaling factor between k and ts
count = 0
#############

rho = 2e07
D = 1.1e-11
c = 0.056
deltx = 50e-06
deltxprime = deltx #*10e6
Dtsx2 = D*ts*deltxprime**-2
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
    
def diffVModel(j):         #Differential model for extracellular virus
    g = vindexy[:,j]
    vtemp[g[1],g[0],1] = vtemp[g[1],g[0],2]
    
    #only store virus data after a certain number of runs
    if f % scale == 0:
        vir[g[1],g[0],(f/scale)] = vtemp[g[1],g[0],1]
        print "divisible"
    
    #is the cell infected? (binary toggle)
    if cells[g[1],g[0]] == 'i':                 
        rho2 = rho
    else:
        rho2 = 0
    
    #unless dictated otherwise, all cells around target cell exist (binary toggle)
    existvrm = 1    
    existvcm = 1
    existvrp = 1
    existvcp = 1    
    
    vr = g[1]       #row coordinate of cell
    vc = g[0]       #column coordinate of cell
    vrm = g[1] -1   #row coordinate above cell
    vcm = g[0] -1   #column coordinate left of cell
    vrp = g[1] +1   #row coordinate below cell
    vcp = g[0] +1   #column coordinate right of cell
    
    if vrm < 0:         #if the cell one row up doesn't exist, it's taken out of the equation
        existvrm = 0
        vrm = 0
    if vcm < 0:         #if the cell one column to the left doesn't exist, it's taken out of the equation
        existvcm = 0
        vcm = 0
    if vrp > Ny-1:      #if the cell one row down doesn't exist, it's taken out of the equation
        existvrp = 0
        vrp = 0
    if vcp > Nx-1:      ##if the cell one column to the right doesn't exist, it's taken out of the equation
        existvcp = 0
        vcp = 0
    #ALL THIS WORK DONE SO THAT THE EQUATION IS GENERALIZED FULLY BELOW#
    vtemp[vr,vc,2] = ((rho2 - 4*Dtsx2)*vtemp[vr,vc,1] + ((2/3)*Dtsx2)(vtemp[vrm,vcm,1]*existvrm*existvcm + vtemp[vrm,vc,1]*existvrm + vtemp[vrm,vcp,1]*existvrm*existvcp + vtemp[vr,vcp,1]*existvcp + vtemp[vrp,vc,1]*existvrp + vtemp[vr,vcm,1]*existvcm))

def cell2cell(params):
    j, iindexy = params
    g = iindexy[:,j]
    
    #unless dictated otherwise, all cells around target cell exist (binary toggle)
    existrm = 1    
    existcm = 1
    existrp = 1
    existcp = 1   
    
    r = g[1]       #row coordinate of cell
    c = g[0]       #column coordinate of cell
    rm = g[1] -1   #row coordinate above cell
    cm = g[0] -1   #column coordinate left of cell
    rp = g[1] +1   #row coordinate below cell
    cp = g[0] +1   #column coordinate right of cell
    
    if rm < 0:         #if the cell one row up doesn't exist, it's taken out of the equation
        existrm = 0
        rm = 0
    if cm < 0:         #if the cell one column to the left doesn't exist, it's taken out of the equation
        existcm = 0
        cm = 0
    if rp > Ny-1:      #if the cell one row down doesn't exist, it's taken out of the equation
        existrp = 0
        rp = 0
    if cp > Nx-1:      ##if the cell one column to the right doesn't exist, it's taken out of the equation
        existcp = 0
        cp = 0
    
    if PG1()>0.5:
        if existcm == 1:
            if cells[r,cm] == 'h':
                cells[r,cm] = 'e'
                ecl[r,cm] = ecl[r,cm] + Te()
    if PG1()>0.5:
        if existcp == 1:
            if cells[r,cp] == 'h':
                cells[r,cp] = 'e'
                ecl[r,cp] = ecl[r,cp] + Te()
    if PG1()>0.5:
        if existrm == 1:
            if cells[rm,c] == 'h':
                cells[rm,c] = 'e'
                ecl[rm,c] = ecl[rm,c] + Te()
    if PG1()>0.5:
        if existrp == 1:
            if cells[rp,c] == 'h':
                cells[rp,c] = 'e'
                ecl[rp,c] = ecl[rp,c] + Te()
    if PG1()>0.5:
        if existrm == 1 and existcm == 1:
            if cells[rm,cm] == 'h':
                cells[rm,cm] = 'e'
                ecl[rm,cm] = ecl[rm,cm] + Te()
    if PG1()>0.5:
        if existrm == 1 and existcp == 1:
            if cells[rm,cp] == 'h':
                cells[rm,cp] = 'e'
                ecl[rm,cp] = ecl[rm,cp] + Te()
                       
    if ut[r,c] > (inf[r,c] + ecl[r,c] + th[r,c]):
        cells[r,c] = 'd'
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


#####################################
# Since vindexy is constant, we     #
# compute it first and refer later  #
#####################################
vindex = numpy.where(cells != 'x')[1]
vindey = numpy.where(cells != 'x')[0]
vindexy = numpy.vstack((vindex, vindey))


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
while cells.all != 'd':
#####################################
#       START PROCESS               #
#####################################  
    next(step)
    count = count+1
#####################################
#       The Universal Time          #
#       is kept here (ut)           #
#####################################
    if next(step) == 'A':
        print "A"
        ut = ut + tsmatrix
#####################################
#       The Healthy Cells' time     #
#       is kept here (th)           #
#####################################
    if next(step) == 'B':
        print "B"
        hindex = numpy.where(cells == 'h', ts, 0)    
        th = th + hindex
        #print th
        #print 'th'
#####################################
#    Eclipse phase -> Infection     #
#                                   #
#####################################
    if next(step) == 'C':
        print "C"
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
    if next(step) == 'D':
        print "D"
        iindex = numpy.where(cells == 'i')[1]
        iindey = numpy.where(cells == 'i')[0]
        iindexy = numpy.vstack((iindex, iindey))
        w = range(len(iindex))
        if __name__ == '__main__':
            params = zip(w, iindexy)
            pool = Pool()
            pool.map(cell2cell, params)
        
        
#####################################
#       Prints status of cells      #
#                                   #
#####################################
    if next(step) == 'E':
        print "E"
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
        print count
        if count > 30000:
            sys.exit('DONE')

#####################################
#       differential equation       #
#       for viral diffusion         #
#####################################
    if next(step) == 'F':
        print "F"
        w = range(len(vindex))
        if __name__ == '__main__':
            pool = Pool()
            pool.map(diffVModel, w)
            #for j in w: distributed

                

#####################################
#       Probability for virus       #
#       to infect (diffusion)       #
#####################################
        if next(step) == 'G':
            print "G"
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
            print "H"
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