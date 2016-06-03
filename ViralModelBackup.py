# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 09:53:35 2015

@author: dylanbarth
"""

##############################
#####   Virus Model     ######
##############################
#
# THIS COPY WORKED WELL AND WAS THOROUGHLY TESTED BEFORE ADAPTIVE TIME STEPS
# ALSO BEFORE A SECOND ATTEMPT AT PARALLELIZATION
# ACCOMPANIED WITH IT IS A GRAPHING TOOL THAT READS "vsum_over_time.txt"
# 
####################################
# Importing modules and packages   #
####################################
import itertools
import sys
import numpy
import math
import gc
import datetime
#from multiprocessing import Pool





####################################
# Defining variables and functions #
####################################
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
def now():
    return str(datetime.datetime.now())


Nx = 10      # Range of cells on x axis
Ny = 10      # Range of cells on y axis
Nt = 100    # data points taken
ts = 5e-3   # Time step for model (hours)
endt = 100 # How long (hours) simulated time == (Nt*scale*ts)
TauI = 12.   # Avg time for infection
TauE = 6.    # Avg time for eclipse
ne = 60     # 
ni = 60     #
probi = 0.5  # Threshold for probablility of cell to cell infection
nmean = 0.5 # Mean of probability of cell to cell infection
scale = endt/(Nt*ts)  #scaling factor between Nt and count
count = 0   # equal to the number of ts elapsed

#############

beta = 1.33e-06
rho = 0.00192
D = 3.96e-8
c = 0.056
deltx = 25.0e-06
deltxprime = deltx#*10e6
Dtsx2 = D*ts*(deltxprime**-2)
dumb = 0
dummy = 0
f = 0
#############################


cells = numpy.empty([Ny,Nx],dtype=str )  
cells.fill('h') #Produces a Nx x Ny matrix of healthy cells (h)
ecl = numpy.zeros((Ny,Nx)) #Produces a time matrix for after eclipse phase (e)
inf = numpy.zeros((Ny,Nx)) #Produces a time matrix for after infection phase (i)
vir = numpy.zeros((Ny,Nx,Nt)) #Produces a matrix that holds amount of virus (vir)
vtemp = numpy.zeros((Ny,Nx,3)) 
th = numpy.zeros((Ny,Nx)) #Produces a time matrix hor healthy cells (t)
ut = numpy.zeros((Ny,Nx)) #Produces a univeral time matrix (ut)
tsmatrix = numpy.empty([Ny,Nx],dtype=float)
tsmatrix.fill(ts) #Produces a matrix filled with value of time step

#####################################
# Since vindexy is constant, we     #
# compute it first and refer later  #
#####################################
vindex = numpy.where(cells != 'x')[1]
vindey = numpy.where(cells != 'x')[0]
vindexy = numpy.vstack((vindex, vindey))

################################################################
#   Infects a random cell, now seen as (i)                     #
################################################################

randx = numpy.random.randint(0,Nx)
randy = numpy.random.randint(0,Ny)
cells[randy,randx] = 'i'
inf[randy,randx] = Ti()


################################################################
#   Deletes past versions of cells_over_time                   #
#   and fills it with a matrix of healthy cells                #
################################################################
deleteContent('cells_over_time.txt')

cells_over_time = open('cells_over_time.txt','w')
for i in cells:
    cells_over_time.write(i)
    dumb = dumb + 1
    if dumb % Nx-1 == 0:
        cells_over_time.write('\n')
        dumb = 0
dumb = 0

################################################################
#   Writes a file with all of our parameters/variables         #
################################################################
varb = open('parameters.txt','w')
varb.write(
"Nx = " + str(Nx) + "\nNy = " + str(Ny) + "\nNt = " + str(Nt) 
+ "\nts = " + str(ts) + "\nscale = " + str(scale) + "\nrho = " + str(rho) 
+ "\nD = " + str(D) + "\ndelta x = " + str(deltxprime) + "\nc = " + str(c)
+ "\nTotal Hours Simulated = " + str(endt) + "\nReal Time Start = " + now()
)
varb.close()

################################################################
#   Checks for Heisenberg status of viral diffusion            #
################################################################
if (D*ts/(deltxprime**2) > 0.5) is True:
    print D*ts/(deltxprime**2)
    sys.exit("CHANGE PARAMETERS TO FIT DIFFUSION LIMITS. VALUE MUST BE UNDER 0.5. VALUE SHOWN ABOVE")

if (scale).is_integer() == False:
    sys.exit("scale must be integer, replace values for Nt, ts, or endt")
    

###################################################################
#                                                                 #
#                        Runs simulation                          #
#                                                                 #
###################################################################

step = itertools.cycle('ABCDEFGHI')
while cells.all != 'd':
#####################################
#       START PROCESS               #
#####################################  
    s=next(step)
#####################################
#       The Universal Time          #
#       is kept here (ut)           #
#####################################
    if s == 'A':
        #print "A"
        ut = ut + tsmatrix
        count = count+1
#####################################
#       The Healthy Cells' time     #
#       is kept here (th)           #
#####################################
    if s == 'B':
        #print "B"
        hindex = numpy.where(cells == 'h', ts, 0)    
        th = th + hindex
#####################################
#    Eclipse phase -> Infection     #
#                                   #
#####################################
    if s == 'C':
        #print "C"
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
    if s == 'D':
        #print "D"
        iindex = numpy.where(cells == 'i')[1]
        iindey = numpy.where(cells == 'i')[0]
        iindexy = numpy.vstack((iindex, iindey))
        if iindex.size != 0:
            for j in range(len(iindex)):
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
                            
                # Adjustment for even/odd rows #
                if g[0] % 2 == 0:
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
                                
                if g[0] % 2 != 0:
                    if PG1()>0.5:
                        if existrm == 1 and existcm == 1:
                            if cells[rp,cm] == 'h':
                                cells[rp,cm] = 'e'
                                ecl[rp,cm] = ecl[rp,cm] + Te()
                                
                    if PG1()>0.5:
                        if existrm == 1 and existcp == 1:
                            if cells[rp,cp] == 'h':
                                cells[rp,cp] = 'e'
                                ecl[rp,cp] = ecl[rp,cp] + Te()
                                
                                   
    
            
#####################################
#       Prints status of cells      #
#                                   #
#####################################
    if s == 'E':
        #print "E"
        
        if count % scale == 0:
            cells_over_time = open('cells_over_time.txt','a')
            for i in cells:
                cells_over_time.write(i)
                dumb = dumb + 1
                if dumb % Nx-1 == 0:
                    cells_over_time.write('\n')
                    dumb = 0
            cells_over_time.close()
        #print count
        

#####################################
#       differential equation       #
#       for viral diffusion         #
#####################################
    if s == 'F':
        #print "F"
        #print cells
        if count % scale == 0:
            dummy = dummy + 1
        #if __name__ == '__main__':
        #    pool = Pool(2)
        #    pool.map(diffVModel, range(len(vindex)))
        #    diffVModel(range(len(vindex)))
        for j in range(len(vindex)):
            g = vindexy[:,j]
#            print(vtemp[g[1],g[0],1])
#            print(g)
#            print(iindexy)
            
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
            
            #is the cell infected? (binary toggle)
            #print vr
            #print vc
            #print cells[vr,vc]
            #print 'CELL'
            if cells[vr,vc] == 'i':                 
                rho2 = rho
            else:
                rho2 = 0
            
            if vrm < 0:         #if the cell one row up doesn't exist, it's taken out of the equation
                existvrm = 0
                vrm = vr
            if vcm < 0:         #if the cell one column to the left doesn't exist, it's taken out of the equation
                existvcm = 0
                vcm = vc
            if vrp > Ny-1:      #if the cell one row down doesn't exist, it's taken out of the equation
                existvrp = 0
                vrp = vr
            if vcp > Nx-1:      ##if the cell one column to the right doesn't exist, it's taken out of the equation
                existvcp = 0
                vcp = vc 
            if g[0] % 2 ==0:
                checknn =  existvcm*existvrm + existvrm + existvcp*existvrm + existvcp + existvrp + existvcm
                NNN = ( vtemp[vrm,vcm,1] + vtemp[vrm,vc,1] + vtemp[vrm,vcp,1] + vtemp[vr,vcp,1] + vtemp[vrp,vc,1] + vtemp[vr,vcm,1] )
            else:
                checknn =  existvcm + existvrm + existvcp + existvcp*existvrp + existvrp + existvcm*existvrp
                NNN = ( vtemp[vr,vcm,1] + vtemp[vrm,vc,1] + vtemp[vr,vcp,1] + vtemp[vrp,vcp,1] + vtemp[vrp,vc,1] + vtemp[vrp,vcm,1] ) 
            
            #ALL THIS WORK DONE SO THAT THE EQUATION IS GENERALIZED FULLY BELOW#
            vprod = rho2*ts
            vleft = 4*Dtsx2*vtemp[vr,vc,1]
            vin = 2*Dtsx2*NNN/3
            vclear = c*vtemp[vr,vc,1]*ts
#            print((vprod,vleft,vin,vclear))
#            vtemp[vr,vc,2] = (rho2*ts +(1 - 4*Dtsx2)*vtemp[vr,vc,1] + 2*Dtsx2*NNN/3 - c*vtemp[vr,vc,1]*ts)
            vtemp[vr,vc,2] = vtemp[vr,vc,1] + vprod - vleft + vin - vclear
            if vtemp[vr,vc,2] < 0:
                vtemp[vr,vc,2] = 0.
                
            
        #only store virus data after a certain number of runs
            if count % scale == 0:
                vir[vr,vc,dummy-1] = vtemp[vr,vc,2]

        vtemp[:,:,1] = vtemp[:,:,2]

#####################################
#       Probability for virus       #
#       to infect (diffusion)       #
#####################################
    if s == 'G':
        #print "G"
        w = range(len(vindex))
        for j in w:
            g = vindexy[:,j]
            pinfect = vtemp[g[1],g[0],2]*D*ts
            if PU1() > pinfect:
                if cells[g[1],g[0]] == 'h':
                    cells[g[1],g[0]] = 'e'
                    ecl[g[1],g[0]] = ecl[g[1],g[0]] + Te()
                    
##########################################
#               kills cells              #
##########################################                    
    if s == 'H':
        iindex = numpy.where(cells == 'i')[1]
        iindey = numpy.where(cells == 'i')[0]
        iindexy = numpy.vstack((iindex, iindey))
        if iindex.size != 0:       
            for j in range(len(iindex)):
                g = iindexy[:,j]
                r = g[1]
                c = g[0]
                if ut[r,c] > (inf[r,c] + ecl[r,c] + th[r,c]):
                    cells[r,c] = 'd'

##########################################
#     Prints output files and exits      #
##########################################                  
    if s == 'I':
        #print "H"
        if count == Nt*scale - 1:
            yaxis = []	
            for i in range(0,Nt-1):
                v = vir[:,:,i]
                yaxis.append(numpy.sum(v))
            virus_time = open('vsum_over_time.txt','w') 
            virus_time.write(str(yaxis))
            print yaxis
            
            varb = open('parameters.txt','a')
            varb.write("\nReal Time End = " + now())
            varb.close()

            v_text = open('virus_over_time.txt','w')            
            eclipse = open('eclipse.txt', 'w')
            healthy = open('healthy.txt', 'w')
            infected = open('infected.txt', 'w')
            time = open('time.txt', 'w')        
            
            numpy.savetxt('eclipse.txt', ecl)
            numpy.savetxt('healthy.txt', th)
            numpy.savetxt('infected.txt', inf)
            numpy.savetxt('time.txt', ut)
            numpy.savetxt('virus_over_time.txt', vir)
            
            eclipse.close()
            healthy.close()
            infected.close()
            time.close()
            v_text.close()
            
            sys.exit("PROGRAM DONE! TIME LIMIT REACHED")
########################################## 
#   garbages unused/unreferenced memory  #
##########################################       
    if s == 'J':
        gc.collect()
        #progress bar#
        if count % Nt*scale*(1/10) == 0.0:
            print(str(count*100/(Nt*scale)) + "% DONE")

                
#Notes:
#D(delta*t)/h^2<1/2
#Heisenburg for diffusion^

#FIX:
#try really hard to exit program when cells die
#adaptive time step DX
#parallelize again..
