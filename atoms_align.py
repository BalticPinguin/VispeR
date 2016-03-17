#!/usr/bin/python2
# filename atoms_align.py
import numpy as np
import random

class align_atoms():
   Angs2Bohr=1/0.52917721092                                  
   CartCoord=[]
   spect=[]
   log=[]
   dim=0

   def __init__(self, method, spect):
      if method in ["moi", "MOI"]:
         self.type='moi'
      if method in ["rmsd" ,"RMSD"]:
         self.type='rmsd'
      self.CartCoord=spect.CartCoord
      self.log=spect.log
      self.spect=spect
      self.dim=spect.dim

   def perform(self):
      if self.type=='rmsd':
         self.RMSD_reorient()
      elif self.type=='rmsd':
         self.MOI_reorient()
      #now, copy the coordinates back to parent-class
      self.spect.CartCoord=self.CartCoord

   def shift(self):
      COM=np.zeros(3)
      X=np.zeros( (2,3,3) )
      diagI=np.zeros( (2,3) )
      #do for initial and final state:
      for i in [0,1]:
         #loop over coordinates:
         for j in [0,1,2]:
            COM[j]=np.sum(self.CartCoord[i][j]*self.spect.mass)
            COM[j]/=np.sum(self.spect.mass) 
         #now it is Cartesian center of mass
         if self.log.level<2:
            if i==0:
               self.log.write("Center of mass (initial) coordinates (Bohr):\n")
            else:
               self.log.write("Center of mass (final) coordinates (Bohr):\n")
            self.log.printVec(COM)
         for j in [0,1,2]:
            #displacement of molecule into center of mass:
            self.CartCoord[i][j]-=COM[j]

   def MOI_reorient(self):
      """This function reorients the final state in space such that
         the moment of inertia frames coincide.
         The function will fail when the order of the moments changes
         or are close to degenerate and hence are rotated towards each other.
         FIXME: I want to add a threshold to make the programme decide itself,
         whether this method is applicable or not.
      """
      #FIRST STEP: move the center of mass (COM) to origin:
      self.shift()

      #SECOND STEP: Calculate the inertia-system.
      # do the following for the initial and final state:
      for i in [0,1]:
      
         #  (MOI: Moment Of Inertia)
         MOI=np.zeros((3,3))# this is Moment Of Inertia
         #loop over coordinates:
         for j in [0,1,2]:
            MOI[j][j]=np.sum(self.spect.mass*self.spect.mass*(self.CartCoord[i][0]*self.CartCoord[i][0]+\
                     self.CartCoord[i][1]*self.CartCoord[i][1]+self.CartCoord[i][2]*self.CartCoord[i][2]-\
                     self.CartCoord[i][j]*self.CartCoord[i][j]))
            for k in range(j):
               MOI[j][k]=np.sum(self.spect.mass*self.spect.mass*(self.CartCoord[i][j]*self.CartCoord[i][k]))
               MOI[k][j]=MOI[j][k]
         #calculate the eigen-system of MOI
         diagI[i],X[i]=np.linalg.eig(MOI) 
         # sort it such, that the rotation is minimal. Therefore, first check
         # that it is not too big; otherwise, this simple code might fail.
         index=np.argsort(diagI[i], kind='heapsort')
         X[i]=(X[i].T[index]).T
         diagI[i]=diagI[i][index]
      #end for i in [0,1].
      
      #output on the information gathered above
      if self.log.level==0:
         self.log.write('\nRotational constants (GHz) in principle axes\n')
         self.log.write('initial state: '+ repr(1/(2*diagI[0].T)*self.Hartree2GHz)+"\n")
         self.log.write('final state: '+ repr(1/(2*diagI[1].T)*self.Hartree2GHz)+"\n")
         self.log.write("Inertia system in Cartesion Coordinates of initial state:\n")
         self.log.printMat(X[0])
         self.log.write("Inertia system in Cartesion Coordinates of final state:\n")
         self.log.printMat(X[1])
      
      #overlap of the moi-systems: gives the rotation of the respective frames
      # but is free in sign; correct this in the following:
      O=X[0].dot(X[1].T) 

      rmsdi=[]
      # now: test all combinations of signs. criterion is the least square of 
      # Cartesian coordinates.
      sign=np.eye(3)
      for i in range(13):
         #indeed, this gives all combinations 
         # after another...
         sign[int((i//3+i)%3)][int((i//3+i)%3)]*=-1
         if i in [4,7,8,10, 11]:
            #they give redundant results.
            continue
         U=sign.dot(O.T)
         rmsdi.append(self.RMSD(self.CartCoord[0]-U.dot(self.CartCoord[1])))

         #print  self.RMSD(self.CartCoord[0]-self.CartCoord[1]), self.RMSD(self.CartCoord[0]-U.dot(self.CartCoord[1]))
      rmsd=RMSD(self.CartCoord[0]-self.CartCoord[1])
      rmsdi.append(rmsd)
      #reassign i 
      i=np.argmin(rmsdi)
      if i==len(rmsdi)-1:
         # no rotation into MOI-frame should be done because initial 
         # geometry is the best. This is the case especially, if there
         # are (almost) degenerate moments of intertia.
         return
      #recover the sing-combination with least squares:
      if i>=4:
         i+=1
         if i>=7:
            i+=1
            if i>=8:
               i+=1
               if i>=10:
                  i+=1
                  if i>=11:
                     i+=1
      sign=np.eye(3)
      for j in range(i):
         sign[int((j//3+j)%3)][int((j//3+j)%3)]*=-1
      U=sign.dot(O)
      #apply this combination to coordinates of final state
      self.apply_change(U)
   
      # finally: print what is done.
      if self.log.level==0:
         self.log.write("Rotation of final state:\n")
         self.log.printMat(U)
         self.log.write("Coordinates after manipulation::\n")
         self.log.write('Cartesian coordinates of initial state: \n')
         self.log.printMat(self.CartCoord[0].T/self.Angs2Bohr)
         self.log.write('Cartesian coordinates of final state: \n')
         self.log.printMat(self.CartCoord[1].T/self.Angs2Bohr)
   
   def apply_change(self, U):
      """ produce respective matrix to rotate Force-constant matrix 
          and gradient (if given) as well.
          This function is only needed by RMSD_reorient and MOI_reorient.
      """
      self.CartCoord[1]=U.dot(self.CartCoord[1])
      Y=np.zeros( (self.dim,self.dim) )
      for j in range(self.dim//3):
         Y[3*j:3*j+3].T[3*j:3*j+3]=U
      # apply the respective rotation:
      self.spect.F[1]=np.dot(Y.dot(self.spect.F[1]),Y.T)
      if any(self.spect.Grad[i]>0 for i in range(len(self.spect.Grad))):
        # temp=np.zeros( (3,len(self.spect.Grad)//3) )
        # for i in [0,1,2]:
        #    for j in range(len(temp)):
        #       temp[i][j]=self.spect.Grad[j*3+i]
        # temp=U.dot(temp)
        # for i in [0,1,2]:
        #    for j in range(len(temp)):
        #       self.spect.Grad[j*3+i]=temp[i][j]

         self.spect.Grad=Y.dot(self.spect.Grad)

        # for j in range(self.dim//3):
        #    self.spect.Grad[3*j:3*j+3]=U.dot(self.spect.Grad[3*j:3*j+3])

   def RMSD_reorient(self):
      """This function reorients the final state in space such that
         the moment of inertia of coordinates is minimized.
         I assume that there might be coordinate flips and/or switches
         but no strong rotations by ~40 degree
      """
      #FIRST STEP: move the center of mass (COM) to origin:
      self.shift()
      
      # test all combinations of signs. and permutations:
      rotated=self.CartCoord[1]
      for j in range(6):
         #loop over all permutations:
         O=np.zeros( (3,3) )
         if j%2==0:
            #for cyclic permutations
            for i in range(3):
               #for cyclic ones, this works...
               O[(i+j//2)%3][(i-j//2)%3]=1
         else:
            #for anti-cyclic permutations
            for i in range(3):
               #for cyclic ones, this works...
               O[-((i+j//2)%3)-1][(i-j//2)%3]=1
         sign=np.eye(3)
         for i in range(13):
            #loop over all sing-combinations.
            #indeed, this gives all combinations 
            # after another...
            sign[int((i//3+i)%3)][int((i//3+i)%3)]*=-1
            if i in [4,7,8,10, 11]:
               #they give redundant results.
               continue
            U=sign.dot(O)
            #print np.shape(rotated) ,np.shape(self.CartCoord[0]), np.shape(U.dot(self.CartCoord[1]))
            if self.RMSD(self.CartCoord[0]-rotated) >self.RMSD(self.CartCoord[0]-U.dot(self.CartCoord[1])):
               #if the current combination is the best, save it.
               U_min=U
               rotated=U_min.dot(self.CartCoord[1])
               #print self.RMSD(self.CartCoord[0]-U.dot(self.CartCoord[1]))
      try:
         #if U_min is known: coordinate system changed.
         # This change is applied to the other quantities as well:
         self.apply_change(U_min)
      except NameError:
         #do nothing
         U=0

      #THIRD STEP: follow some Monte Carlo scheme.
      # comparatively small angles. Do 200 steps.
      U=np.eye(3)
      for i in range(40):
         #chose some angle:
         #for every level: do 5 tests and the refine the search...
         for j in range(50):
            #a rotation-matrix with small rotations:
            theta=0.063*(random.random()-.5)/(i+1)
            eta=0.063*(random.random()-.5)/(i+1)
            nu=0.063*(random.random()-.5)/(i+1)
            R_x=np.matrix([ [1,0,0],
                  [0,np.cos(theta),-np.sin(theta)],
                  [0,np.sin(theta),np.cos(theta)] ], dtype=float)
            R_y=np.matrix([[np.cos(eta),0,-np.sin(eta)],
                  [0,1,0],
                  [np.sin(eta),0,np.cos(eta)]], dtype=float)
            R_z=np.matrix([[np.cos(nu),-np.sin(nu),0],
                  [np.sin(nu),np.cos(nu),0],
                  [0,0,1] ], dtype=float)
            R=R_x.dot(R_y).dot(R_z)
            test=R.A.dot(U.dot(self.CartCoord[1]))
            if self.RMSD(self.CartCoord[0]-U.dot(self.CartCoord[1]))> self.RMSD(self.CartCoord[0]-test):
               #if it gets better: apply the change.
               #self.CartCoord[1]=test
               U=R.A.dot(U)
      
      #FOURTH STEP: apply the rotation.

      self.apply_change(U)
   
      # finally: print what is done.
      if self.log.level==0:
         self.log.write("Rotation of final state::\n")
         self.log.printMat(U)
         self.log.write("Coordinates after manipulation::\n")
         self.log.write('Cartesian coordinates of initial state: \n')
         self.log.printMat(self.CartCoord[0].T/self.Angs2Bohr)
         self.log.write('Cartesian coordinates of final state: \n')
         self.log.printMat(self.CartCoord[1].T/self.Angs2Bohr)
   # END OF FUNCTION DEFINITIONS

   def RMSD(self,Delta):
      """This function calculates the RMSD of a matrix (intended
         for self.CartCoords)
      """
      rmsd=0
      for i in range(len(Delta)):
         for j in range(len(Delta[0])):
            rmsd+=Delta[i][j]*Delta[i][j]
      return rmsd

#version 0.0
# End of atoms_align.py
