#!/usr/bin/python2
# filename: DR_spects.py

#include [[Spect.py]]
import Spect

#include [[Btree.py]]
import Btree as bt

from copy import deepcopy
import numpy as np
from ctypes import *
import re, mmap, math, os

#  CHANGELOG 
# ===========
#in version 1.0:  
#  1) remove dist_FCfOPA(); the model behind is inconsistent.   
#  2) renamed resortFCfOPA to calcspect and made the interface
#     accordingly.  
#  3) SDR_spect class runs now
#  4) Changed makeLine in SDR_spect to be compatible with the MultiPart
#     class. The model is inconsistent but may be still better.
#

class URDR_spect(Spect.Spect):
   """The class to calculate full Duschinsky-spectra. Due to its computational structure,
      integral prescreaning is not possible, therefore really all transitions need to be 
      computed.
   """
   Hartree2cm_1=219474.63 
   Threshold=3e-10

   m = 10
   A=[]
   b=[]
   C=[]
   d=[]
   E=[]
   Gamma =[]
   Gammap=[]

   def __init__(self, f): 
      """In addition to the initialisation of
         Spect( see first line of code), here the number
         of modes is evaluated and the matrices/vectors
         needed for evaluation of the recursive equations are
         calculated.
      """
      # first, initialise as done for all spectra:
      self.type='URDR'
      Spect.Spect.__init__(self, f)
      # now, calculate additional quantities for the iteration:
      self.__GetQuants()
      #get m from opt. (-> number of states considered)
      modes=re.findall(r"(?<=modes\=)[\d]+",self.opt, re.I)
      if modes!=[]:
         self.m=modes[-1]

   def __GetQuants(self):
      """This function computes the matrices and vectors needed interally for the computation
         of intensities.
      """
      self.Gamma=np.diag(self.f[0])   # in atomic units. It is equivalent to 4pi^2/h f_i
      self.Gammap=np.diag(self.f[1])  # for initial state
      sqGamma=np.diag(np.sqrt(self.f[0]))
      sqGammap=np.diag(np.sqrt(self.f[1]))
      unity=np.eye(len(self.Gamma))

      TMP=np.linalg.inv(self.nm.J.T.dot(self.Gammap).dot(self.nm.J) + self.Gamma)
      self.A=2.*sqGammap.dot(self.nm.J).dot(TMP).dot(self.nm.J.T).dot(sqGammap) -unity
      self.b=2.*sqGammap.dot( unity - self.nm.J.dot(TMP).dot(self.nm.J.T).dot(self.Gammap) ).dot(self.nm.K)
      self.C=2.*sqGamma.dot(TMP).dot(sqGammap) -unity
      self.d=-2.*sqGamma.dot(TMP).dot(self.nm.J.T.dot(self.Gammap.dot(self.nm.K)))
      self.E=4.*sqGamma.dot(TMP).dot(self.nm.J.T).dot(sqGammap)

   def calcspect(self):
      """This function prepares all variables according to the given 
         options to calculate the spectrum.
      """
      #first: resort the elements of J, K, f to make J most closely to unity
      resort=np.zeros(np.shape(self.nm.J))
      for i in xrange(len(self.nm.J)):
         j=np.argmax(self.nm.J[i])
         k=np.argmin(self.nm.J[i])
         if self.nm.J[i][j]>-self.nm.J[i][k]:
           resort[i][j]=1
         else:
            resort[i][k]=-1
      self.nm.J=resort.dot(self.nm.J.T)
      self.nm.K=resort.dot(self.nm.K.T)
      for i in xrange(len(resort)):
         k=np.argmin(resort[i])
         if resort[i][k]==-1:
            resort[i][k]=1 #use absolute value only.
      self.f[1]=resort.dot(self.f[1].T)
      self.log.write("before truncation:\n")
      self.log.write("Duschinsky-Matrix:\n")
      self.log.printMat(self.nm.J)
      self.log.write("Displacement vector:\n")
      self.log.printVec(self.nm.K)
   
      #truncate only, if this really is truncation.
      if self.m<len(self.nm.J): 
         index=np.argsort(np.abs(self.nm.K), kind="heapsort")[::-1]
         ind=index[:self.m]
         self.nm.K=self.nm.K[ind]
         self.nm.J=self.nm.J[ind].T[ind]
         f=np.zeros(( 2,self.m))
         f[1]=self.f[1][ind]
         f[0]=self.f[0][ind]
         self.f=f
         # now recalculate the matrices A,C,E,...
         self.__GetQuants()
      self.log.write("After truncation:\n")
      self.log.write("Duschinsky-Matrix:\n")
      self.log.printMat(self.nm.J)
      self.log.write("Displacement vector:\n")
      self.log.printVec(self.nm.K)
   
      # finally, calculate the stick spectrum in this picture
      self.spect=self.FCf(self.states1+ self.states2, self.T)
   
   def FCf(self,N, T):
      """Calculates the intensities including the Duschinsky-effect. 
                     
         **PARAMETERS:**  
         N:      Max. number of excitation quanta state considered  
         T:      Temperature of the system
            
         All parameters are obligatory.
         
         **RETURNS:**  
         linespectrum 
      """
      def CalcI00():
         """This function calculates the overlap-integral for zero vibrations 
         """
   
         Tree=bt.Tree(2*len(self.nm.K))
         Tree.fill(0)
         Zero=np.zeros(2*len(self.nm.K))
         # correct for vibrational groundstates
         E=abs(self.Energy[0]-self.Energy[1])
         for i in range(len(self.Gammap)):
            E+=.5*(self.Gammap[i][i]-self.Gamma[i][i])
         Tree.insert(Zero, np.array([1.0, (self.Energy[0]-self.Energy[1])*self.Hartree2cm_1, 0]) )
         #I_00 transition-probability [[Btree.py#extract]]
         #this is done using implicit side effects
         lines.append(1.0)
         freqs.append(E*self.Hartree2cm_1)
         initF.append(0) #needed for boltzmann-weighing
         return Tree
      
      def iterate(L1, L2, i):
         """This is the main function calculating the intensities iteratively.
            It fills the tree L3 using the previous trees L1 and L2.
            Each tree contains all transitions with a given total sum (initial + final)
            of vibrational quanta (i).
         """
         def freq( E, Gamma, Gammap):
            """Calculates the frequency of respective transition including 
            vibrational frequency
   
            *PARAMETERS:*
            E:      energy-difference of states
            Gamma:  vector of vibrational frequencies in inital state 
                    (in atomic units)
            Gammap: vector of vibrational frequencies in final state 
                    (in atomic units)
   
            *RETURNS;*
            frequency of respective transition
            """
            return (E+sum(Gammap-Gamma))*self.Hartree2cm_1 
      
         def FirstNonzero(n): 
            """Find first non-zero elements in first and second half of array n """
            leng=len(n)//2
            m=leng+1 #this means there is no excitation in this state
            mp=leng+1
            ni=n[:leng] #interger division (python3-compatible)
            nf=n[leng:]
            for j in range(leng):
               if ni[j]>0:
                  m=j
                  break
            for j in range(leng):
               if nf[j]>0:
                  mp=j
                  break
            return m, mp
         
         #initialize new tree
         alpha=2*len(self.b)
         L3=bt.Tree(i)
         L3.fill(alpha)
         States=states(alpha, i)           # States are all possible
         
         Energy=self.Energy[0]-self.Energy[1]
         sgnE=np.sign(Energy)
         if sgnE==0:
            sgnE=1
         npsqrt=np.sqrt
         leng=len(self.b) # = alpha//2
         for n in States: #for each possible state, described by n(vector)
            # index of first-non-zero element of (initial, final) state
            m, mp= FirstNonzero(n)
            # if the 'first' excited state is in initial state: 
            #          need first iteration formula
            if m<=mp:
               # need first iteration-formula
               n_m=n[m]
               n[m]-=1 #n[m] is at least 1
               Ps=L2.getState(n)[0]
               I_nn=self.b[m]*Ps                                     # first term 
               if n[m]>0:
                  n[m]-=1
                  Ps=L1.getState(n)[0]
                  I_nn+=npsqrt(2.0*(n_m-1))*self.A[m][m]*Ps          # second term
                  n[m]+=1
               #now, the two sums follow:
               if n[m+leng]>0:   #first summand in second sum:
                  n[m+leng]-=1
                  Ps=L1.getState(n)[0]
                  n[m+leng]+=1
                  I_nn+=npsqrt(n[m+leng]*0.5)*(self.E[m][m])*Ps # second term
               for k in range(m+1,leng):                           # rest of the sums
                  #first summand
                  if n[k]>0:    
                     n[k]-=1
                     Ps=L1.getState(n)[0]
                     n[k]+=1
                     I_nn+=npsqrt(n[k]*0.5)*(self.A[k][m]+self.A[m][k])*Ps # second term
                  #second summand
                  if n[k+leng]>0:
                     n[k+leng]-=1
                     Ps=L1.getState(n)[0]
                     n[k+leng]+=1
                     I_nn+=npsqrt(n[k+leng]*0.5)*(self.E[k][m])*Ps # second term
   
               n[m]+=1
            #else: need the other iteration-formula
            else: 
               n_m=n[mp+leng]
               n[mp+leng]-=1
               Ps=L2.getState(n)[0]
               I_nn=self.d[mp]*Ps                                       # first term 
               if n[mp+leng]>0:
                  n[mp+leng]-=1
                  Ps=L1.getState(n)[0]
                  I_nn+=npsqrt(2.0*(n_m-1.))*self.C[mp][mp]*Ps           # second term
                  n[mp+leng]+=1
               #compute terms in the two sums:
               for k in range(mp+1, leng):
                  if n[k]>0: #first sum
                     n[k]-=1
                     Ps=L1.getState(n)[0]
                     n[k]+=1
                     I_nn+=npsqrt(.5*n[k])*(self.C[mp][k]+self.C[k][mp])*Ps
                  if n[k+leng]>0: #second sum
                     n[k+leng]-=1
                     Ps=L1.getState(n)[0]
                     n[k+leng]+=1
                     I_nn+=npsqrt(.5*(n[k+leng]-1.))*self.E[mp][k]*Ps  

               n[mp+leng]+=1
            I_nn/=npsqrt(2.*n_m)
            assert not math.isnan(I_nn) ,"I_nn is not a number! I_nn:"+\
                                         " {0}\n, n:{1}\n:".format(I_nn, n)
            # don't mess with too low intensities (<1e-16)
            if np.abs(I_nn)>1e-8: 
               L3.insert(n, [I_nn,\
                freq(sgnE*Energy, sgnE*self.f[1]*n[:leng],\
                             sgnE*self.f[0]*n[leng:]),\
                         freq(0, 0, self.f[0]*n[leng:])])
         return L2, L3
      
      def states(alpha, n): 
         """This function creates all possible states having a total number of n 
            excitations in alpha different states
      
            **PARAMETERS:**   
            alpha: number of degrees of freedom  
            n:     number of excitations  
      
            **RETURNS:**  
            a vector containing the states.
         """
   
         #needed for 'states'
         def unlabeled_balls_in_labeled_boxes(balls, box_sizes): 
            """
            These functions are part of python-package: 'combinatorics' 
            (download from https://pypi.python.org/pypi/Combinatorics)
            unlabeled_balls_in_labeled_boxes(balls, box_sizes): This function 
            returns a generator that produces all distinct distributions of 
            indistinguishable balls
            among labeled boxes with specified box sizes (capacities). This is 
            a generalization of the most common formulation of the problem, 
            where each box is
            sufficiently large to accommodate all of the balls, and is an 
            important 
            example of a class of combinatorics problems called 
            'weak composition' problems.
         
            OVERVIEW
         
            This function returns a generator that produces all distinct 
            distributions of
            indistinguishable balls among labeled boxes with specified box sizes
            (capacities).  This is a generalization of the most common 
            formulation of the
            problem, where each box is sufficiently large to accommodate all of the
            balls, and is an important example of a class of combinatorics problems
            called 'weak composition' problems.
      
      
            CONSTRUCTOR INPUTS
         
            n: the number of balls
            
            box_sizes: This argument is a list of length 1 or greater.  The length of
            the list corresponds to the number of boxes.  `box_sizes[i]` is a positive
            integer that specifies the maximum capacity of the ith box.  If
            `box_sizes[i]` equals `n` (or greater), the ith box can accommodate all `n`
            balls and thus effectively has unlimited capacity.
      
      
            ACKNOWLEDGMENT
      
            I'd like to thank Chris Rebert for helping me to convert my prototype
            class-based code into a generator function.
            """
            def _unlabeled_balls_in_labeled_boxes(balls, box_sizes): #needed for 'unlabeled_balls_in_labeled_boxes' needed for 'states'
               """
               This recursive generator function was designed to be returned by
               `unlabeled_balls_in_labeled_boxes`.
               """
               # If there are no balls, all boxes must be empty:
               if not balls:
                  yield len(box_sizes) * (0,)
            
               elif len(box_sizes) == 1:
         
                  # If the single available box has sufficient capacity to store the balls,
                  # there is only one possible distribution, and we return it to the caller
                  # via `yield`.  Otherwise, the flow of control will pass to the end of the
                  # function, triggering a `StopIteration` exception.
                  if box_sizes[0] >= balls:
                     yield (balls,)
         
               else:
                  # Iterate over the number of balls in the first box (from the maximum
                  # possible down to zero), recursively invoking the generator to distribute
                  # the remaining balls among the remaining boxes.
                  for balls_in_first_box in range( min(balls, box_sizes[0]), -1, -1 ):
                     balls_in_other_boxes= balls - balls_in_first_box
                     for distribution_other in _unlabeled_balls_in_labeled_boxes(
                     balls_in_other_boxes, box_sizes[1:]):
                        yield (balls_in_first_box,) + distribution_other
               # end three alternative blocks
     
            
            #skip error-checks: balls is not integer, <0 or not a non-empty list
            #capacity= 0
            #for size in box_sizes:
            #   capacity+= size
            #
            #if capacity < balls:
            #   raise ValueError("The total capacity of the boxes is less than the "
            #   "number of balls to be distributed.")
      
            return _unlabeled_balls_in_labeled_boxes(balls, box_sizes)
            # end def unlabeled_balls_in_labeled_boxes(balls, box_sizes)
      
         #States=np.zeros((math.factorial(n+alpha-1)/(math.factorial(n)*
                                 #math.factorial(alpha-1)),alpha))
         States2=[]
         a=np.ones(alpha).tolist()
         for i in range(len(a)):
            a[i]=n #create the needed list
         for distributions in unlabeled_balls_in_labeled_boxes(n, a):
            States2.append(np.array(distributions, dtype=np.int8)) #save memory!
         return States2
   
      lines=[]
      freqs=[]
      initF=[]
      lineapp=lines.append
      freqapp=freqs.append
      initapp=initF.append
   
      L2=CalcI00()
      #both trees can be considered to coincide for first state. 
      L1=L2 
      for i in range(1,N+1):
         L1, L2=iterate(L1, L2, i)
         if L1==0 and L2==0:
            #only by assert: MemoryError
            break # kill calculation
         spect=L2.extract()
         for j in xrange(len(spect)):
            #I_nn is stored; the intensity is I_nn*I_nn
            lineapp(spect[j][0]*spect[j][0])  
            freqapp(spect[j][1])             #energy of respective transition 
            initapp(spect[j][2])             #for thermal population: frequency in init. state
      result=np.zeros((3, len(lines) ))
      #print freqs
      result[0]=freqs
      T=self.T*self.Hartree2cm_1
      for i in range(len(result[0])):
         #arbitrary but constant number for mode
         result[2][i]=42
         result[1][i]=lines[i]*np.exp(-initF[i]/T) #thermally weighting of transitions
      return result
   
class SDR_spect(Spect.Spect):
   """The class to calculate Duschinsky-spectra in a pseudo one-particle approximation. 
      This makes the computation much faster than in URDR-spect, but in fact the model behind
      is inconsistent. If the Duschinsky effect is large, this model is no more applicable.
   """
   m = 10
   Threshold=3e-7

   def __init__(self, f): 
      """In addition to the initialisation of
         Spect( see first line of code), here the number
         of modes is evaluated and the matrices/vectors
         needed for evaluation of the recursive equations are
         calculated.
      """
      # first, initialise as done for all spectra:
      self.type='URDR'
      Spect.Spect.__init__(self, f)
      #get m from opt. (-> number of states considered)
      modes=re.findall(r"(?<=modes\=)[\d]+",self.opt, re.I)
      if modes!=[]:
         self.m=modes[-1]
      else:
         self.m=len(self.nm.J)
   
   def __GetQuants(self):
      """This function computes the matrices and vectors needed interally for the computation
         of intensities.
      """
      self.Gamma=np.diag(self.f[0])   # in atomic units. It is equivalent to 4pi^2/h f_i
      self.Gammap=np.diag(self.f[1])  # for initial state
      sqGamma=np.diag(np.sqrt(self.f[0]))
      sqGammap=np.diag(np.sqrt(self.f[1]))
      unity=np.eye(len(self.Gamma))

      TMP=np.linalg.inv(self.nm.J.T.dot(self.Gammap).dot(self.nm.J) + self.Gamma)
      self.A=2.*sqGammap.dot(self.nm.J).dot(TMP).dot(self.nm.J.T).dot(sqGammap) -unity
      self.b=2.*sqGammap.dot( unity - self.nm.J.dot(TMP).dot(self.nm.J.T).dot(self.Gammap) ).dot(self.nm.K)
      self.C=2.*sqGamma.dot(TMP).dot(sqGammap) -unity
      self.d=-2.*sqGamma.dot(TMP).dot(self.nm.J.T.dot(self.Gammap.dot(self.nm.K)))
      self.E=4.*sqGamma.dot(TMP).dot(self.nm.J.T).dot(sqGammap)

   def __simpleFCfOPA(self, J, K, f, Energy, N, T, E0=0):
      """Calculates the FC-factors for given Duschinsky-effect. No restriction to OPA
      
         **PARAMETERS:**  
         J:      Duschisky-matrix  
         K:      Displacement-Vector  
         f:      frequency: two-dim array (freq_initial, freq_final)  
         Energy: Energy-difference of minima  
         N:      Max. number of excitation quanta state considered  
         T:      temperature of the system (in atomic units)  
         E0:     relative energy with respect to lowest triplet state    
         
         All arguments are obligatory.  
      
         **RETURNS:**  
         linespectrum 
      """
   
      def CalcI00( E):
         """This function calculates the overlap-integral for zero vibrations """
         opa=OPA(self.dim,0) #call OPA-class, initialize object for 0->0 trasition.
         zeros=np.zeros(2*self.dim)
         #set intensity of pure electronic transition (scaled arbitrarily)
         inten=1.00 
         # insert the transition
         opa.insert(zeros, np.sqrt(inten)) 
         return opa
   
      def iterate(L1, L2, i, f, J, K):
         """ Calculates the Franck-Condon factors of an eletronic transition using the lower levels L1 and L2   
            **PARAMETERS:**  
            L1:     binary tree where i-2 quanta are excited (structure: [[Btree.py]]  
            L2:     binary tree where i-1 quanta are excited  
            i:      number of excitation-quanta  
            f:      (2xN) frequencies of both states  
            J:      Duschisky-rotation matrix  
            K:      Displacement-vector  
   
            **RETURNS:**  
            L2:     input-parameter (needed for next iteration)  
            L3:     new binary tree   
         """
         def states( alpha, n): 
            """This function creates all possible states having a total number of n excitations in alpha different states  
               **PARAMETERS:**  
               alpha: number of degrees of freedom  
               n:     number of excitations  
   
               **RETURNS:**
            """
            States=[]
            distributions=np.zeros(2*alpha, dtype=np.int8)
            for i in range(alpha):
               for j in range(n+1):
                  distributions[i]=n-j
                  distributions[i+alpha]=j
                  States.append(np.array(distributions)) #save memory!
                  distributions[i]=0
                  distributions[i+alpha]=0
            return States
   
         #initialize new OPA-object
         leng=len(self.b)
         L3=OPA(leng, i)                   # initialize root-node
         States=states(leng, i)           # States are all possible
         npsqrt=np.sqrt
   
         for n in States: 
            m=np.argmax(n[:leng])  #if there is no excitation: it returns 0
            if n[m]!=0:
               # need first iteration-formula
               n_m=n[m]
               n[m]-=1 #n[m] is at least 1
               Ps=L2.getState(n)
               I_nn=self.b[m]*Ps                                     # first term 
               if n[m]>0:
                  n[m]-=1
                  Ps=L1.getState(n)

                  I_nn+=npsqrt(2.0*(n_m-1))*self.A[m][m]*Ps          # second term
                  n[m]+=1
               if n[m+leng]>0:
                  n[m+leng]-=1
                  Ps=L1.getState(n)
                  n[m+leng]+=1
                  I_nn+=npsqrt(n[m+leng]*0.5)*(self.E[m][m])*Ps # second term
   
               n[m]+=1
            #else: need the other iteration-formula
            else: 
               m=np.argmax(n[leng:])  #in this range must be an excitation.
               n_m=n[m+leng]
               n[m+leng]-=1
               Ps=L2.getState(n)
               I_nn=self.d[m]*Ps                                       # first term 
               if n[m+leng]>0:
                  n[m+leng]-=1
                  Ps=L1.getState(n)
                  I_nn+=npsqrt(2.0*(n_m-1.))*self.C[m][m]*Ps           # second term
                  n[m+leng]+=1
               n[m+leng]+=1
            I_nn/=npsqrt(2.*n_m)
            assert not math.isnan(I_nn) ,"I_nn is not a number! I_nn:"+\
                                         " {0}\n, n:{1}\n:".format(I_nn, n)
            # don't mess with too low intensities (<1e-16)
            if np.abs(I_nn)>1e-8: 
               L3.insert(n, I_nn)
         return L2, L3

      def makeLine(intens, E0, T, index, ex, Gamma, Gammap, E, n):
         """This function computes the actual spectrum .

            ===PARAMETERS===
            intens   a vector containing the square root of intensities of
                     the transitions
            E0       offset in energy compared to lowest electronic state
                     (always 0, will be removed in later versions)
            T        temperature of system, in a.u.
            index    contains the number of mode that changes in the 
                     respective transition (is a vector)
            ex       vector containing the number of quanta in the initial 
                     vibronic state
            Gamma    vector of frequencies in final state
            Gammap   vector of frequencies in initial state
            E        electronic energy gap (minimum to minimum)
            n        total number of vibrational quanta considered

            ===RETURNS===
            F   Matrix containing the vibronic transitions in the format
                energy    intensity   mode nr.
                Due to the model used, these transitions all are in OPA,
                the use of MultiPart is inconsistent but can be used to enhance
                the spectrum.
         """
         #F=np.zeros(( len(index),3 ))
         F=[]
         for i in xrange(len(index)):
            indi=index[i]
            if intens[i]*intens[i]*np.exp(-(Gammap[indi]*ex[i]+E0)/T) >self.Threshold:
               F.append([(-np.sign(E)*Gamma[indi]*(n-ex[i])
                         +np.sign(E)*Gammap[indi]*(ex[i])+np.abs(E))*self.Hartree2cm_1,
                        intens[i]*intens[i]*np.exp(-(Gammap[indi]*ex[i]+E0)/T) ,
                        indi+1])
              # if F[-1][1]>0.0001:
              #    print index[i], ex[i], n-ex[i], F[-1][1], F[-1][0]-E*self.Hartree2cm_1
   
         if F==[]: # no transitions with enough intensity occured.
            F=[0,0,7]
         return np.matrix(F)
      
      dimen=0
   
      # correct for vibrational groundstates:
      Energy+=(sum(f[0])-sum(f[1]))*.5
   
      linspect=[]

      #append the 0->0 transition to the linspect-array
      inten=1.00 
      linspect.append(np.matrix([abs(Energy)*self.Hartree2cm_1, inten, 7])) 

      L2=CalcI00(Energy)
      #this is already extracted to linspect (using side-effects)
      L1=L2 
      for i in range(1, N+1):
         L1, L2=iterate(L1, L2, i, f, J, K)
         intens, index, excitation, N=L2.extract()
         linspect.append(makeLine(intens, E0, T, index, excitation, f[1], f[0], Energy, i))
         dimen+=len(linspect[i])+1
      spect=np.zeros((dimen,3))
      spect[:len(linspect[0])]=linspect[0]
      dimen=len(linspect[0])
      for i in xrange(1,len(linspect)):
         spect[dimen:dimen+len(linspect[i])]=linspect[i]
         dimen+=len(linspect[i])
      return spect.T
   
   def calcspect(self):
      """This function prepares all variables according to the given 
         options to calculate the spectrum.
      """
      #first: resort the elements of J, K, f to make J most closely to unity
      resort=np.zeros(np.shape(self.nm.J))
      for i in xrange(len(self.nm.J)):
         j=np.argmax(self.nm.J[i])
         k=np.argmin(self.nm.J[i])
         if self.nm.J[i][j]>-self.nm.J[i][k]:
           resort[i][j]=1
         else:
            resort[i][k]=-1

      self.nm.J=resort.dot(self.nm.J.T)
      self.nm.K=resort.dot(self.nm.K.T)
      for i in xrange(len(resort)):
         k=np.argmin(resort[i])
         if resort[i][k]==-1:
            resort[i][k]=1 #use absolute value only.
      self.f[1]=resort.dot(self.f[1].T)
      self.log.write("before truncation:\n")
      self.log.write("Duschinsky-Matrix:\n")
      self.log.printMat(self.nm.J)
      self.log.write("Displacement vector:\n")
      self.log.printVec(self.nm.K)
   
      #truncate only, if this really is truncation.
      if self.m<len(self.nm.J): 
         index=np.argsort(np.abs(self.nm.K), kind="heapsort")[::-1]
         ind=index[:self.m]
         self.nm.K=self.nm.K[ind]
         self.nm.J=self.nm.J[ind].T[ind]
         f=np.zeros(( 2,self.m))
         f[1]=self.f[1][ind]
         f[0]=self.f[0][ind]
         self.f=f

      #print truncated quantites:
      self.log.write("After truncation:\n")
      self.log.write("Duschinsky-Matrix:\n")
      self.log.printMat(self.nm.J)
      self.log.write("Displacement vector:\n")
      self.log.printVec(self.nm.K)
   
      # finally, calculate the stick spectrum in this picture
      E=self.Energy[0]-self.Energy[1]
      self.__GetQuants()
      self.spect=self.__simpleFCfOPA(self.nm.J, self.nm.K, self.f, E, self.states1+ self.states2, self.T, 0)
   
class OPA:
   """ This class containes the functions and objects for the calculation of vibronic spectra 
       with the SDR_spect class.

       ==OBJECTS==
       L   number of excited states (in total: initial+final)
       mat matrix containing FC-factors: first index: number of mode being excited, second index is exc. number
   """
   Threshold=3e-10

   def __init__(self, alpha, L): #should it be __new__?
      """ initialises the class
         alpha: number of vibrational modes
         L:     number of states excited
      """
      assert alpha>0, "There must be at least one degree of freedom!"
      self.L=L
      # first index in mat:number of mode, second: excitation number
      self.mat=np.zeros((2*alpha,L+1)) #change elements to some float32 or so...
   
   def insert(self, N, FC):
      exc=int(N[len(N)//2:].max())
      if exc>0:
         index=np.matrix(N[len(N)//2:]).argmax()
         self.mat[index][exc]=FC
      else:
         index=np.matrix(N[:len(N)//2]).argmax()
         self.mat[index][0]=FC

   def getState(self, N): 
      exc=int(max(N[len(N)//2:]))
      if exc>0:
         index=np.matrix(N[len(N)//2:]).argmax()
      else:
         index=np.matrix(N[:len(N)//2]).argmax()
      return self.mat[index][exc]

   def extract(self): #extract all elements 
      selfmat=self.mat
      intens=[]
      ind=[]
      excs=[]
      intapp=intens.append
      indapp=ind.append
      excapp=excs.append
      for index in range(len(selfmat)):
         for exc in range(len(selfmat[0])):
            if selfmat[index][exc]>self.Threshold:
               intapp(selfmat[index][exc])
               indapp(index)
               excapp(exc)
            elif selfmat[index][exc]<-self.Threshold:
               intapp(selfmat[index][exc])
               indapp(index)
               excapp(exc)

      return intens, ind, excs, self.L #I need squares as intensities

version='1.0'
# End of DR_spects.py
