#!/usr/bin/python
import numpy as np, math, logging
from copy import deepcopy

Hartree2cm_1=219474.63 
Threshold=3e-10

class OPA:
   def __init__(self,alpha,L): #should it be __new__?
      """ initializes the class"""
      assert alpha>0, "There must be at least one degree of freedom!"
      self.mat=np.zeros((alpha,L+1)) #change elements to some float32 or so...
   
   def insert(self, N, FC):
      ni=N[len(N)//2:]
      nf=N[:len(N)//2]
      #exc=np.max(ni)
      exc=exc[np.argmax(ni)]
      if exc>0:
	 index=np.argmax(ni)
	 self.mat[index][exc]=FC
      else:
	 index=np.argmax(nf)
	 self.mat[index][0]=FC

   def getState(self, N): 
      ni=N[len(N)//2:]
      nf=N[:len(N)//2]
      exc=exc[np.argmax(ni)]
      if exc>0:
	 index=np.argmax(ni)
	 return self.mat[index][exc]
      else:
	 index=np.argmax(nf)
	 return self.mat[index][0]

   def extract(self): #extract all elements 
      intens=[]
      ind=[]
      excs=[]
      for index in range(len(self.mat)):
	 for exc in range(len(self.mat[0])):
	    if self.mat[index][exc]>Threshold:
   	       intens.append(self.mat[index][exc])
   	       ind.append(index)
   	       excs.append(exc)

	       #print "diff", self.mat[index][exc],self.mat[index][-exc] #this should, in principle, coincide --> does!
      return intens, ind, excs #I need squares as intensities

def FCfOPA(J, K, f, Energy, N, T,E0):
   """Calculates the FC-factors for given Duschinsky-effect. No restriction to OPA
   
   *PARAMETERS:*
   J:      Duschisky-matrix
   K:      Displacement-Vector
   f:      frequency: two-dim array (freq_initial, freq_final)
   Energy: Energy-difference of minima
   N:      Max. number of excitation quanta state considered
   T:      temperature of the system (in atomic units)
   E0:     relative energy with respect to lowest triplet state
     
   All arguments are obligatory.

   *RETURNS:*
   linespectrum 
   """
   def CalcI00(dim, E):
      """This function calculates the overlap-integral for zero vibrations """

      opa=OPA(dim,1) #is this clear to be the respective class?
      zeros=np.zeros(2*dim)
      for i in range(len(zeros)): #insert this transition into all values... important for recursion...
	 zeros[i]=1
	 opa.insert(zeros, 10) #sum(sum()) due to matrix
	 zeros[i]=0
      linspect.append([10, E*Hartree2cm_1]) 
      return opa

   def iterate(L1, L2, Energy, i, f, J, K):
      """ Calculates the Franck-Condon factors of an eletronic transition using the lower levels L1 and L2
   
      *PARAMETERS:*
      L1:     binary tree where i-2 quanta are excited (structure: [[Btree.py]]
      L2:     binary tree where i-1 quanta are excited
      Energy: Energy-difference between the states (minimal energy)
      i:      number of excitation-quanta
      f:      (2xN) frequencies of both states
      J:      Duschisky-rotation matrix
      K:      Displacement-vector

      *RETURNS:*
      L2:     input-parameter (needed for next iteration)
      L3:     new binary tree 
      """
      def states(alpha, n): 
	 """This function creates all possible states having a total number of n excitations in alpha different states
      
       	 *PARAMETERS:*
	 alpha: number of degrees of freedom
	 n:     number of excitations
      
       	 *RETURNS:*
	 """
	 States=[]
	 distributions=np.zeros(2*alpha)
	 for i in range(alpha):
	    for j in range(n+1):
	       distributions[i]=n-j
	       distributions[i+alpha]=j
	       States.append(np.array(distributions, dtype=np.int8)) #save memory!
	       distributions[i]=0
	       distributions[i+alpha]=0
	 return States

      #quantities for the iterative spectrum-calculation
      Gamma=np.diag(f[0])              	# in atomic units. It is equivalent to 4pi^2/h f_i
      Gammap=np.diag(f[1])             	# for final state
      sqGamma=np.diag(np.sqrt(f[0]))   
      sqGammap=np.diag(np.sqrt(f[1]))  
      unity=np.eye(len(Gamma))
   
      C=np.linalg.inv(J.T.dot(Gammap).dot(J)+Gamma) #C is only temporary matrix here
      A=J.dot(np.dot(C,J.T)) 
      A=2*np.dot(sqGammap,A.dot(sqGammap))-unity
      TMP=J.dot(C).dot(J.T).dot(Gammap)
      b=2*sqGammap.dot(unity-TMP).dot(K)
      d=-2*sqGamma.dot(C).dot(J.T).dot(Gammap).dot(K)
      E=4*sqGamma.dot(C).dot(J.T).dot(sqGammap)
      C=2*sqGamma.dot(C).dot(sqGamma)-unity 		#this is 'real' C-matrix
   
      #initialize new OPA-object
      alpha=len(b)
      L3=OPA(alpha,i)    	  	# initialize root-node
      States=states(alpha, i) 		# States are all possible
      for n in States:			#for each possible state, described by n(vector)
	 # index of excited elements
	 m=np.argmax(n[:len(n)//2])	#if there is no excitation: it returns 0
	 I_nn=0
	 #need first iteration formula
	 if n[m]!=0:
	    n_m=n[m]
	    ntemp=deepcopy(n)
	    ntemp[m]-=1 #n[m] is at least 1
	    Ps=L2.getState(ntemp)
	    if not math.isnan(Ps) and abs(Ps)>1e-8:
	       I_nn=b[m]*Ps					# first term 
	    if ntemp[m]>0:
	       ntemp[m]-=1
	       Ps=L1.getState(ntemp)
	       if abs(Ps)>1e-8 and not math.isnan(Ps):
		  I_nn+=np.sqrt(2*(n_m-1))*A[m][m]*Ps		# second term
	    if n[m+len(n)//2]>0:
	       ntemp=deepcopy(n)
	       ntemp[m]-=1
	       ntemp[m+len(n)//2]-=1
     	       Ps=L1.getState(ntemp)
	       if not math.isnan(Ps) and abs(Ps)>1e-8:
		  I_nn+=np.sqrt(n[m+len(n)//2]*0.5)*E[m][m]*Ps	# second term
	 #else: need the other iteration-formula
	 else: 
	    m=np.argmax(n[len(n)//2:])				# index of excited elements
	    n_m=n[m+len(n)//2]
	    ntemp=deepcopy(n)
	    ntemp[m+len(n)//2]-=1
	    Ps=L2.getState(ntemp)
	    if not math.isnan(Ps) and abs(Ps)>1e-8:
	       I_nn=d[m]*Ps					# first term 
	    if ntemp[m+len(n)//2]>0:
	       ntemp[m+len(n)//2]-=1
	       Ps=L1.getState(ntemp)
	       if not math.isnan(Ps) and abs(Ps)>1e-8:
		  I_nn+=np.sqrt(2*(n_m-1))*C[m][m]*Ps        	# second term; all other terms vanish in OPA...
     	 I_nn/=np.sqrt(2*n_m)
	 L3.insert(n, I_nn)
      return L2, L3

   def makeLine(intens,E0, T, index, ex, Gamma, Gammap,E, n):
      F=np.zeros(( len(index),2 ))
      log=open("calculation.log", "a")
      for i in range(len(index)): 
	 F[i][1]=(Gammap[index[i]][index[i]]*(ex[i]+0.5)-
   		   Gamma[index[i]][index[i]]*(n-ex[i]+0.5)+E)*Hartree2cm_1
	 F[i][0]=intens[i]*intens[i]*np.exp((-Gamma[index[i]][index[i]]*ex[i]+E0)/T)
	 log.write(u"{0}".format(repr(F[i][1])+"  "+repr(F[i][0])+"  "+repr(index[i])))
	 #logging.critical(repr(F[i][1])+"  "+repr(F[i][0])+"  "+repr(index[i]))
      return np.array(F)
   
   np.shape(f)
   Gamma=np.diag(f[0])
   Gammap=np.diag(f[1]) # for final state

   linspect=[]
   logging.critical("frequency,           intensity ,      mode")
   L2=CalcI00(len(K), Energy)
   #this is already extracted to linspect (using side-effects)
   L1=L2 
   for i in range(1,N+1):
      L1, L2=iterate(L1, L2, Energy, i, f, J, K)
      intens, index, excitation=L2.extract()
      for lines in makeLine(intens,E0, T, index, excitation, Gamma, Gammap, Energy, i):
	 linspect.append(lines)
   spect=np.zeros(( 3,len(linspect) ))
   for i in range(len(linspect)):
      for j in range(2):
         spect[1-j][i]=linspect[i][j]
      spect[2][i]=2 ######################## this has to be changed!
   return spect 

version=0.1
