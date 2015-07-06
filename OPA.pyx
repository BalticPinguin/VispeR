#!/usr/bin/python
import math
import numpy as np
from copy import deepcopy

Hartree2cm_1=219474.63 
Threshold=3e-10

class OPA:
   """ This class containes the functions and objects for the calculation of vibronic spectra 
       in Duschinski-rotated systems.
   """
   def __init__(self, int alpha, int L): #should it be __new__?
      """ initializes the class
      alpha: number of vibrational modes
      L:     number of states excited
      """
      assert alpha>0, "There must be at least one degree of freedom!"
      self.L=L
      self.mat=np.zeros((2*alpha,L+1)) #change elements to some float32 or so...
   
   def insert(self, N, double FC):
      cdef int exc=int(N[len(N)//2:].max())
      cdef int index
      if exc>0:
         index=np.matrix(N[len(N)//2:]).argmax()
         self.mat[index][exc]=FC
      else:
         index=np.matrix(N[:len(N)//2]).argmax()
         self.mat[index][0]=FC

   def getState(self, double[:] N): 
      cdef int exc=int(max(N[len(N)//2:]))
      cdef int index
      if exc>0:
         index=np.matrix(N[len(N)//2:]).argmax()
      else:
         index=np.matrix(N[:len(N)//2]).argmax()
      return self.mat[index][exc]

   def extract(self): #extract all elements 
      cdef int index
      cdef int exc
      selfmat=self.mat
      intens=[]
      ind=[]
      excs=[]
      intapp=intens.append
      indapp=ind.append
      excapp=excs.append
      for index in range(len(selfmat)):
         for exc in range(len(selfmat[0])):
            if selfmat[index][exc]>Threshold:
               intapp(selfmat[index][exc])
               indapp(index)
               excapp(exc)
            elif selfmat[index][exc]<-Threshold:
               intapp(selfmat[index][exc])
               indapp(index)
               excapp(exc)

      return intens, ind, excs, self.L #I need squares as intensities

#def simpleFCfOPA(logging, double[:,:] J, double[:] K,double[:,:] f,double Energy, int N, float T, float E0):
def simpleFCfOPA(logging, J, K, double[:,:] f, double Energy, int N, float T, float E0):
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

   def CalcI00(int dim, float E, Gamma, Gammap,J):
      """This function calculates the overlap-integral for zero vibrations """
      cdef int i
      opa=OPA(dim,0) #is this clear to be the respective class?
      zeros=np.zeros(2*dim)
      #explicit calculation of this intensity is important for this part...
      #unity=np.eye(dim)
      #inten=Gammap.dot(K)
      #inten=(Gammap.dot(J).dot(np.linalg.inv(J.T.dot(Gammap.dot(J))+Gamma).dot(J.T)-unity)).dot(inten)
      #inten=np.exp(K.T.dot(inten))
      #inten*=np.linalg.det(Gamma)
      #inten/=np.linalg.det(J.dot(J.T.dot(Gammap.dot(J))+Gamma))
      #inten*=2^dim*100
      inten=10
      opa.insert(zeros, np.sqrt(inten)) #sum(sum()) due to matrix

      linspect.append(np.matrix([E*Hartree2cm_1, inten, 0]))
      return opa

   def iterate(L1, L2, float Energy, int i, f, J, K):
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
      def states(int alpha, int n): 
         """This function creates all possible states having a total number of n excitations in alpha different states
         *PARAMETERS:*
         alpha: number of degrees of freedom
         n:     number of excitations

         *RETURNS:*
         """
         cdef int i
         cdef int j
         States=[]
         cdef double[:] distributions=np.zeros(2*alpha)
         for i in range(alpha):
            for j in range(n+1):
               distributions[i]=n-j
               distributions[i+alpha]=j
               States.append(np.array(distributions)) #save memory!
               distributions[i]=0
               distributions[i+alpha]=0
         return States

      cdef int leng
      #quantities for the iterative spectrum-calculation
      Gamma=np.diag(f[1])               # in atomic units. It is equivalent to 4pi^2/h f_i
      Gammap=np.diag(f[0])              # for final state
      sqGamma=np.diag(np.sqrt(f[1]))   
      sqGammap=np.diag(np.sqrt(f[0]))
      unity=np.eye(len(Gamma))
   
      C=np.linalg.inv(J.T.dot(Gammap).dot(J)+Gamma) #C is only temporary matrix here
      A=J.dot(np.dot(C,J.T))
      A=2*np.dot(sqGammap,A.dot(sqGammap))-unity
      TMP=J.dot(C).dot(J.T).dot(Gammap)
      b=2*(sqGammap.dot((unity-TMP).dot(K)))
      d=-2*sqGamma.dot(C.dot(J.T.dot(Gammap.dot(K))))
      E=4*sqGamma.dot(C).dot(J.T).dot(sqGammap)
      C=2*sqGamma.dot(C).dot(sqGamma)-unity             #this is 'real' C-matrix
      #print "J", J, Gamma
      #print "matrices:"
      #print A, '\n', b, '\n', C, '\n', d, "\n", E

      #initialize new OPA-object
      leng=len(b)
      L3=OPA(leng, i)                   # initialize root-node
      States=states(leng, i)           # States are all possible
      npsqrt=np.sqrt

      for n in States: #for each possible state, described by n(vector)
         I_nn=0
         #need first iteration formula
         m=np.argmax(n[:leng])  #if there is no excitation: it returns 0
         # if there excists excitation in initial state: need first iteration formula
         if n[m]!=0:
            # need first iteration-formula
            n_m=n[m]
            n[m]-=1 #n[m] is at least 1
            Ps=L2.getState(n)
            #if not math.isnan(Ps):
            I_nn=b[m]*Ps                                     # first term 
            if n[m]>0:
               n[m]-=1
               Ps=L1.getState(n)
               #if not math.isnan(Ps) and abs(Ps)>1e-8:
               I_nn+=npsqrt(2.0*(n_m-1))*A[m][m]*Ps         # second termA
               n[m]+=1
            if n[m+leng]>0:
               n[m+leng]-=1
               Ps=L1.getState(n)
               n[m+leng]+=1
               #if not math.isnan(Ps) and abs(Ps)>1e-8:
               I_nn+=npsqrt(n[m+leng]*0.5)*(E[m][m])*Ps # second term

            n[m]+=1
         #else: need the other iteration-formula
         else: 
            m=np.argmax(n[leng:])  #if there is no excitation: it returns 0
            n_m=n[m+leng]
            n[m+leng]-=1
            Ps=L2.getState(n)
            I_nn=d[m]*Ps                                    # first term 
            if n[m+leng]>0:
               n[m+leng]-=1
               Ps=L1.getState(n)
               #if not math.isnan(Ps) and abs(Ps)>1e-8:
               I_nn+=npsqrt(2.0*(n_m-1))*C[m][m]*Ps         # second term
               n[m+leng]+=1
            n[m+leng]+=1
         I_nn/=npsqrt(2*n_m)
         assert not math.isnan(I_nn) ,"I_nn is not a number! I_nn: {0}\n, n:{1}\n:".format(I_nn, n)
         L3.insert(n, I_nn)
         print n, I_nn
      return L2, L3

   #def makeLine(logging, double[:] intens, double E0,float T, int[:] index, int[:] ex, double[:,:] Gamma, double[:,:] Gammap,float E,int  n):
   def makeLine(logging, intens, double E0, double T, index, ex, Gamma, Gammap, float E, int n):
      cdef int indi
      F=np.zeros(( len(index),3 ))
      for i in xrange(len(index)):
         indi=index[i]
         F[i][2]= indi
         F[i][0]=(Gammap[indi][indi]*(ex[i]+0.5)-
                   Gamma[indi][indi]*(n-ex[i]+0.5)+E)*Hartree2cm_1
         F[i][1]=intens[i]*intens[i]*np.exp(-(Gamma[indi][indi]*ex[i]+E0)/T) 
   #      logging[1].write(u"%f   %f    %f\n"%(F[i][1], F[i][0], indi))
      return F
   
   cdef int dimen=0
   cdef int i

   cdef double[:,:] Gamma=np.diag(f[0])
   cdef double[:,:] Gammap=np.diag(f[1]) # for final state

   linspect=[]
   #logging[1].write("frequency,           intensity ,      mode\n")
   L2=CalcI00(len(K), Energy, Gamma, Gammap, J)
   #this is already extracted to linspect (using side-effects)
   L1=L2 
   for i in range(1, N+1):
      L1, L2=iterate(L1, L2, Energy, i, f, J, K)
      intens, index, excitation, N=L2.extract()
      linspect.append(makeLine(logging, intens, E0, T, index, excitation, Gamma, Gammap, Energy, i))
   for i in xrange(len(linspect)):
      dimen+=len(linspect[i])
   spect=np.zeros((dimen,3))
   spect[:len(linspect[0])]=linspect[0]
   dimen=len(linspect[0])
   for i in xrange(1,len(linspect)):
      spect[dimen:dimen+len(linspect[i])]=linspect[i]
      dimen+=len(linspect[i])
   return spect.T

def resortFCfOPA(logging, J, K, f, Energy,int N,float T,float E0):
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
   resort=np.zeros(np.shape(J))
   for i in range(len(J)):
      j=np.argmax(J[i])
      k=np.argmin(J[i])
      if J[i][j]>-J[i][k]:
         resort[i][j]=1
      else:
         resort[i][k]=-1
   J=resort.dot(J)
   #J=np.identity(len(J))
   K=resort.dot(K)
   for i in range(len(resort)):
      k=np.argmin(resort[i])
      if resort[i][k]==-1:
         resort[i][k]=1 #use absolute value only.
   f[1]=resort.dot(f[1])
   spect=simpleFCfOPA(logging, J, K, f, Energy, N, T,E0)
   return spect

def distFCfOPA(logging, J, K, f,float Energy,int N,float T,float E0,int threshold):
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

   resort=np.zeros(np.shape(J))
   for i in range(len(J)):
      j=np.argmax(J[i])
      k=np.argmin(J[i])
      if J[i][j]>-J[i][k]:
         resort[i][j]=1
      else:
         resort[i][k]=-1
   #J=resort.dot(J.T)
   J=resort
   K=resort.dot(K.T)

   for i in range(len(resort)):
      k=np.argmin(resort[i])
      if resort[i][k]==-1:
         resort[i][k]=1 #use absolute value only.
   f[1]=resort.dot(f[1].T)
   spect2=[]
   spect2.append(simpleFCfOPA(logging, J, K, f, Energy, N, T,E0))
   #print simpleFCfOPA(logging, J, K, f, Energy, N, T,E0).T

   resort=np.eye(len(resort))
   if threshold>len(resort)/2:
      #which is due to full OPA-scheme
      threshold=len(resort)/2

   Jo=deepcopy(J)
   Ka=deepcopy(K)
   f1=deepcopy(f[1])
   f0=deepcopy(f[0])
   for i in range(threshold):
      vec=deepcopy(resort[0])
      for j in range(len(resort)-1):
         resort[j]=resort[j+1]
      resort[-1]=vec

      J=resort.dot(Jo.T)
      K=resort.dot(Ka.T)
      f[1]=resort.dot(f1.T)
      f[0]=resort.dot(f0.T)
      #print "J", J,"\n"
      temp=simpleFCfOPA(logging, J, K, f, Energy, N, T,E0)
      spect2.append(temp[:].T[1:].T)# do not take the 0-0 transititon a second time
      #print (simpleFCfOPA(logging, J, K, f, Energy, N, T,E0)[:].T[1:])

   resort=np.eye(len(resort))
   for i in range(threshold):
      vec=deepcopy(resort[-1])
      for j in range(len(resort)-1,-1,-1):
         resort[j]=resort[j-1]
      resort[0]=vec

      J=resort.dot(Jo.T)
      K=resort.dot(Ka.T)
      f[1]=resort.dot(f1.T)
      f[0]=resort.dot(f0.T)
      #print "J", J,"\n"
      temp=simpleFCfOPA(logging, J, K, f, Energy, N, T,E0)
      spect2.append(temp[:].T[1:].T)# do not take the 0-0 transititon a second time
      #spect2.append(temp[1:])# do not take the 0-0 transititon a second time
      #print (simpleFCfOPA(logging, J, K, f, Energy, N, T,E0)[:].T[1:])

   dim=0
   for i in range(len(spect2)):
      dim+=len(spect2[i][0])
   spect=np.zeros((dim, 3))
   leng=0
   for i in range(len(spect2)):
      spect[leng:leng+len(spect2[i][0])]=spect2[i].T
      leng+=len(spect2[i][0])

   return spect.T

version=0.5
