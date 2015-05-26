#!/usr/bin/python
# filename: Dusch_unrest.pyx
import numpy as np
import math
cimport numpy as np
#include [[Btree.pyx]]
import Btree as bt

Hartree2cm_1=219474.63 
Threshold=3e-10

def unrestricted(logging, J, K, F, Energy, N, T, E0, m):
   """This function 
   **PARAMETERS**
   logging: object having as first element the mode and as second the (already opened) log-file
   J:        Duschisky-matrix
   K:        Displacement-Vector
   f:        frequency: two-dim array (freq_initial, freq_final)
   Energy:   Energy-difference of minima
   N:        Max. number of excitation quanta state considered
   T:        Temperature of the system
   E0:       frequency of the purely electronic transition
   m:        number of vibrational modes taken into account for the spectum calculation

   **RETURNS**
   linespectrum: vibrational line-spectrum
   """
   #first: resort the elements of J, K, f to make J most closely to unity
   resort=np.zeros(np.shape(J))
   for i in xrange(len(J)):
      j=np.argmax(J[i])
      k=np.argmin(J[i])
      if J[i][j]>-J[i][k]:
        resort[i][j]=1
      else:
         resort[i][k]=-1
   J=resort.dot(J.T)
   K=resort.dot(K.T)
   for i in xrange(len(resort)):
      k=np.argmin(resort[i])
      if resort[i][k]==-1:
         resort[i][k]=1 #use absolute value only.
   F[1]=resort.dot(F[1].T)

   # --> change the following: not size of K but off-diagonal-elements of J!!! 
   # --> if J[i][j] is large, add indices i and j to ind if they are not contained already...
   if m<len(J): #this model might be more realistic!!
      index=np.argsort(np.abs(K), kind="heapsort")
      ind=index[:m]
      k=K[ind]
      j=J[ind].T[ind]
      f=np.zeros(( 2,len(ind) ))
      f[1]=F[1][ind]
      f[0]=F[0][ind]
   else:
      k=K
      j=J
      f=np.zeros(( 2,len(F[0]) ))
      f[1]=F[1]
      f[0]=F[0]

   # finally, calculate the Duschinsky-rotated line spectrum in this picture
   linspect=FCf(logging, j, k, f, Energy, N, T, E0)
   return linspect #3-dimensional array

def FCf(logging, J, K, f, Energy, N, T, E0):
   """Calculates the FC-factors for given Duschinsky-effect. No restriction to OPA
   
   *PARAMETERS:*
   logging: object having as first element the mode and as second the (already opened) log-file
   J:      Duschisky-matrix
   K:      Displacement-Vector
   f:      frequency: two-dim array (freq_initial, freq_final)
   Energy: Energy-difference of minima
   N:      Max. number of excitation quanta state considered
   T:      Temperature of the system
   E0:     frequency of the purely electronic transition
     
   All parameters are obligatory.

   *RETURNS:*
   linespectrum 
   """
   def CalcI00(J, K, Gamma, Gammap, E):
      """This function calculates the overlap-integral for zero vibrations """

      Tree=bt.Tree(2*len(K))
      Tree.fill(0)
      Zero=np.zeros(2*len(K))
      Tree.insert(Zero, np.array([10, (E+sum(sum(Gammap-Gamma))*0.5)*Hartree2cm_1, 0]) ) #sum(sum()) due to matrix
      #I_00 transition-probability [[Btree.py#extract]]
      #linespect=np.array(Tree.extract())
      #this is done using implicit side effects
      lines.append(10)
      freqs.append((E+sum(sum(Gammap-Gamma))*0.5)*Hartree2cm_1)
      initF.append(0) #needed for boltzmann-weighing
      return Tree
   
   def iterate(L1, L2,float Energy, int i, f, J, K):
      """ Calculates the Franck-Condon factors of an eletronic transition using the lower levels L1 and L2
   
      *PARAMETERS:*
      L1:     binary tree where i-2 quanta are excited (structure: [[Btree.py]]
      L2:     binary tree where i-1 quanta are excited
      Energy: Energy-difference between the states (minimal energy)
      i:      number of excitation-quanta
      f:      (2xN) frequencies of both states
      J:           Duschisky-rotation matrix
      K:           Displacement-vector

      *RETURNS:*
      L2:     input-parameter (needed for next iteration)
      L3:     new binary tree 
      """
   
      #quantities for the iterative spectrum-calculation
      cdef int leng, m, mp, alpha
      cdef double I_nn, Ps
      #quantities for the iterative spectrum-calculation
      npdiag=np.diag
      Gamma=npdiag(f[1])               # in atomic units. It is equivalent to 4pi^2/h f_i
      Gammap=npdiag(f[0])              # for final state
      sqGamma=npdiag(np.sqrt(f[1]))   
      sqGammap=npdiag(np.sqrt(f[0]))
      unity=np.eye(len(Gamma))
   
      S=np.linalg.inv(J.T.dot(Gammap).dot(J)+Gamma)
      A=2*sqGammap.dot(J).dot(S).dot(J.T).dot(sqGammap)-unity
      C=2*sqGamma.dot(S).dot(sqGamma)-unity
      E=4*sqGamma.dot(S).dot(J.T).dot(sqGammap)
      b=2*( sqGammap.dot(K)-sqGammap.dot(J.dot(S.dot((J.T).dot(Gammap.dot(K))))) )
      d=-2*sqGamma.dot(S.dot((J.T).dot(Gammap.dot(K))))

      #initialize new tree
      alpha=2*len(b)
      L3=bt.Tree(i)
      L3.fill(alpha)
      States=states(alpha, i)           # States are all possible

      def freq(float E, Gamma, Gammap):
         """Calculates the frequency of respective transition including vibrational frequency

         *PARAMETERS:*
         E:      energy-difference of states
         Gamma:  vector of vibrational frequencies in inital state (in atomic units)
         Gammap: vector of vibrational frequencies in final state (in atomic units)

         *RETURNS;*
         frequency of respective transition
         """
         return (E+sum(Gammap-Gamma))*Hartree2cm_1 
   
      def FirstNonzero(n): 
         """Find first non-zero elements in first and second half of array n """
         cdef int leng=len(n)//2
         cdef int m=leng+1 #this means there is no excitation in this state
         cdef int mp=leng+1
         cdef int j
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
      
      npsqrt=np.sqrt
      for n in States: #for each possible state, described by n(vector)
         m, mp= FirstNonzero(n)# index of first-non-zero element of (initial, final) state
         # if the 'first' excited state is in initial state: need first iteration formula
         leng=len(n)//2
         if m<=mp:
            # need first iteration-formula
            n_m=n[m]
            n[m]-=1 #n[m] is at least 1
            Ps=L2.getState(n)[0]
            #if not math.isnan(Ps) and abs(Ps)>1e-8:
            I_nn=b[m]*Ps                                     # first term 
            if n[m]>0:
               n[m]-=1
               Ps=L1.getState(n)[0]
               I_nn+=npsqrt(2.0*(n_m-1))*A[m][m]*Ps           # second term
               n[m]+=1
            for i in range(m+1, leng):
               if n[i]>0:
                  n[i]-=1
                  Ps=L1.getState(n)[0]
                  n[i]+=1
                  I_nn+=npsqrt(n[i]*0.5)*(A[m][i]+A[i][m])*Ps   # second term
            for i in range(mp+leng, len(n)):                    # sum over respective final states
               if mp>leng:                                      # that means: there are no excited vibrations
                  break
               if n[i]>0:
                  n[i]-=1
                  Ps=L1.getState(n)[0]
                  n[i]+=1
                  I_nn+=npsqrt(n[i]*0.5)*E[i-leng][m]*Ps        # second term
            n[m]+=1
         #else: need the other iteration-formula
         else: 
            n_m=n[mp+leng]
            n[mp+leng]-=1
            Ps=L2.getState(n)[0]
            I_nn=d[mp]*Ps                                       # first term 
            if n[mp+leng]>0:
               n[mp+leng]-=1
               Ps=L1.getState(n)[0]
               I_nn+=npsqrt(2*(n_m-1.0))*C[mp][mp]*Ps           # second term
               n[mp+leng]+=1
            for i in range(mp+1+leng, len(n)):
               if n[i]>0:
                  n[i]-=1
                  Ps=L1.getState(n)[0]
                  n[i]+=1
                  I_nn+=npsqrt(n[i]*0.5)*(C[mp][i-leng]+        # second term
                              C[i-leng][mp])*Ps 
            for i in range(m, leng):                            #sum over respective final states
               if m>leng:                                       # that means: there are no excited vibrations
                  break                                         #actually not needed, right?
               if n[i]>0:
                  n[i]-=1
                  Ps=L1.getState(n)[0]
                  n[i]+=1
                  #if not math.isnan(Ps) and abs(Ps)>1e-8:
                  I_nn+=npsqrt(n[i]*0.5)*(E[mp][i-leng])*Ps # second term
            n[mp+leng]+=1
         I_nn/=npsqrt(2*n_m)
         assert not math.isnan(I_nn) ,"I_nn is not a number! I_nn: {0}\n, n:{1}\n:".format(I_nn, n)
         #threshold for insertion: saves memory, since int insead of float is used
         if np.abs(I_nn)>1e-8:
            L3.insert(n, [I_nn, freq(Energy, f[0]*n[:leng], f[1]*n[leng:]), freq(0, 0, f[1]*n[leng:]) ])
            #if I_nn>=10:
               #print n, I_nn, freq(Energy, f[0]*n[:leng], f[1]*n[leng:])
      return L2, L3

   def states(int alpha,int  n): 
      """This function creates all possible states having a total number of n excitations in alpha different states
   
      *PARAMETERS:*
      alpha: number of degrees of freedom
   
      *RETURNS:*
      """

      def unlabeled_balls_in_labeled_boxes(balls, box_sizes): #needed for 'states'
         """
         These functions are part of python-package: 'combinatorics' 
         (download from https://pypi.python.org/pypi/Combinatorics)
         unlabeled_balls_in_labeled_boxes(balls, box_sizes): This function 
         returns a generator that produces all distinct distributions of indistinguishable balls
         among labeled boxes with specified box sizes (capacities). This is 
         a generalization of the most common formulation of the problem, where each box is
         sufficiently large to accommodate all of the balls, and is an important 
         example of a class of combinatorics problems called 'weak composition' problems.
      
         OVERVIEW
      
         This function returns a generator that produces all distinct distributions of
         indistinguishable balls among labeled boxes with specified box sizes
         (capacities).  This is a generalization of the most common formulation of the
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
               for balls_in_first_box in xrange( min(balls, box_sizes[0]), -1, -1 ):
                  balls_in_other_boxes= balls - balls_in_first_box
                  for distribution_other in _unlabeled_balls_in_labeled_boxes(
                  balls_in_other_boxes, box_sizes[1:]):
                     yield (balls_in_first_box,) + distribution_other
            # end three alternative blocks
   
         if not isinstance(balls, int):
               raise TypeError("balls must be a non-negative integer.")
         if balls < 0:
            raise ValueError("balls must be a non-negative integer.")
      
         if not isinstance(box_sizes,list):
            raise ValueError("box_sizes must be a non-empty list.")
      
         capacity= 0
         for size in box_sizes:
            if not isinstance(size, int):
               raise TypeError("box_sizes must contain only positive integers.")
            if size < 1:
               raise ValueError("box_sizes must contain only positive integers.")
            capacity+= size
      
         if capacity < balls:
            raise ValueError("The total capacity of the boxes is less than the "
            "number of balls to be distributed.")
   
         return _unlabeled_balls_in_labeled_boxes(balls, box_sizes)
         # end def _unlabeled_balls_in_labeled_boxes(balls, box_sizes)
   
      #States=np.zeros((math.factorial(n+alpha-1)/(math.factorial(n)*math.factorial(alpha-1)),alpha))
      States2=[]
      a=np.ones(alpha).tolist()
      for i in range(len(a)):
         a[i]=n*int(a[i]) #create the needed list
      i=0
      for distributions in unlabeled_balls_in_labeled_boxes(n,a):
         #States[i]=np.matrix(distributions)
         try:
            States2.append(np.array(distributions, dtype=np.int8)) #save memory!
         except MemoryError: 
            #if the memory is not enough: don't add further states and use what is availibleA
            print('Memory is full for state',n,'having only', len(States2),'objects in it. Use these states only.')
            break #in principle: do more stuff here!
         i+=1
      return States2

   Gammap=np.diag(f[0]) # for initial state
   Gamma=np.diag(f[1]) #in atomic units. It is equivalent to 4pi^2/h f_i
   lines=[]
   freqs=[]
   initF=[]
   lineapp=lines.append
   freqapp=freqs.append
   initapp=initF.append

   L2=CalcI00(J, K, Gamma, Gammap, Energy)
   #both trees can be considered to coincide for first state. 
   L1=L2 
   for i in range(1,N+1):
      L1, L2=iterate(L1, L2, Energy, i, f, J, K)
      if L1==0 and L2==0:
         #only by assert: MemoryError
         break # kill calculation
      spect=L2.extract()
      for j in xrange(len(spect)):
         lineapp(spect[j][0])
         freqapp(spect[j][1])
         initapp(spect[j][2])
   result=np.zeros((3, len(lines) ))
   result[0]=freqs
   T*=Hartree2cm_1
   result[2]=[42 for i in range(len(result[2]))]
   #result[2]=map(42, result[1]) # this may be even faster 
   #result[2]=42 #does this work? how fast is it??
   result[1]=[lines[i]*lines[i]*np.exp(-initF[i]/T)/10 for i in range(len(result[1]))]  # devide by ten since I_00 is 10 not 100 (should be square)
   #for i in range(len(result[0])):
      ##arbitrary but constant number for mode
      #result[2][i]=42
      ## intensities are proportional to square of the transition rates
      #result[1][i]=lines[i]*lines[i]*np.exp(-initF[i]/T)/10  # devide by ten since I_00 is 10 not 100 (should be square)
   return result

version=0.5
# End of Dusch_unrest.py
