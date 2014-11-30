#!/usr/bin/python
# filename: Dusch_unrest.py
import numpy as np, math
from copy import deepcopy
#include [[Btree.py]]
import Btree as bt

Hartree2cm_1=219474.63 
Threshold=3e-10

def unrestricted(logging, J, K, F, Energy, N, T, E0, m):
   """This function 
      **PARAMETERS**
      logging: object having as first element the mode and as second the (already opened) log-file
      J:	Duschisky-matrix
      K:	Displacement-Vector
      f:	frequency: two-dim array (freq_initial, freq_final)
      Energy:	Energy-difference of minima
      N:	Max. number of excitation quanta state considered
      T:	Temperature of the system
      E0:	frequency of the purely electronic transition
      m:	number of vibrational modes taken into account for the spectum calculation

      **RETURNS**
      linespectrum: vibrational line-spectrum
   """
   #first: resort the elements of J, K, f to make J most closely to unity
   resort=np.zeros(np.shape(J))
   for i in range(len(J)):
      j=np.argmax(J[i])
      k=np.argmin(J[i])
      if J[i][j]>-J[i][k]:
	 resort[i][j]=1
      else:
	 resort[i][k]=-1

   J=resort.dot(J.T)
   K=resort.dot(K.T)
   for i in range(len(resort)):
      k=np.argmin(resort[i])
      if resort[i][k]==-1:
	 resort[i][k]=1 #use absolute value only.
   F[1]=resort.dot(F[1].T)

   ## change the following: not size of K but off-diagonal-elements of J!!! 
   ## if J[i][j] is large, add indices i and j to ind if they are not contained already...
   Jtemp=np.zeros((len(J)-1)*(len(J)-1))
   Jtemp[:len(J)-1]=J[0][1:]
   for i in range(1,len(Jtemp)):
      for j in range(i):
	 Jtemp[i*len(Jtemp)+j]=J[i][i+j]
      for j in range(len(Jtemp)-i):
	 Jtemp[i*len(Jtemp)+j+i]=J[i][i+j+1]
   ##now: find biggest elements in Jtemp: get indices of it -> ready.
   index=np.argsort(np.abs(Jtemp), kind="heapsort")
   print index, Jtemp[index]

   # index K by the size of its elements and truncate K and J due to this.
   index=np.argsort(np.abs(K), kind="heapsort")
   if m<len(index):
      ind=index[:m]
   else:
      ind=index
   k=K[ind]
   j=J[ind].T[ind]
   f=np.zeros(( 2,len(ind) ))
   f[1]=F[1][ind]
   f[0]=F[0][ind]
   
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
   T:	   Temperature of the system
   E0:	   frequency of the purely electronic transition
     
   All parameters are obligatory.

   *RETURNS:*
   linespectrum 
   """
   def CalcI00(J, K, Gamma, Gammap, E):
      """This function calculates the overlap-integral for zero vibrations """
      #pref=math.pow(2,len(Gamma))*np.linalg.det(Gamma)
      #TMP=J.dot(J.T.dot(Gammap).dot(J)+Gamma)
      #pref/=np.linalg.det(TMP)

      ########error!!!! This has to be positive but is not always!!
      #pref=np.sqrt(np.abs(pref))
      #TMP=J.T.dot(Gammap).dot(J)+Gamma
      #TMP=Gammap.dot(J).dot(np.linalg.inv(TMP)).dot(J.T)-np.eye(len(J))
      #exp=np.exp(0.5*K.T.dot(TMP).dot(Gammap).dot(K))

      Tree=bt.Tree(2*len(K))
      Tree.fill(0)
      Zero=np.zeros(2*len(K))
      #Tree.insert(Zero, [pref*exp, (E+sum(sum(Gammap-Gamma))/2)*Hartree2cm_1] ) #sum(sum()) due to matrix
      Tree.insert(Zero, [10, (E+sum(sum(Gammap-Gamma))/2)*Hartree2cm_1] ) #sum(sum()) due to matrix
      #I_00 transition-probability [[Btree.py#extract]]
      linespect=np.array(Tree.extract())
      ### this is done using implicit side effects
      lines.append(linespect[0][0])
      freqs.append(linespect[0][1])
      return Tree
   
   def iterate(L1, L2, Energy, i, f, J, K):
      """ Calculates the Franck-Condon factors of an eletronic transition using the lower levels L1 and L2
   
      *PARAMETERS:*
      L1:     binary tree where i-2 quanta are excited (structure: [[Btree.py]]
      L2:     binary tree where i-1 quanta are excited
      Energy: Energy-difference between the states (minimal energy)
      i:      number of excitation-quanta
      f:      (2xN) frequencies of both states
      J:	   Duschisky-rotation matrix
      K:	   Displacement-vector

      *RETURNS:*
      L2:     input-parameter (needed for next iteration)
      L3:     new binary tree 
      """
   
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
      b=2*(sqGammap.dot((unity-TMP).dot(K)))
      d=-2*sqGamma.dot(C.dot(J.T.dot(Gammap.dot(K))))
      E=4*sqGamma.dot(C).dot(J.T).dot(sqGammap)
      C=2*sqGamma.dot(C).dot(sqGamma)-unity 		#this is 'real' C-matrix

      #initialize new tree
      alpha=2*len(b)
      L3=bt.Tree(i)    	       		# initialize root-node
      L3.fill(alpha)         		# initialize tree
      States=states(alpha, i) 		# States are all possible

      def freq(E, Gamma, Gammap):
         """Calculates the frequency of respective transition including vibrational frequency

	 *PARAMETERS:*
	 E:	 energy-difference of states
	 Gamma:	 vector of vibrational frequencies in inital state (in atomic units)
	 Gammap: vector of vibrational frequencies in final state (in atomic units)

	 *RETURNS;*
	 frequency of respective transition
	 """
	 return (E+sum(Gammap-Gamma))*Hartree2cm_1 
   
      def FirstNonzero(n): 
	 """Find first non-zero elements in first and second half of array n """
	 leng=len(n)//2
	 ni=n[:leng] #interger division (python3-compatible)
	 nf=n[leng:]
	 m=leng+1 #this means there is no excitation in this state
	 mp=leng+1
	 for j in range(leng):
	    if ni[j]>0:
	       m=j
	       break
	 for j in range(leng):
	    if nf[j]>0:
	       mp=j
	       break
	 return m, mp

      for n in States: #for each possible state, described by n(vector)
	 m, mp= FirstNonzero(n)# index of first-non-zero element of (initial, final) state
	 # if the 'first' excited state is in initial state: need first iteration formula
	 I_nn=0
	 leng=len(n)//2
	 if m<=mp:
	    # need first iteration-formula
	    n_m=n[m]
	    ntemp=deepcopy(n)
	    ntemp[m]-=1 #n[m] is at least 1
	    Ps=L2.getState(ntemp)[0]
	    if not math.isnan(Ps) and abs(Ps)>1e-8:
	       I_nn=b[m]*Ps					# first term 
	    if ntemp[m]>0:
	       ntemp[m]-=1
	       Ps=L1.getState(ntemp)[0]
	       if not math.isnan(Ps) and abs(Ps)>1e-8:
		  I_nn+=np.sqrt(2*(n_m-1))*A[m][m]*Ps		# second term
	    for i in range(m+1, leng):
	       if n[i]>0:
		  ntemp=deepcopy(n)
		  ntemp[m]-=1
		  ntemp[i]-=1
		  Ps=L1.getState(ntemp)[0]
		  if not math.isnan(Ps) and abs(Ps)>1e-8:
		     I_nn+=np.sqrt(float(n[i])*0.5)*(A[m][i]+A[i][m])*Ps	# second term

	    for i in range(mp+leng, len(n)): 			# sum over respective final states
	       if mp>len(n)//2:					# that means: there are no excited vibrations
		  break
	       if n[i]>0:
		  ntemp=deepcopy(n)
		  ntemp[m]-=1
		  ntemp[i]-=1
		  Ps=L1.getState(ntemp)[0]
		  if not math.isnan(Ps) and abs(Ps)>1e-8:
		     I_nn+=np.sqrt(float(n[i])*0.5)*(E[i-leng][m])*Ps		# second term
	 #else: need the other iteration-formula
	 else: 
	    n_m=n[mp+leng]
	    ntemp=deepcopy(n)
	    ntemp[mp+leng]-=1
	    Ps=L2.getState(ntemp)[0]
	    if not math.isnan(Ps) and abs(Ps)>1e-8:
	       I_nn=d[mp]*Ps					# first term 
	    if ntemp[mp+leng]>0:
	       ntemp[mp+leng]-=1
	       Ps=L1.getState(ntemp)[0]
	       if not math.isnan(Ps) and abs(Ps)>1e-8:
		  I_nn+=np.sqrt(2*(n_m-1))*C[mp][mp]*Ps         # second term
	    for i in range(mp+1, len(n)):
	       if n[i]>0:
		  ntemp=deepcopy(n)
		  ntemp[mp+leng]-=1
		  ntemp[i]-=1
		  Ps=L1.getState(ntemp)[0]
		  if not math.isnan(Ps) and abs(Ps)>1e-8:
		     I_nn+=np.sqrt(n[i]/2)*(C[mp][i-len(n)//2]+ # second term
			      C[i-leng][mp])*Ps	
	    for i in range(m, leng): 				#sum over respective final states
	       if m>len(n)//2:					# that means: there are no excited vibrations
		  break						#actually not needed, right?
   	       if n[i]>0:
   		  ntemp=deepcopy(n)
   		  ntemp[mp+leng]-=1
   		  ntemp[i]-=1
   		  Ps=L1.getState(ntemp)[0]
   		  if not math.isnan(Ps) and abs(Ps)>1e-8:
   		     I_nn+=np.sqrt(n[i]/2)*(E[mp][i-len(n)//2])*Ps # second term
     	 I_nn/=np.sqrt(2*n_m)
	 #threshold for insertion: saves memory, since int insead of float is used
	 if np.abs(I_nn)>1e-8:
	    try:
	       L3.insert(n, [I_nn, freq(Energy, f[0]*n[:leng], f[1]*n[leng:]) ])
	       print n, I_nn*I_nn, freq(Energy, f[0]*n[:leng], f[1]*n[leng:])
	    except MemoryError: 
	       logging[1].write('memory-error by inserting data. Finishing calculation.')
	       logging[1].writelines("%s\n" % item  for item in L2.extract())
	       #linspect=open('/tmp/linspect', "a")
	       #linspect.writelines("%s\n" % item  for item in L2.extract())
	       #linspect.close()
	       return 0,0
      return L2, L3

   def states(alpha, n): 
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


   Gamma=np.diag(f[0]) #in atomic units. It is equivalent to 4pi^2/h f_i
   Gammap=np.diag(f[1]) # for final state
   lines=[]
   freqs=[]

   L2=CalcI00(J, K, Gamma, Gammap, Energy)
   #both trees can be expected to coincide for first state. 
   L1=L2 
   for i in range(1,N+1):
      L1, L2=iterate(L1, L2, Energy, i, f, J,K)
      if L1==0 and L2==0:
	 #only by assert: MemoryError
	 break # kill calculation
      spect=L2.extract()
      for j in range(len(spect)):
	 lines.append(spect[j][0])
	 freqs.append(spect[j][1])
   result=np.zeros((3, len(lines) ))
   result[0]=freqs
   for i in range(len(result[0])):
      #arbitrary but constant number for mode
      result[2][i]=42
      # intensities are proportional to square of the transition rates
      result[1][i]=lines[i]*lines[i] 
   return result

version=0.2
# End of Dusch_unrest.py
