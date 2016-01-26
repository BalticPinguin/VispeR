#!/usr/bin/python2
import Spect
from copy import deepcopy

class URDR-spect(Spect):
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
      ##K*=np.sqrt(np.pi)/2. # --> do it in functions_smsc.Duschinsky now
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
   
      #truncate only, if this really is truncation.
      if m<len(J): 
         index=np.argsort(np.abs(K), kind="heapsort")[::-1]
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
   
      # finally, calculate the Duschinsky-rotated stick spectrum in this picture
      linspect=FCf(logging, j, k, f, Energy, N, T, E0)
      return linspect #3-dimensional array
   
   def FCf(logging, J, K, f, Energy, N, T, E0):
      """Calculates the FC-factors for given Duschinsky-effect. 
       No restriction to OPA
      
      *PARAMETERS:*
      logging: object having as first element the mode and as second the 
              (already opened) log-file
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
         #E+=(sum(Gammap)-sum(Gamma))*.5    # correct for vibrational groundstates
         #linespect=np.array(Tree.extract())
         Tree.insert(Zero, np.array([1.0, (E)*Hartree2cm_1, 0]) ) #sum(sum()) due to matrix
         #I_00 transition-probability [[Btree.py#extract]]
         #linespect=np.array(Tree.extract())
         #this is done using implicit side effects
         lines.append(1.0)
         freqs.append((E)*Hartree2cm_1)
         initF.append(0) #needed for boltzmann-weighing
         return Tree
      
      def iterate(L1, L2, Energy, i, f, J, K):
         """ Calculates the Franck-Condon factors of an eletronic transition 
             using the lower levels L1 and L2
      
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
            return (E+sum(Gammap-Gamma))*Hartree2cm_1 
      
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
         
         Gamma=np.diag(f[0])   # in atomic units. It is equivalent to 4pi^2/h f_i
         Gammap=np.diag(f[1])  # for initial state
         sqGamma=np.diag(np.sqrt(f[0]))
         sqGammap=np.diag(np.sqrt(f[1]))
         unity=np.eye(len(Gamma))
      
         TMP=np.linalg.inv(J.T.dot(Gammap).dot(J) + Gamma)
         A=2.*sqGammap.dot(J).dot(TMP).dot(J.T).dot(sqGammap) -unity
         b=2.*sqGammap.dot( unity - J.dot(TMP).dot(J.T).dot(Gammap) ).dot(K)
         C=2.*sqGamma.dot(TMP).dot(sqGammap) -unity
         d=-2.*sqGamma.dot(TMP).dot(J.T.dot(Gammap.dot(K)))
         E=4.*sqGamma.dot(TMP).dot(J.T).dot(sqGammap)
      
         #initialize new tree
         alpha=2*len(b)
         L3=bt.Tree(i)
         L3.fill(alpha)
         States=states(alpha, i)           # States are all possible
   
         npsqrt=np.sqrt
         leng=len(States[0])//2
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
               I_nn=b[m]*Ps                                     # first term 
               #print "1", b[m]*Ps 
               if n[m]>0:
                  n[m]-=1
                  Ps=L1.getState(n)[0]
                  I_nn+=npsqrt(2.0*(n_m-1))*A[m][m]*Ps          # second term
                  #print "2", npsqrt(2.0*(n_m-1))*A[m][m]*Ps
                  n[m]+=1
               if n[m+leng]>0:
                  n[m+leng]-=1
                  Ps=L1.getState(n)[0]
                  n[m+leng]+=1
                  I_nn+=npsqrt(n[m+leng]*0.5)*(E[m][m])*Ps # second term
                  #print "3", npsqrt(n[m+leng]*0.5)*(E[m][m])*Ps 
   
               n[m]+=1
            #else: need the other iteration-formula
            else: 
               n_m=n[mp+leng]
               n[mp+leng]-=1
               Ps=L2.getState(n)[0]
               I_nn=d[mp]*Ps                                       # first term 
               #print "1a", d[mp]*Ps
               if n[mp+leng]>0:
                  n[mp+leng]-=1
                  Ps=L1.getState(n)[0]
                  I_nn+=npsqrt(2.0*(n_m-1.))*C[mp][mp]*Ps           # second term
                  #print "2a", npsqrt(2.0*(n_m-1.))*C[mp][mp]*Ps
                  n[mp+leng]+=1
               n[mp+leng]+=1
            I_nn/=npsqrt(2.*n_m)
            #print "4", npsqrt(2.*n_m)
            assert not math.isnan(I_nn) ,"I_nn is not a number! I_nn:"+\
                                         " {0}\n, n:{1}\n:".format(I_nn, n)
            if np.abs(I_nn)>1e-8: # don't mess with too low intensities (<1e-16)
               L3.insert(n, [I_nn,\
                freq(Energy, np.sign(Energy)*f[1]*n[:leng],\
                             np.sign(Energy)*f[0]*n[leng:]),\
                         freq(0, 0, f[0]*n[leng:]) ])
            #print "   ", n
            #m=min(m,mp)
            #print m, n[m], n[m+leng], I_nn, I_nn*I_nn ,\
            #      (f[1].T.dot(n[:leng])-f[0].T.dot(n[leng:]))*Hartree2cm_1
         return L2, L3
      
      def states( alpha, n): 
         """This function creates all possible states having a total number of n 
            excitations in alpha different states
      
         *PARAMETERS:*
         alpha: number of degrees of freedom
         n:     number of excitations
      
         *RETURNS:*
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
      
         #States=np.zeros((math.factorial(n+alpha-1)/(math.factorial(n)*
                                 #math.factorial(alpha-1)),alpha))
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
   
      def states2( alpha, n): 
         """This function creates all possible states having a total number of n 
             excitations in alpha DIFFERENT states.
             Hence it is for OPA-spectrum; here, for debugging only.
         *PARAMETERS:*
         alpha: number of degrees of freedom
         n:     number of excitations
   
         *RETURNS:*
         """
         alpha=alpha//2
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
   
      Gammap=np.diag(f[0]) # for initial state
      Gamma=np.diag(f[1]) #in atomic units. It is equivalent to 4pi^2/h f_i
      lines=[]
      freqs=[]
      initF=[]
      lineapp=lines.append
      freqapp=freqs.append
      initapp=initF.append
   
      L2=CalcI00(J, K, f[1], f[0], Energy)
      #both trees can be considered to coincide for first state. 
      L1=L2 
      for i in range(1,N+1):
         L1, L2=iterate(L1, L2, Energy, i, f, J, K)
         if L1==0 and L2==0:
            #only by assert: MemoryError
            break # kill calculation
         spect=L2.extract()
         for j in xrange(len(spect)):
            #I_nn is stored; the intensity is I_nn*I_nn
            lineapp(spect[j][0]*spect[j][0])  
            freqapp(spect[j][1])             #energy of respective transition 
            initapp(spect[j][2])             #for thermal population: frequency in init. state
           # if lines[-1] > 0.0001:
           #    print lines[-1], freqs[-1]
      result=np.zeros((3, len(lines) ))
      result[0]=freqs
      T*=Hartree2cm_1
      for i in range(len(result[0])):
         #arbitrary but constant number for mode
         result[2][i]=42
         result[1][i]=lines[i]*np.exp(-initF[i]/T) #thermally weighting of transitions
      #print result[0][0]-Energy*Hartree2cm_1, result[1][0]
      return result
   
class SDR-spect(Spect):
   Hartree2cm_1=219474.63 
   Threshold=3e-10

   class OPA:
      """ This class containes the functions and objects for the calculation of vibronic spectra 
          in Duschinski-rotated systems.
          ==OBJECTS==
          L   number of excited states (in total: initial+final)
          mat matrix containing FC-factors: first index: number of mode being excited, second index is exc. number
      """
      def __init__(self, alpha, L): #should it be __new__?
         """ initializes the class
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
               if selfmat[index][exc]>Threshold:
                  intapp(selfmat[index][exc])
                  indapp(index)
                  excapp(exc)
               elif selfmat[index][exc]<-Threshold:
                  intapp(selfmat[index][exc])
                  indapp(index)
                  excapp(exc)
   
         return intens, ind, excs, self.L #I need squares as intensities
   
   def simpleFCfOPA(logging, J, K, f, Energy, N, T, E0):
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
         opa=OPA(dim,0) #call OPA-class, initialize object for 0->0 trasition.
         zeros=np.zeros(2*dim)
         inten=1.00 #set intensity (scaled arbitrarily
         opa.insert(zeros, np.sqrt(inten)) # insert the intensity (transition strength)
   
         linspect.append(np.matrix([E*Hartree2cm_1, inten, 0])) #append the 0->0 transition to the linspect-array
         return opa
   
      def iterate(L1, L2, Energy,i, f, J, K):
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
         def states( alpha, n): 
            """This function creates all possible states having a total number of n excitations in alpha different states
            *PARAMETERS:*
            alpha: number of degrees of freedom
            n:     number of excitations
   
            *RETURNS:*
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
   
         #K*=np.sqrt(np.pi)/2
         #Gamma=np.diag(f[1])*2*np.pi               # in atomic units. It is equivalent to 4pi^2/h f_i
         #Gammap=np.diag(f[0])*2*np.pi              # for final state
         #sqGamma=np.diag(np.sqrt(f[1]))*np.sqrt(2*np.pi)
         #sqGammap=np.diag(np.sqrt(f[0]))*np.sqrt(2*np.pi)
         Gamma=np.diag(f[0])               # in atomic units. It is equivalent to 4pi^2/h f_i
         Gammap=np.diag(f[1])              # for initial state
         sqGamma=np.diag(np.sqrt(f[0]))
         sqGammap=np.diag(np.sqrt(f[1]))
         unity=np.eye(len(Gamma))
      
      #   C=np.linalg.inv(J.T.dot(Gammap).dot(J)+Gamma) #C is only temporary matrix here
      #   A=J.dot(np.dot(C,J.T))
      #   A=2*np.dot(sqGammap,A.dot(sqGammap))-unity
      #   TMP=J.dot(C).dot(J.T).dot(Gammap)
      #   b=2*sqGammap.dot((unity-TMP).dot((K)))
      #   d=-2*sqGamma.dot(C.dot(J.T.dot(Gammap.dot(K))))
      #   E=4*sqGamma.dot(C).dot(J.T).dot(sqGammap)
         TMP=np.linalg.inv(J.T.dot(Gammap).dot(J) + Gamma)
         A=2.*sqGammap.dot(J).dot(TMP).dot(J.T).dot(sqGammap) -unity
         b=2.*sqGammap.dot( unity - J.dot(TMP).dot(J.T).dot(Gammap) ).dot(K)
         C=2.*sqGamma.dot(TMP).dot(sqGammap) -unity
         d=-2.*sqGamma.dot(TMP).dot(J.T.dot(Gammap.dot(K)))
         E=4.*sqGamma.dot(TMP).dot(J.T).dot(sqGammap)
      
   
         #print "J", J, Gamma
         #print TMP.dot(Gammap), J
       #simple versions of the above matrices if J=unity and Gamma=Gammap
       #  A=np.zeros( (len(A),len(A)) )
       #  C=np.zeros( (len(A),len(A)) )
       #  E=np.ones( (len(A),len(A)) )*2.
       #  b=sqGammap.dot(K)
       #  d=-sqGamma.dot(K)
       #  print "matrices:"
       #  print A, '\n', b, '\n', C, '\n', d, "\n", E
   
         #for j in range(len(b)):
         #   print b[j], d[j], sqGamma[j][j]*K[j], (b[j]+d[j])/b[j]
   
         #initialize new OPA-object
         leng=len(b)
         L3=OPA(leng, i)                   # initialize root-node
         States=states(leng, i)           # States are all possible
         npsqrt=np.sqrt
   
         for n in States: #for each possible state, described by n(vector)
            #print "   ", n
            I_nn=0
            #need first iteration formula
            m=np.argmax(n[:leng])  #if there is no excitation: it returns 0
            # if there excists excitation in initial state: need first iteration formula
            #if n[m]!=0:
            #   print b[m]*b[m]*0.5, f[0][m]*Hartree2cm_1
            if n[m]!=0:
               # need first iteration-formula
               n_m=n[m]
               n[m]-=1 #n[m] is at least 1
               Ps=L2.getState(n)
               #if not math.isnan(Ps):
               I_nn=b[m]*Ps                                     # first term 
               print "1", b[m]*Ps 
               #print "   ", b[m], Ps
               if n[m]>0:
                  n[m]-=1
                  Ps=L1.getState(n)
                  #if not math.isnan(Ps) and abs(Ps)>1e-8:
                  I_nn+=npsqrt(2.0*(n_m-1))*A[m][m]*Ps         # second termA
                  print "2", npsqrt(2.0*(n_m-1))*A[m][m]*Ps
                  #print "   ",npsqrt(2.0*(n_m-1))*A[m][m], Ps
                  n[m]+=1
               if n[m+leng]>0:
                  n[m+leng]-=1
                  Ps=L1.getState(n)
                  n[m+leng]+=1
                  #if not math.isnan(Ps) and abs(Ps)>1e-8:
                  I_nn+=npsqrt(n[m+leng]*0.5)*(E[m][m])*Ps # second term
                  print "3", npsqrt(n[m+leng]*0.5)*(E[m][m])*Ps 
                  #print "    ",npsqrt(n[m+leng]*0.5)*(E[m][m]), Ps # second term
   
               n[m]+=1
            #else: need the other iteration-formula
            else: 
               m=np.argmax(n[leng:])  #if there is no excitation: it returns 0
               n_m=n[m+leng]
               n[m+leng]-=1
               Ps=L2.getState(n)
               I_nn=d[m]*Ps                                    # first term 
               print "1a", d[m]*Ps
               #print "   ", b[m], Ps
               if n[m+leng]>0:
                  n[m+leng]-=1
                  Ps=L1.getState(n)
                  #if not math.isnan(Ps) and abs(Ps)>1e-8:
                  I_nn+=npsqrt(2.0*(n_m-1))*C[m][m]*Ps         # second term
                  print "2a", npsqrt(2.0*(n_m-1.))*C[m][m]*Ps
                  #print "   ", npsqrt(2.0*(n_m-1))*C[m][m], Ps
                  n[m+leng]+=1
               n[m+leng]+=1
            I_nn/=npsqrt(2.*n_m)
            print "4", npsqrt(2.*n_m)
            assert not math.isnan(I_nn) ,"I_nn is not a number! I_nn: {0}\n, n:{1}\n:".format(I_nn, n)
            #if m==7 and n[m]>0:
               #print  I_nn*I_nn, sum(Gamma.dot(n[:leng])-Gamma.dot(n[leng:]))*Hartree2cm_1, m
            ##print sum(Gamma.dot(n[:leng])-Gamma.dot(n[leng:]))*Hartree2cm_1 ,I_nn*I_nn*np.exp(-sum(Gamma.dot(n[leng:]))/T) 
            #if m==7:
         #   if abs(I_nn)>1e-3:
         #      print m, n[m], n[m+leng], I_nn*I_nn , (f[1].T.dot(n[:leng])-f[0].T.dot(n[leng:]))*Hartree2cm_1
            print "   ", n
            print m, n[m], n[m+leng], I_nn, I_nn*I_nn , (f[1].T.dot(n[:leng])-f[0].T.dot(n[leng:]))*Hartree2cm_1
            L3.insert(n, I_nn)
         return L2, L3
   
      def makeLine(logging, intens, E0, T, index, ex, Gamma, Gammap, E, n):
         #F=np.zeros(( len(index),3 ))
         F=[]
         for i in xrange(len(index)):
            indi=index[i]
            if intens[i]*intens[i]*np.exp(-(Gammap[indi]*ex[i]+E0)/T) >Threshold:
               F.append([(-np.sign(E)*Gamma[indi]*(n-ex[i])
                         +np.sign(E)*Gammap[indi]*(ex[i])+E)*Hartree2cm_1,
                        intens[i]*intens[i]*np.exp(-(Gammap[indi]*ex[i]+E0)/T) ,
                        0])
               if F[-1][1]>0.0001:
                  print index[i], ex[i], n-ex[i], F[-1][1], F[-1][0]-E*Hartree2cm_1
   
         if F==[]: # no transitions with enough intensity occured.
            F=[0,0,0]
         return np.matrix(F)
      
      dimen=0
   
      linspect=[]
   
####   print Energy, sum(f[0])-sum(f[1])
      # correct for vibrational groundstates:
      #Energy+=(sum(f[0])-sum(f[1]))*.5
      L2=CalcI00(len(K), Energy)
      #this is already extracted to linspect (using side-effects)
      L1=L2 
      for i in range(1, N+1):
         L1, L2=iterate(L1, L2, Energy, i, f, J, K)
         intens, index, excitation, N=L2.extract()
         #print linspect
         #print "\n"
         #print N, N, N, N
         #print intens
         #print "\n\n\n\n"
         linspect.append(makeLine(logging, intens, E0, T, index, excitation, f[1], f[0], Energy, i))
         dimen+=len(linspect[i])+1
      spect=np.zeros((dimen,3))
      spect[:len(linspect[0])]=linspect[0]
      dimen=len(linspect[0])
      for i in xrange(1,len(linspect)):
         spect[dimen:dimen+len(linspect[i])]=linspect[i]
         dimen+=len(linspect[i])
      return spect.T
   
   def resortFCfOPA(logging, J, K, f, Energy, N, T,E0):
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
      #quantities for the iterative spectrum-calculation
      #f[0]=f[1]
   
      #J=np.eye(len(f[0]))
      #K*=np.sqrt(np.pi)/2.  # --> do it in functions_smsc.Duschinsky now.
   
      resort=np.zeros(np.shape(J))
      for i in range(len(J)):
         j=np.argmax(J[i])
         k=np.argmin(J[i])
         if J[i][j]>-J[i][k]:
            resort[i][j]=1
         else:
            resort[i][k]=-1
      J=resort.dot(J)
      K=resort.dot(K)
      for i in range(len(resort)):
         k=np.argmin(resort[i])
         if resort[i][k]==-1:
            resort[i][k]=1 #use absolute value only.
      f[1]=resort.dot(f[1])
      spect=simpleFCfOPA(logging, J, K, f, Energy, N, T,E0)
      return spect
   
   def dist_FCfOPA(logging, J, K, f, Energy, N, T, E0, threshold):
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
   
      K*=np.sqrt(np.pi)/2.  
         
      resort=np.zeros(np.shape(J))
      for i in range(len(J)):
         j=np.argmax(J[i])
         k=np.argmin(J[i])
         if J[i][j]>-J[i][k]:
            resort[i][j]=1
         else:
            resort[i][k]=-1
      J=resort.dot(J)
      K=resort.dot(K)
      for i in range(len(resort)):
         k=np.argmin(resort[i])
         if resort[i][k]==-1:
            resort[i][k]=1 #use absolute value only.
      f[1]=resort.dot(f[1])
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
         print resort
   
         J=resort.dot(Jo)
         K=resort.dot(Ka)
         f[1]=resort.dot(f1)
         f[0]=resort.dot(f0)
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
