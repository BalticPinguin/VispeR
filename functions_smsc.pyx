#!/usr/bin/python
# filename: functions_smsc.pyx
import numpy as np, re, mmap, os.path, math, sys
from scipy.linalg.lapack import dsyev as dsyev #is this faster?
#import scipy.linalg.lapack as LA
#for python-3 compatibility
from io import open 
   #debug, info, warning, error, critical
# Below are the conversion factors and fundamental constant
AMU2au=1822.88839                                          
Angs2Bohr=1/0.52917721092                                  
Hartree2GHz=6.579684e6                                     
Hartree2cm_1=219474.63 

def calcspect(logging, HR, freq, E, E0, N, M, T):
   """This is used to calculate the line spectrum assuming no mode mixing (shift only) 
   and coinciding frequencies in both electronic states.

   **PARAMETERS:**
   HR:     Huang-Rhys factors
   n:      number of modes that are considered here (with biggest HR)
   freq:   frequencies (have to be in the same order as HR
   E:      energy difference of energy surfaces
   N,M:    are the numbers of vibrational quanta can be in the modes
   All arguments are neccesary.

   **RETURNS:**
   nothing (output into /tmp/linspect)
   """
   # M,N are maximal numbers of vibrational modes (+1, since they should be arrived really; count from 0)
   N+=1
   M+=1
   def FCeqf( float Deltag, int M, int N):
      """Calculate Franck-Condon factors under assumption of equal frequencies 
      for only one vibrational mode

      PARAMETERS:
      Deltag: HR-factor of respective state
      N:      excitation number of initial state
      M:      excitation number of final state

      RETURNS:
      Franck-Condon factor of the respective transition
      """
      cdef float exg=np.exp(-Deltag/2) #actually Deltag should be >0, but is not always due to negative frequencies
      cdef float faktNM=math.factorial(M)*math.factorial(N)
      cdef float FC=0
      for x in range(int(min(N,M))+1):
         FC+=exg*math.pow(-1,N-x)*math.pow(np.abs(Deltag),(M+N)*0.5-x)/(math.factorial(M-x)*math.factorial(N-x))*\
               math.sqrt(faktNM)/math.factorial(x)
      return FC*FC
   
   def unifSpect(intens, freqs, E, FC00):
      """ Calculation of the line spectrum respecting only shift of minima (no Duschinsky rotation) 
      and assuming coinciding frequencies for initial and final state

      **PARAMETERS:**
      intens: matrix of intensities of transitions
      freqs:  matrix of respective energies

      **RETURNS**
      a 2-dimensional array with energies of transition (1. column) and their rate (2. column)
      """
      if logging[0]<1:
         logging[1].write('Spectrum\n {0}  {1}  {2}  {3}\n'.format(N, M, len(intens), len(intens[0])))
      J=len(intens[0])
      spect=np.zeros((3,len(intens)*len(intens[0])+1))
      #first transition: 0-0 transition. 
      spect[1][0]=FC00 #0->0 transition
      spect[0][0]=E
      spect[2][0]=0
      for i in range(len(intens)):#make this easier: reshapeing of intens, freqs
         for j in range(J):
            spect[2][i*J+j+1]=i+1
            spect[1][i*J+j+1]=intens[i][j]
            spect[0][i*J+j+1]=freqs[i][j]
      return spect

   n=len(HR) #=len(freq)
   assert n>0, "There is no Huang-Rhys factor larger than the respective threshold. No mode to be calculated."
   FC=np.zeros((n,M*N-1))
   uency=np.zeros((n,M*N-1)) #freqUENCY
   if approx=="TPA":
      FC2=np.zeros((1,((n+1)*(n+2)//2)*(M*N-1)*(N*M-1)))
      uency2=np.zeros((1,((n+1)*(n+2)//2)*(M*N-1)*(N*M-1)))
   #calculate 0->0 transition
   FC00=1
   for a in range(n):
      FC00*=FCeqf(HR[a],0,0)
   FC00/=N*0.97 
   #since there are N FC-like progressions, each normalised to 1
   uency00=E*Hartree2cm_1 #zero-zero transition
   if logging[0]<1:
      if approx=="OPA":
         logging[1].write("Line-spectrum in One-Particle approximation:\n")
      elif approx=="TPA":
         logging[1].write("Line-spectrum in Two-Particle approximation:\n")
      logging[1].write(u"frequency     intensity  \n")
      #logging[1].write(u" {0}\n".format(repr(E*Hartree2cm_1)+" "+repr(FC00)))
   print "0-0 trans:",FC00
   for a in range(n):
      #=first, calculate all FC-factors and finally normalise them (sum over all FC-factors is one)
      HRf0=FCeqf(HR[a], 0, 0)
      for j in range(N):
         HRf=np.zeros( int((HR[a]+1)*(M+j)) )
         s=0
         i=0
         while s<=0.96:
            HRf[i]=FCeqf(HR[a], i, j)
            s+=HRf[i]
            i+=1
         print a,j , HR[a],s, '\n', HRf
         i-=1
         for k in range(i):
            if k==0 and j==0: 
               #skip 0-0 transitions
               continue
            tmp=HRf[k]/HRf0
            FC[a][j*M+k-1]=tmp*FC00*np.exp(-(E0+freq[a]*j)/T)
            uency[a][j*M+k-1]=(E+freq[a]*(j-k))*Hartree2cm_1
            print FC[a][j*M+k-1], tmp
            if k>=M-1:
               ############# this is due to an old feature; make this more flexible in the future!!
               break
   FC00*=np.exp(-E0/T)
   spect=unifSpect(FC, uency, E*Hartree2cm_1, FC00)
   print "spectrum (first)\n", spect.T
   return spect

def CalculationHR(logging, initial, final, opt):
   """ This function gathers most essential parts for calculation of HR-factors from g09-files"""
   assert len(initial)>0, 'no initial state found!'
   assert len(final)>0, 'no final state found!'
   for i in range(len(initial)):
      assert os.path.isfile(initial[i]) and os.access(initial[i], os.R_OK),\
            initial[i]+' is not a valid file name or not readable.'
   for i in range(len(final)):
      assert os.path.isfile(final[i]) and os.access(final[i], os.R_OK),\
            final[i]+' is not a valid file name or not readable.'
   for i in range(len(initial)):
      #read coordinates, force constant, binding energies from log-files and calculate needed quantities
      dim, Coord, mass, B, A, E=ReadLog(logging, initial[i])
      if i is 0:# do only in first run
         F, CartCoord, X, P, Energy=quantity(logging, dim, len(initial)+len(final)) #creates respective quantities (empty)
         if logging[0]==0:
            logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2))
      X[i],F[i],Energy[i]=B, A, E
      CartCoord[i]=Coord
      P[i]=GetProjector(logging, X[i], dim, mass, CartCoord[i])
      if logging[0]<2:
         logging[1].write('Projector onto internal coordinate subspace\n{0}\n'.format(P[i]))
   for i in range(len(final)):
      dim, Coord, mass, B, A, E=ReadLog(logging, final[i]) 
      X[len(initial)+i], F[len(initial)+i], Energy[len(initial)+i]=B, A, E
      CartCoord[len(initial)+i]=Coord
      P[len(initial)+i]=GetProjector(logging, X[len(initial)+i], dim, mass, CartCoord[len(initial)+i])
      if logging[0]<2:
         logging[1].write('Projector onto internal coordinate subspace\n{0}\n'.format(P[len(initial)+i]))
   if logging[0]<3:
      logging[1].write('difference of minimum energy between states:'
                       ' Delta E= {0}\n'.format((Energy[0]-Energy[1])*Hartree2cm_1))
      if logging[0]<2:
         logging[1].write('Cartesion coordinates of initial state: \n{0}\n'.format( CartCoord[0].T))
         logging[1].write('Cartesion coordinates of final state: \n{0}\n Forces:\n'.format( CartCoord[1].T))
         logging[1].write('initial state: \n{0}\n'.format(F[0]))
         logging[1].write('final state: \n {0}\n'.format(F[1]))

   #Calculate Frequencies and normal modes
   f, Lsorted=GetL(logging, dim, mass,F, P)
   J, K=Duschinsky(logging, Lsorted, mass, dim, CartCoord)
   #Gauf=gaussianfreq(logging, initial, final, dim) 
   
   print K.T
   #calculate HR-spect
   HR, funi= HuangR(logging, K, f)
   if (re.search(r"makeLog", opt, re.I) is not None) is True:  
      for i in range(len(initial)): #### this needs to be enhanced
         replace(logging, initial[i], f[i], Lsorted[i])
   return HR, funi, Energy, J, K, f

def Duschinsky(logging, L, mass, dim, x):
   """
   **PARAMETERS:**
   L:    Matrix having mass-weighted normal modes as column-vectors
   mass: array of square-roots of nuclear masses (length: N)
   dim:  dimensionality of the problem: 3*N
   x:    ??

   **RETURN:**
   J:    Duschinsky-rotation matrix
   K:    displacement-vector of energy-minima in normal coordinates
   """
   J=np.zeros((len(L)-1,dim-6,dim-6))
   K=np.zeros((len(L)-1,dim-6))
   M=np.zeros((dim,dim))
   DeltaX=np.zeros((len(L)-1,dim))

   for i in range(dim):
      M[i][i]=mass[i/3] #square root of masses
   for i in range(len(J)):
      J[i]=np.dot(L[0].T, np.linalg.pinv(L[i+1].T)) ################ check: use L[i+1] instead L.-T

   for i in range(len(DeltaX)):
      DeltaX[i]=np.array(x[0]-x[i+1]).flatten('F')
      if logging[0] <1:
         logging[1].write('changes of Cartesian coordinates:(state'+repr(i)+')\n'+repr(DeltaX[i])+'\n')
      K[i]=np.dot(L[i+1].T.dot(M),DeltaX[i].T) #at the moment: mass-weighted

   np.set_printoptions(suppress=True)
   np.set_printoptions(precision=5, linewidth=138)
   if logging[0]<2:
      for i in range(len(J)):
         logging[1].write('Duschinsky rotation matrix, state '+repr(i)+\
               '  :\n'+ repr(J[i])+'  :\n'+ repr(J[i][:4].T[11:25])+\
               '\nDuschinsky displacement vector:\n'+ repr(K[i])+'\n')
   return J, K 

def gaussianfreq(logging, initial, final , dim):
   """Extraction of the frequencies from g09-log file"""
   ##read frequencies calculated by g09 from file
   #Gauf=of.gaussianfreq(ContntInfo, dim) 
   #Gauf/=Hartree2cm_1  #convert to atomic units
   f=np.zeros((len(initial)+len(final), dim-6))
   for i in range(len(initial)): 
      files=open(initial[i], "r") #open file and map it for better working
      mapedlog=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ) # i-th file containing freq calculations
      files.close
      f1=re.findall(r" Frequencies -- [\d .-]+", mapedlog, re.M)# 
      f2=[re.findall(r"[- ]\d+.\d+", f1[j]) for j in range(len(f1))]
      s=0
      for j in range(len(f2)):
         f[i][s:s+len(f2[j])]=f2[j]
         s+=len(f2[j])
   for i in range(len(final)): 
      files=open(final[i], "r") #open file and map it for better working
      mapedlog=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ) # i-th file containing freq calculations
      files.close
      f1=re.findall(r" Frequencies -- [\d .-]+", mapedlog, re.M)# 
      f2=[re.findall(r"[- ]\d+.\d+", f1[j]) for j in range(len(f1))]
      s=0
      for j in range(len(f2)):
         f[i][s:s+len(f2[j])]=f2[j]
         s+=len(f2[j])
   return f

def GetL(logging, dim, mass, F, D):
   """ Function that calculates the frequencies and normal modes from force constant matrix 
   with and without projection onto internal degrees of freedom

   **argumets**
   1. The dimensions of force-constant matrix
   2. square-root of masses dim/3-dimensional array
   3. force-constant matrix
   4. projection-matrix onto vibrations

   **return**
   1. matrix of all normal modes (including global modes) 
   2. matrix of vibr. normal modes (including global modes) 
   3. frequencies 
   4. massweighted L for comparison with the respective matrix from the g09 log-file
   """
   # Defining arrays
   L=np.zeros(( len(F), len(F[0]), len(F[0])-6 )) 
   Ltest=np.zeros(( len(F), len(F[0]), len(F[0])-6 )) 
   Lsorted=np.zeros(( len(F), len(F[0]), len(F[0])-6 ))
   f=np.zeros(( len(F), len(F[0])-6 ))
   Ltemp=np.zeros(( len(F[0]), len(F[0])-6 ))
   ftemp=np.zeros(len(F[0]-6))

   for i in range(len(F)):
      #### the condition number of F[i] is some millions...
      #ftemp,Ltemp=np.linalg.eig(F[i])
      #ftemp,Ltemp=np.linalg.eigh(F[i])
      ftemp,Ltemp,info=dsyev(F[i]) #this seems to be the best function

      #assert np.any(ftemp< 0) or np.any(np.imag(ftemp)!=0),\
#       'Frequencies smaller than 0 occured. Please check the input-file!! {0}'.format(ftemp)
      if np.any(ftemp<0):
         if logging[0]<4:
            logging[1].write('Frequencies smaller than 0 occured. The absolute'
                        ' values are used in the following.\n{0}\n'.format(ftemp))
         ftemp=np.abs(ftemp)

      M=np.zeros((dim,dim))
      for j in range(0,dim):
         M[j,j]=1/mass[j/3]
      #Lcart=np.dot(M,np.dot(D[i],np.real(Ltemp)))
      Lcart=M.dot(np.real(Ltemp)) # there is no need for this projector!!!???
      #Lcart=np.real(Ltemp)
      for j in range(0,dim):
         norm=np.sum(Lcart.T[j]*Lcart.T[j])
         if np.abs(norm)>1e-12:
            Lcart.T[j]/=np.sqrt(norm)

      index=np.argsort(np.real(ftemp),kind='heapsort') # ascending sorting f
      if logging[0]<1:
         logging[1].write("Frequencies (cm-1) \n"+ repr(np.sqrt(np.abs(ftemp[index]))*Hartree2cm_1))
      f[i]=np.real(ftemp[index]).T[:].T[6:].T
      Lsorted[i]=(Lcart.T[index].T)[:].T[6:].T
      L[i]=(Ltemp.T[index].T)[:].T[6:].T
      if logging[0]<1:
         logging[1].write("Normalized, sorted and truncated Lcart\n"+ repr(L[i]))

      for j in range(len(f[i])):
         f[i][j]=np.sign(f[i][j])*np.sqrt(np.abs(f[i][j]))
      if logging[0]<2:
         logging[1].write("After projecting onto internal coords subspace\n"+"Frequencies (cm-1)\n"+\
               repr(f[i]*Hartree2cm_1)+"\nL-matrix \n"+ repr(L[i]))
   return f, Lsorted

def GetProjector(logging, X, dim, m, Coord):
   D=np.zeros((dim,6))
   for k in range(3):# first three rows in D: The translational vectors
      for j in range(dim/3):
         D[3*j+k][k]=m[j]
   for k in range(dim):# next three rows in D: The rotational vectors
      D[k][3:6]=(np.cross(np.dot(X,Coord)[:].T[k/3],X[:].T[k%3]))*m[k/3]
   if logging[0]<1:
      logging[1].write("Original translational and rotational displacement vectors:\n{0}\n".format(D[3:13].T))
   AOE=gs(np.array(D)) #orhogonalize it
   ones=np.identity(dim)
   one_P=ones-np.dot(AOE,AOE.T)
   prob_vec=(AOE.T[1]+AOE.T[4]+AOE.T[0]+AOE.T[5]).T #what is this actually??
   assert not np.any(np.abs(prob_vec-np.dot(np.dot(AOE,AOE.T),prob_vec))>0.00001), \
            'Translations and rotations are affected by projection operator.'+\
            repr(np.abs(prob_vec-np.dot(np.dot(AOE,AOE.T),prob_vec)))
   assert not  np.any(np.abs(np.dot(one_P,prob_vec))>0.00001), \
            "Projecting out translations and rotations from probe vector"
   return one_P

def gradientHR(logging, initial, final, opt):
   """ This function gathers most essential parts for calculation of HR-factors from g09-files"""
   assert len(initial)>0, 'no initial state found!'
   assert len(final)>0, 'no final state found!'
   initial=initial[0]
   final=final[0]
   assert os.path.isfile(initial) and os.access(initial, os.R_OK),\
            initial+' is not a valid file name or not readable.'
   assert os.path.isfile(final) and os.access(final, os.R_OK),\
            final+' is not a valid file name or not readable.'
   dim, Coord, mass, B, A, E=ReadLog(logging, initial)
   F, CartCoord, X, P, Energy=quantity(logging, dim, 2 ) #creates respective quantities (empty)
   if logging[0]==0:
      logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2))
   X[0],F[0],Energy[0]=B, A, E
   F[1]=F[0] #force constant matrix in both states coincides
   Grad, E=ReadLog2(logging, final) 
   Energy[1]=E
   #read coordinates, force constant, binding energies from log-files and calculate needed quantities
   if logging[0]<3:
      logging[1].write('difference of minimum energy between states:'
                       ' Delta E= {0}\n'.format((Energy[0]-Energy[1])*Hartree2cm_1))
      if logging[0]<2:
         logging[1].write('initial state: \n{0}\n'.format(F[0]))

   #Calculate Frequencies and normal modes
   f, Lsorted=GetL(logging, dim, mass,F, P)
   K=GradientShift(logging, Lsorted, mass, Grad)
      ##################################################
   
   #calculate HR-spect
   HR, funi= HuangR(logging, K, f)
   if (re.search(r"makeLog", opt, re.I) is not None) is True:  
      for i in range(len(initial)): #### this needs to be enhanced
         replace(logging, initial[i], f[i], Lsorted[i])
   return HR, funi, Energy, K, f

def GradientShift(logging, L, mass, Grad):
   K=L[1].T.dot(Grad).T
   for i in range(len(K)):
     K[i]*=mass[i//3]
   print "K", K.T
   return K

def gs(A):
   """This function does row-wise Gram-Schmidt orthonormalization of matrices. 
   code for Gram-Schmidt adapted from iizukak, see https://gist.github.com/iizukak/1287876
   """
   def proj(v1, v2):
      return map(lambda x : x *(np.dot(v2,
               v1) / np.dot(v1, v1)) , v1)

   X=A.T # I want to orthogonalize row-wise
   Y = []
   for i in range(len(X)):
      temp_vec = X[i]
      for inY in Y :
         proj_vec = proj(inY, X[i])
         temp_vec = map(lambda x, y : x - y, temp_vec, proj_vec)
      Y.append( temp_vec/np.linalg.norm(temp_vec)) # normalise vectors
   return np.matrix(Y).T # undo transposition in the beginning

def HuangR(logging, K, f): #what is with different frequencies???
   """ Function that calculates the Huang-Rhys factors for all vibrational states

   **Arguments**
   1. The displacements of minima in internal coordinates
   2. The frequencies of respective modes

   **returns**
   return sortuni, funi, sortmulti, sortfG, sortfE
   1. HR-factors for coinciding frequencies sorted by size (decreasing)
   2. respective frequencies for 1. (same order)
   3. HR-factors for different frequencies of involved electronic states sorted by size (decreasing)
   4. respecivp frequencies of initial state for 3 (same order)
   5. respecivp frequencies of final state for 3 (same order)
   """
   sortuni=np.zeros((len(K),len(K[0])))
   funi=np.zeros((len(K),len(K[0])))
   uniHRall=[]
   uniFall=[]
   for i in range(len(K)):
      unif=K[i]*K[i]*f[i+1]*0.5
      index=np.argsort(unif, kind='heapsort')
      sortuni[i]=unif[index]
      funi[i]=f[i+1][index]
      if np.any(funi)<0:
         if logging[0]<4:
            logging[1].write('ATTENTION: some HR-factors are <0.\
                  In the following their absolute value is used.')
         funi[i]=np.abs(funi[i])
      uniHR=[]
      uniF=[]

      logging[1].write(u'HR-fact           freq\n')
      for j in range(len(sortuni[i])):
         #select all 'big' HR-factors 
         if sortuni[i][-j]>0.001:
            uniHR.append(sortuni[i][-j])
            uniF.append(funi[i][-j])
            logging[1].write(u"{0}   {1}\n".format(sortuni[i][-j], funi[i][-j]*Hartree2cm_1))
      uniHRall.append(uniHR)
      uniFall.append(uniF)
   return uniHRall, uniFall

def quantity(logging, dim, num_of_files):
   F=np.zeros((num_of_files, dim, dim)) 
   CartCoord=np.zeros((num_of_files, 3, dim/3))
   X=np.zeros((num_of_files, 3,3))
   P=np.zeros((num_of_files, dim,dim))
   Energy=np.zeros(num_of_files)
   return F, CartCoord, X, P, Energy

def ReadHR(logging, HRfile):
   """ This function reads the HR-factors and electronic transition energy from a given file and brings them into a 
   similar structure as they are used in the 'smallscript'.

   **PARAMETERS**
   logging:     array containing the mode (how detailed the printed information is) (first element) and the file-object
   HRfile:      the file where the information is found

   **RETURNS**
   initial:     a dummy-array that originally contains information about the inital states. Here at the moment only one
            is allowed and only its length is relevant in the further programme.
   HRm: a 2-dimensional array containing all Huang-Rhys-factors
   freqm:       analogously to HRm, containing the respective frequencies
   Energy:      the energy-difference of the electronic states. (in atomic units)

   """
   assert os.path.isfile(HRfile) and os.access(HRfile, os.R_OK),\
            HRfile+' is not a valid file name or not readable.'
   fi=open(HRfile, "r")
   f=mmap.mmap(fi.fileno(), 0, prot=mmap.PROT_READ)
   fi.close()
   Energy=re.findall(r"(?<=Delta E=)[ \d\.\-]*", f, re.I)
   assert len(Energy)==1, "Please specify only one energy-difference!"
   Energy=float(Energy[0])
   Energy/=Hartree2cm_1
   HR=[]
   funi=[]
   HRfreq=re.findall(r"HR-fact[\s]*freq[\s]*\n[\n\d\.\s]*", f, re.I)
   assert len(HRfreq)==1, "The file-format could not be read. exit now"
   HRf=re.findall(r"(?<=\n)[\d\.]*[\s]+[\d\.]*", HRfreq[0], re.I)
   for i in range(len(HRf)):
      line=re.findall(r"[\d.]+",HRf[i], re.I)
      HR.append(float(line[0]))
      funi.append(float(line[1])/Hartree2cm_1)
   initial=['excited']
   #the following is just to be consistent with structure of HR calculated in first part
   HRm=np.zeros((1,len(HR)))
   HRm[0]=HR
   freqm=np.zeros((1,len(HR)))
   freqm[0]=funi
   return initial, HRm, freqm, Energy

def ReadLog(logging, fileN):
   # Mapping the log file
   files=open(fileN, "r")
   log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   files.close

   # Determine atomic masses in a.u. Note mass contains sqrt of mass!!!
   temp=[]
   temp=re.findall(r' Number     Number       Type             X           Y           Z[\n -.\d]+', log)
   tmp=re.findall(r'[ -][\d]+.[\d]+', temp[-1])

   atmwgt=re.findall(r"AtmWgt= [\d .]+",log)
   mtemp=[]
   foonum=0
   for j in range(len(atmwgt)/2): # because atomic masses are printed twize in log-files...
      mtemp.append(re.findall(r'[\d.]+',atmwgt[j]))
   dim=0
   for j in range(len(mtemp)):
         dim+=len(mtemp[j]) # dim will be sum over all elements of temp
   dim*=3
   if dim!=len(tmp):
      # this is necessary since they are not always printed twice...
      atmwgt=re.findall(r"AtmWgt= [\d .]+",log)
      mtemp=[]
      foonum=0
      for j in range(len(atmwgt)): 
         mtemp.append(re.findall(r'[\d.]+',atmwgt[j]))
      dim=0
      for j in range(len(mtemp)):
            # dim will be sum over all elements of temp
            dim+=len(mtemp[j]) 
      dim*=3
   #if still something is wrong with the dimensionality:
   assert len(tmp)==dim, 'Not all atoms were found! Something went wrong...{0}, {1}'.format(len(tmp),dim)

   mass=np.zeros(dim/3) # this is an integer since dim=3*N with N=atomicity
   for j in range(len(mtemp)):
      for k in range(len(mtemp[j])):
         mass[k+foonum]=np.sqrt(float(mtemp[j][k])*AMU2au) #elements in m are sqrt(m_i) where m_i is the i-th atoms mass
      foonum+=len(mtemp[j])
   assert not np.any(mass==0) , "some atomic masses are zero. Please check the input-file! {0}".format(mass)
   if logging[0]<2:
      logging[1].write("Number of atoms: {0}\nNumber of vibrational "
                        "modes: {1} Sqrt of masses in a.u. as read from log file\n{2}\n".format(dim/3,dim,mass))
   # Reading Cartesian coordinates
   Coord=np.zeros((3, dim/3))
   MassCenter=np.zeros(3)
   for j in range(len(tmp)):
      Coord[j%3][j/3]=tmp[j]
   for j in range(3):
      MassCenter[j]=np.sum(Coord[j]*mass)
      MassCenter[j]/=np.sum(mass) #now it is cartesian center of mass
   if logging[0]<2:
      logging[1].write("Cartesian (Angstrom) coordinates before alignment to center of "
                       "mass\n {0} \nCenter of mass coordinates (Angstrom):\n{1}\n".format(Coord.T, MassCenter))
   
   for j in range(3):#displacement of molecule into center of mass:
      Coord[j]-=MassCenter[j] # if commented we get rotational constants in agreement with Gaussian log
   Coord*=Angs2Bohr
   if logging[0]<2:
      logging[1].write("Cartesian coordinates (a.u.) in center of mass system\n{0}\n".format(Coord.T))

   # Getting tensor of inertia, transforming to principlas axes
   moi=np.zeros((3,3))# this is Moment Of Inertia
   for j in range(3):
      for k in range(3):
         if k is j:
            moi[j][k]=np.sum(mass*mass*(Coord[0]*Coord[0]+\
                     Coord[1]*Coord[1]+Coord[2]*Coord[2]-\
                     Coord[j]*Coord[k]))
         else:
            moi[j][k]=np.sum(mass*mass*(Coord[j]*Coord[k]))
   if logging[0]<1:
      logging[1].write("Moments of intertia as read from log file\n,{0}".format(moi))
   diagI,X=np.linalg.eig(moi) # this can be shortened of course!
   index=np.argsort(diagI, kind="heapsort")
   #X=np.matrix(X[index]) #sorting by eigenvalues
   X=np.matrix(X) #sorting by eigenvalues
   diagI=diagI[index]
   if logging[0]<2:
      logging[1].write("Moments of inertia (a.u.) in principle axes\n {0}\nRotational "
                       "constants (GHz) in principle axes\n {1} Rotation matrix\n{2}"
                       .format(diagI.T,1/(2*diagI.T)*Hartree2GHz, X))
   # Reading of Cartesian force constant matrix  
   f=re.findall(r"Force constants in Cartesian coordinates: [\n\d .+-D]+", log, re.M)
   if f==[]:
      #if Freq was not specified in Gaussian-file:
      f=re.findall(r"The second derivative matrix:[XYZ\n\d .-]+", log, re.M)
      #try to find matrix from option "Force"
      assert f!=[], 'The input-file does not contain information on the force-constants!'
      #else: error message. The matrix is really needed...
      f_str=str([f[-1]])#[2:-2]
      lines=f_str.strip().split("\\n")
      F=np.zeros((dim,dim))
      n=0
      k=0
      line=0
      for i in range(2,len(lines)):
         if i == dim+k-5*n+2: 
            #these are those lines where no forces are written to
            k=i-1
            n+=1
            line=n
            continue
         elements=lines[i].split()[1:] #don't use the first element
         #line=int(re.findall(r"[\d]+", lines[i].split()[0])[0])+(i-2-n)
         line+=1
         for j in range(len(elements)):
            F[line-1][j+5*n]=float(elements[j])
            F[j+5*n][line-1]=float(elements[j])
   else:
      f_str=str([f[-1]])#[2:-2]
      lines=f_str.strip().split("\\n")
      F=np.zeros((dim,dim))
      n=0
      k=0
      for i in range(2,len(lines)):
         if i == dim+k-5*n+2: 
            k=i-1
            n+=1
            continue
         elements=lines[i].replace('D','e').split()
         for j in range(1,len(elements)):
            F[int(elements[0])-1][j-1+5*n]=float(elements[j])
            F[j-1+5*n][int(elements[0])-1]=float(elements[j])
   if logging[0]<1:
      logging[1].write('F matrix as read from log file\n{0} \n'.format(F))
   for i in range(0,dim):
      for j in range(0,dim):
         F[i][j]/= (mass[i/3]*mass[j/3]) 

   Etemp=re.findall(r'HF=-[\d.\n ]+', log, re.M)
   assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
   if re.search(r'\n ', Etemp[-1]) is not None:
      Etemp[-1]=Etemp[-1].replace("\n ", "") 
   if logging[0]<=1:
      logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
   E=-float(re.findall(r'[\d.]+', Etemp[-1])[0]) #energy is negative (bound state)
   return dim, Coord, mass, X, F, E

def ReadLog2(logging, final):
   files=open(final, "r")
   log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   print "in readlog2"
   files.close
   Etemp=re.findall(r'HF=-[\d.\n ]+', log, re.M)
   assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
   if re.search(r'\n ', Etemp[-1]) is not None:
      Etemp[-1]=Etemp[-1].replace("\n ", "") 
   if logging[0]<=1:
      logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
   E=-float(re.findall(r'[\d.]+', Etemp[-1])[0]) #energy is negative (bound state)
  
   #grad=re.findall(r"Final forces over variables, Energy=[\-\+ :\.\d D\n]+", log)
   #Grad=re.findall(r"(?<=\:\n)[\-\+\.\d D\n]+", grad[0])
   #Grad=re.findall(r"[\-\d\.]+D[-+][\d]{2}", Grad[0])
   #grad=np.zeros((len(Grad),1))
   #for i in range(len(Grad)):
      #element=Grad[i].replace('D','e')
      #grad[i]=float(element)
   grad=re.findall(r"Number     Number              X              Y              Z\n [\-\.\d \n]+",log)
   Grad=re.findall(r"[\-\d\.]+[\d]{9}", grad[0])
   grad=np.zeros((len(Grad),1))
   for i in range(len(Grad)):
      grad[i]=float(Grad[i])
   return grad, E

def replace(logging, files, freq, L):
   """ This function creates a new file (determined by files, ending with 
   ".rep" and copies the log-file (files) into it, replacing the frequencies and 
   normal modes by those calculated by smallscript.
   The function is suited to test, whether these results coincide qualitatively with the gaussian's.

   ** PARAMETERS: **
   logging: specifies the level of logging-outputs
   log:   i
   files: file taken as basis ('.rep' added to be used for replacements)
   freq:  frequencies to be inserted
   L:     normal modes to be inserted  

   no return-statements (results are written to file)

   **NOTE:**
   The code is originally from stevaha (http://stackoverflow.com/questions/1597649/replace-strings-in-files-by-python)
   """
   freq*=Hartree2cm_1
   with open(files) as f:
      out_fname = files + ".rep"
      out = open(out_fname, "w")
      s=0
      t=0
      u=-3
      for line in f:
         if re.search(r'Frequencies -- [\d .-]+', line) is not None:
            t=0 #reset t, when frequencies occur
            u+=3 #
          #  if logging[0]<1:
          #     logging[1].write('frequencies not yet written to file:'+ repr(len(freq[s:].T))+ repr(freq[s:].T))
            if len(freq[s:].T)> 2: # there are at least three more frequencies
               out.write(re.sub(r'Frequencies -- [\d .-]+',
                     'Frequencies --'+'    '+str("%.4f" % freq[s])+'              '\
                     +str("%.4f" % freq[s+1])+'               '+str("%.4f" % freq[s+2]), line))
            elif len(freq[s:].T)== 2: # there are only two frequencies left
               out.write(re.sub(r'Frequencies -- [\d .-]+',
                     'Frequencies --'+'    '+str("%.4f" % freq[s])+'               '\
                     +str("%.4f" % freq[s+1]), line))
            elif len(freq[s:].T)== 1: # there is just one additional freq
               out.write(re.sub(r'Frequencies -- [\d .-]+',
                     'Frequencies --'+'   '+str("%.4f" % freq[s]), line))
            s+=3
         elif re.search(r'[ ]+\d+[ ]+\d+[ -]+\d.\d\d[ -]+\d.\d\d+[ \d.-]+', line) is not None:
            if len(L[t][u:].T)> 2: # there are at least three more frequencies
               out.write(re.sub(r'[\d .-]+', '     '+repr(t/3+1)+'    '+repr(s)+'   '+
                  '  '+str("%.6f" % L[t+0][u+0])+'  '+str("%.6f" % L[t+1][u+0])+' '+
                       str("%.6f" % L[t+2][u+0])+
                  '  '+str("%.6f" % L[t+0][u+1])+'  '+str("%.6f" % L[t+1][u+1])+' '+
                       str("%.6f" % L[t+2][u+1])+
                  '  '+str("%.6f" % L[t+0][u+2])+'  '+str("%.6f" % L[t+1][u+2])+' '+
                       str("%.6f" % L[t+2][u+2]), line))
            elif len(L[t][u:].T)== 2:
               out.write(re.sub(r'[\d .-]+', '     '+repr(t/3+1)+'    '+repr(s)+'   '+
                  '  '+str("%.6f" % L[t+0][u+0])+'  '+str("%.6f" % L[t+1][u+0])+' '+
                       str("%.6f" % L[t+2][u+0])+
                  '  '+str("%.6f" % L[t+0][u+1])+'  '+str("%.6f" % L[t+1][u+1])+' '+
                       str("%.6f" % L[t+2][u+1]), line))
            elif len(L[t][u:].T)== 1:
               out.write(re.sub(r'[\d .-]+', '     '+repr(t/3+1)+'    '+repr(s)+'   '+
                  '  '+str("%.6f" % L[t+0][u+0])+'  '+str("%.6f" % L[t+1][u+0])+' '+
                       str("%.6f" % L[t+2][u+0]), line))
            t+=3 # 
         else: 
            out.write(re.sub('replace nothing','by nothing', line)) #just write line as it is
      out.close()

version=1.2
# End of functions_smsc.pyx
