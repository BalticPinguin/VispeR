#!/usr/bin/python
# filename: functions_smsc.py
import numpy as np
import re, mmap, os.path, math, sys
import readLogs as rl
#for python-3 compatibility
from io import open 
# Below are the conversion factors and fundamental constant
AMU2au=1822.88839                                          
Angs2Bohr=1/0.52917721092                                  
Hartree2GHz=6.579684e6                                     
Hartree2cm_1=219474.63 

def calcspect(logging, HR, freq, E, E0, N, M, T):
   """This is used to calculate the line spectrum assuming no mode mixing (shift only) 
   and coinciding frequencies in both electronic states.

   **PARAMETERS:**
   logging:This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
           and 4- only main information) and logging[1] is the file, already opened, to write the information in.
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
   def FCeqf( Deltag, M, N):
      """Calculate Franck-Condon factors under assumption of equal frequencies 
      for only one vibrational mode

      PARAMETERS:
      Deltag: HR-factor of respective state
      N:      excitation number of initial state
      M:      excitation number of final state

      RETURNS:
      Franck-Condon factor of the respective transition
      """
      exg=np.exp(-Deltag/2) 
      fact=math.factorial
      faktNM=fact(M)*fact(N)
      FC=0
      for x in range(int(min(N,M))+1):
         FC+=exg*math.pow(-1,N-x)*math.pow(Deltag,(M+N)*0.5-x)/(fact(M-x)*fact(N-x))*\
               math.sqrt(faktNM)/fact(x)
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
      J=len(intens[0])
      spect=np.zeros((3,len(intens)*J+1))
      #first transition: 0-0 transition. 
      spect[1][0]=FC00 #0->0 transition
      spect[0][0]=E
      spect[2][0]=0
      for i in range(len(intens)):#make this easier: reshapeing of intens, freqs
         for j in xrange(J):
            spect[2][i*J+j+1]=i+1
            spect[1][i*J+j+1]=intens[i][j]
            spect[0][i*J+j+1]=freqs[i][j]
      return spect

   n=len(HR) #=len(freq)
   assert n>0, "There is no Huang-Rhys factor larger than the respective threshold. No mode to be calculated."
   FC=np.zeros((n,M*N-1))
   uency=np.zeros((n,M*N-1)) #freqUENCY
   #calculate 0->0 transition
   FC00=10
   uency00=E*Hartree2cm_1 #zero-zero transition
   loggingwrite=logging[1].write #avoid dots!
   if logging[0]<1:
      loggingwrite("Line-spectrum in One-Particle approximation:\n")
      loggingwrite(u" %f   %f  %f\n"%(E*Hartree2cm_1, FC00, 0))
   npexp=np.exp #avoiding dots accelerates python quite a lot
   #here a goes over all modes
   for a in xrange(n):
      temp=FCeqf(HR[a],0,0)
      for j in range(N):  
         s=0 # this is a temp variable for effinciency purpose
         for i in range(M):
            if i==0 and j==0:
               ##skip 0-0 transitions
               continue
            tmp=FCeqf(HR[a], i, j)/temp
            s+=tmp*temp # s+=FCeqf(HR[a], i, j) but more efficient than calling it twice
            if s>=0.96:
               # than all transitions contributing in intensity are done
               break
            FC[a][j*M+i-1]=tmp*FC00*npexp(-(E0+freq[a]*i)/T)
            uency[a][j*M+i-1]=(E+freq[a]*(i-j))*Hartree2cm_1
            #if logging[0]<1:
               #logging[1].write(u" {0}\n".format(repr(uency[a][j*M+i-1])+" "+repr(FC[a][j*M+i-1])+" "+repr(a)))
   FC00*=npexp(-E0/T)
   spect=unifSpect(FC, uency, E*Hartree2cm_1, FC00)
   return spect

def CalculationHR(logging, initial, final, opt, HRthresh):
   """ This function gathers most essential parts for calculation of HR-factors from g09-files.
   That is: read neccecary information from the  g09-files and calculate HR-factors as well as the 
   Duschinsky-rotation matrix and the shift between minima (needed later if the option Duschinsky is specified)

   **PARAMETERS**
   logging: This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
            and 4- only main information) and logging[1] is the file, already opened, to write the information in.
   initial: name of the file containing the initial state's geometry
   final:   name of the file containing the initial state's geometry
   opt:     options that are given to this calculation; especially it is of interest, whether there should be frequencies and/or normal
            modes to be read from the g09-files.
   HRthresh threshold for HR-factors (only those larger than this are taken into account)

   **RETURNS**
   HR:      Huang-Rhys factors, sorted by size
   funi:    vibrational frequencies sorted same as HR 
   Energy:  Energy difference between the two states
   J:       Duschinsky-rotation matrix
   K:       shift between the states (in normal coordinates)
   f:       frequencies of the vibrational modes, sorted by size of the frequencies

   """                                                                                             
   assert len(initial)==1, 'there must be one initial state'
   assert len(final)==1, 'there must be one final state'
   assert os.path.isfile(initial[0]) and os.access(initial[0], os.R_OK),\
            initial[0]+' is not a valid file name or not readable.'
   assert os.path.isfile(final[0]) and os.access(final[0], os.R_OK),\
            final[0]+' is not a valid file name or not readable.'
   # if the input-files are G09-formcheck-files (recognised by ending):
   if ( ".fchk" in final[0]) and (".fchk" in initial[0]):
      dim, Coord, mass, A, E=rl.ReadGO9_fchk(logging, initial[0])
      F, CartCoord, P, Energy=quantity(logging, dim, 2) #creates respective quantities (empty)
      if logging[0]==0:
         logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2)+"\n")
      F[0],Energy[0]=A, E
      CartCoord[0]=Coord
      dim, Coord, mass, A, E=rl.ReadGO9_fchk(logging, final[0]) 
      F[1], Energy[1]=A, E
      CartCoord[1]=Coord
   #else, test what kind of file was given: G09, GAMESS or NWChem
   else:
      with open(initial[0], "r+b") as f: #open file as mmap
         mapping = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
         for line in iter(mapping.readline, ""): #go through every line and test characteristic part
            if "GAMESS" in line: #if it is found: read important quantities from file
               dim, Coord, mass, A, E=rl.ReadGAMESS(logging, initial[0])
               F, CartCoord, P, Energy=quantity(logging, dim, 2) #creates respective quantities (empty)
               if logging[0]==0:
                  logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2)+"\n")
               F[0],Energy[0]=A, E
               CartCoord[0]=Coord
               dim, Coord, mass, A, E=rl.ReadGAMESS(logging, final[0]) 
               F[1], Energy[1]=A, E
               CartCoord[1]=Coord
               break
            elif "Northwest Computational Chemistry Package (NWChem)" in line:
               dim, Coord, mass, A, E=rl.ReadNWChem(logging, initial[0])
               F, CartCoord, P, Energy=quantity(logging, dim, 2) #creates respective quantities (empty)
               if logging[0]==0:
                  logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2)+"\n")
               F[0],Energy[0]=A, E
               CartCoord[0]=Coord
               dim, Coord, mass, A, E=rl.ReadNWChem(logging, final[0]) 
               F[1], Energy[1]=A, E
               CartCoord[1]=Coord
               break
            elif "Gaussian(R)" in line:
               dim, Coord, mass, A, E=rl.ReadG09(logging, initial[0])
               F, CartCoord, P, Energy=quantity(logging, dim, 2) #creates respective quantities (empty)
               if logging[0]==0:
                  logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2)+"\n")
               F[0],Energy[0]=A, E
               CartCoord[0]=Coord
               dim, Coord, mass, A, E=rl.ReadG09(logging, final[0]) 
               F[1], Energy[1]=A, E
               CartCoord[1]=Coord
               break
         else: # is there some error??  probably this message is printed also if files were recognised
            print "file type not recognised"
   #read coordinates, force constant, binding energies from log-files and calculate needed quantities
   if logging[0]<3:
      logging[1].write('difference of minimum energy between states:'
                       ' Delta E= {0}\n'.format((Energy[0]-Energy[1])*Hartree2cm_1))
      if logging[0]<2:
         logging[1].write('Cartesion coordinates of initial state: \n{0}\n'.format( CartCoord[0].T/Angs2Bohr))
         logging[1].write('Cartesion coordinates of final state: \n{0}\n Forces:\n'.format( CartCoord[1].T/Angs2Bohr))
         logging[1].write('initial state: \n{0}\n'.format(F[0]))
         logging[1].write('final state: \n {0}\n'.format(F[1]))

   #Calculate Frequencies and normal modes
   f, Lsorted, Lmassw=GetL(logging, dim, mass, F)
   J, K=Duschinsky(logging, Lsorted, mass, dim, CartCoord)
   extra=re.findall(r"g09Vectors",opt, re.I)
   if extra!=[]:
      g09L=rl.getGaussianL(final, mass, dim)
      g09f=rl.getGaussianf(final,dim)
      #rl.replace(logging, initial[0], g09f, g09L)
   elif re.search(r"g09Vector",opt, re.I) is not None:
      g09L=rl.getGaussianL(final, mass, dim)
      g09f=rl.getGaussianf(final,dim)
      #rl.replace(logging, initial[0], g09f, g09L)

   #comparet  f, g09f  and Lsorted/Lcart with g09L --> which is which??
   
   #calculate HR-spect
   HR, funi= HuangR(logging, K, f, HRthresh)
   if (re.search(r"makeLog", opt, re.I) is not None) is True:  
      rl.replace(logging, initial[0], f[0], Lmassw[0])
   return HR, funi, Energy, J, K, f

def Duschinsky(logging, L, mass, dim, x):
   """
   This function calculates the shift between two electronic states (whose geometry is known, see x) as well as the
   Duschinsky-rotation matrix.
   **PARAMETERS:**
   logging: This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
            and 4- only main information) and logging[1] is the file, already opened, to write the information in.
   L:       Matrix having mass-weighted normal modes as column-vectors
   mass:    array of square-roots of nuclear masses (length: N)
   dim:     dimensionality of the problem: 3*N
   x:       cartesian coordinates of the states of interest

   **RETURN:**
   J:    Duschinsky-rotation matrix
   K:    displacement-vector of energy-minima in normal coordinates
   """
   J=np.zeros((len(L)-1,dim-6,dim-6))
   K=np.zeros((len(L)-1,dim-6))
   M=np.zeros((dim,dim))
   DeltaX=np.zeros((len(L)-1,dim))

   for i in range(dim):
      M[i][i]=mass[i//3] #square root of masses
   for i in range(len(J)):
      #Jtemp=np.dot(L[0].T, np.linalg.pinv(L[i+1].T)) ################ check: use L[i+1] instead L.-T
      J[i]=np.dot(L[0].T, L[i+1]) 

   DeltaX[0]=np.array(x[0]-x[1]).flatten('F')
   if logging[0] <1:
      logging[1].write('changes of Cartesian coordinates:\n'\
            +repr(DeltaX[0])+'\n')
   K[0]=(DeltaX[0].dot(M)).dot(L[0])*0.8676157 ##what factor is this??
   #K[0]=(DeltaX[0].dot(M)).dot(L[0])/1.63 ##what factor is this??
   #K[0]=L[0].T.dot(M).dot(DeltaX[0])
   
   #np.set_printoptions(suppress=True)
   #np.set_printoptions(precision=5, linewidth=138)
   if logging[0]<2:
      for i in range(len(J)):
         logging[1].write('Duschinsky rotation matrix, state '+repr(i)+':\n')
         k=range(0,dim)
         s=0
         t=min(s+5,dim)
         while s<dim:
            for temp in range(s,t):
               logging[1].write("               %d "%(k[temp]+1))
            logging[1].write("\n")
            for j in range(len(J[i])):
               logging[1].write(" %03d"%(j+1))
               for temp in range(s,t):
                  logging[1].write("    %+.6e"%(J[i][j][k[temp]]))
               logging[1].write("\n")
            s+=t
            t=min(s+5,dim-6)
         logging[1].write('\nDuschinsky displacement vector:\n')

         for j in range(len(K[i])):
            logging[1].write("  %d    %e\n"%(j+1, K[i][j]))
         #logging[1].write('\nDuschinsky displacement vector:\n'+ repr(K[i])+'\n')

   return J, K 

def GetL(logging, dim, mass, F):
   """ Function that calculates the frequencies and normal modes from force constant matrix 
   with and without projection onto internal degrees of freedom

   **argumets**
   logging: This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
            and 4- only main information) and logging[1] is the file, already opened, to write the information in.
   dim      The dimensions of force-constant matrix
   mass     square-root of masses dim/3-dimensional array
   F        force-constant matrix

   **return**
   return f, Lsorted, Lcart
   f        two vertors (f[0] and f[1]) with the frequencies
   Lsorted  matrix of vibr. normal modes (transformation matrix)
   Lcart    massweighted L for comparison with the respective matrix from the g09 log-file
   """
   # Defining arrays
   L=np.zeros(( len(F), len(F[0]), len(F[0])-6 )) 
   Lmass=np.zeros(( len(F), len(F[0]), len(F[0])-6 ))
   f=np.zeros(( len(F), len(F[0])-6 ))

   for i in range(len(F)):
      # here one can choose between the methods: result is more or less independent
      #ftemp,Ltemp=np.linalg.eig(F[i])
      ftemp,Ltemp=np.linalg.eigh(F[i])
      #ftemp,Ltemp,info=dsyev(F[i]) #this seems to be the best function

      index=np.argsort(np.real(ftemp),kind='heapsort') # ascending sorting f
      f[i]=np.real(ftemp[index]).T[:].T[6:].T
      L[i]=(Ltemp.T[index].T)[:].T[6:].T

      #the frequencies are square root of the eigen values of F
      for j in range(len(f[i])):
         f[i][j]=np.sign(f[i][j])*np.sqrt(np.abs(f[i][j]))

      if np.any(f[i]<0):
         logging[1].write('imaginary frequencies occured. The absolute'
                        ' values are used in the following.\n{0}\n'.format(f[i]))
         f[i]=np.abs(f[i])
      if logging[0]<2:
         logging[1].write("Frequencies (cm-1)\n"+\
               repr(f[i]*Hartree2cm_1)+"\nL-matrix \n"+ repr(L[i])+"\n")

      M=np.eye(dim)
      for j in range(0,dim):
         M[j,j]/=mass[j//3]
      Lmass[i]=M.dot(L[i])
      #renormalise
      for j in range(0,dim-6):
         norm=np.sum(Lmass.T[j]*Lmass.T[j])
         if np.abs(norm)>1e-12:
            Lmass.T[j]/=np.sqrt(norm)
   return f, L, Lmass

def gradientHR(logging, initial, final, opt, HRthresh):
   """ This function gathers most essential parts for calculation of HR-factors from g09-files

   **PARAMETERS**
   logging: This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
            and 4- only main information) and logging[1] is the file, already opened, to write the information in.
   initial: name of the file containing the initial state's geometry
   final:   name of the file containing the initial state's geometry
   opt:     options that are given to this calculation; especially it is of interest, whether there should be frequencies and/or normal
            modes to be read from the g09-files.
   HRthresh threshold for HR-factors (only those larger than this are taken into account)

   **RETURNS**
   HR:      Huang-Rhys factors, sorted by size
   funi:    vibrational frequencies sorted same as HR 
   Energy:  Energy difference between the two states
   J:       Duschinsky-rotation matrix
   K:       shift between the states (in normal coordinates)
   f:       frequencies of the vibrational modes, sorted by size of the frequencies

   """
   assert len(initial)>0, 'no initial state found!'
   assert len(final)>0, 'no final state found!'
   initial=initial[0]
   final=final[0]
   assert os.path.isfile(initial) and os.access(initial, os.R_OK),\
            initial+' is not a valid file name or not readable.'
   assert os.path.isfile(final) and os.access(final, os.R_OK),\
            final+' is not a valid file name or not readable.'
   dim, Coord, mass, A, E=rl.ReadG09(logging, initial)
   F, CartCoord, P, Energy=quantity(logging, dim, 2 ) #creates respective quantities (empty)
   if logging[0]==0:
      logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2))
   F[0],Energy[0]=A, E
   F[1]=F[0] #force constant matrix in both states coincides
   Grad, E=rl.ReadG092(logging, final) 
   Energy[1]=E
   #read coordinates, force constant, binding energies from log-files and calculate needed quantities
   if logging[0]<3:
      logging[1].write('difference of minimum energy between states:'
                       ' Delta E= {0}\n'.format((Energy[0]-Energy[1])*Hartree2cm_1))
      if logging[0]<2:
         logging[1].write('initial state: \n{0}\n'.format(F[0]))

   #Calculate Frequencies and normal modes
   f, Lsorted, Lcart=GetL(logging, dim, mass,F)
   K, J=GradientShift(logging, Lsorted, mass, Grad, f[0])
   extra=re.findall(r"g09Vectors",opt, re.I)
   if extra!=[]:
      g09L=getGaussianL([initial], dim)
      rl.replace(logging, final, g09f[0], g09L)
   elif re.search(r"g09Vector",opt, re.I) is not None:
      g09L=getGaussianL([initial], dim)
      rl.replace(logging, final, g09f[0], g09L)
   
   #calculate HR-spect
   HR, funi= HuangR(logging, K, f, HRthresh)
   if (re.search(r"makeLog", opt, re.I) is not None) is True:  
      rl.replace(logging, initial, f[0], Lcart[0])
   return HR, funi, Energy, J, K, f

def GradientShift(logging, L, mass, Grad, Freq):
   """ This function calculates the 'shift' between excited state and ground state from the gradient of the excited state 
   at ground state geometry assuming coinciding frequencies and harmonic potentials.
   """
   #get dimensionality of the problem and initialise quantitios
   dim=len(mass)*3
   J=np.zeros(( len(L)-1,dim-6,dim-6 ))
   M=np.zeros((dim,dim)) 
   for j in range(0,dim):
      M[j,j]=1/mass[j//3]
   #K becomes now gradient in massweighted internal coordinates
   K=Grad.T.dot(M).dot(L[0])
   # scale consistently: Now it is really the shift in terms of normal modes
   K/=Freq*Freq*np.sqrt(2)  ##??? wtf
   #calculate Duschinsky-matrix
   for i in range(len(J)):
      J[i]=np.dot(L[0].T, np.linalg.pinv(L[i+1].T)) 
   return K, J

def gs(A):
   """This function does row-wise Gram-Schmidt orthonormalization of matrices. 
   code for Gram-Schmidt adapted from iizukak, see https://gist.github.com/iizukak/1287876
   """
   #def proj(v1, v2):
      #return map(lambda x : x *(np.dot(v2,v1) / np.dot(v1, v1)) , v1)

   X=A.T # I want to orthogonalize row-wise
   Y = []
   npdot=np.dot
   for i in range(len(X)):
      temp_vec = X[i]
      for inY in Y :
         #proj_vec = proj(inY, X[i])
         proj_vec = map(lambda x : x *(npdot(X[i],inY) / npdot(inY, inY)) , inY)
         temp_vec = map(lambda x, y : x - y, temp_vec, proj_vec)
      Y.append( temp_vec/np.linalg.norm(temp_vec)) # normalise vectors
   return np.matrix(Y).T # undo transposition in the beginning

def HuangR(logging, K, f, HRthresh): #what is with different frequencies???
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
   sortHR=np.zeros(len(K[0]))
   HR=np.zeros(len(K[0]))
   fsort=np.zeros(len(K[0]))
   uniHRall=[]
   uniFall=[]
   #HR=[K[0][j]*K[0][j]*0.5*f[0][j] for j in range(len(K[0]))]
   for j in range(len(K[0])):
      HR[j]=K[0][j]*K[0][j]*f[1][j]*.5
      #HR[j]=K[0][j]*K[0][j]*f[1][j]
   index=np.argsort(HR, kind='heapsort')
   sortHR=HR[index]
   fsort=f[1][index]
   if np.any(fsort)<0:
      logging[1].write('ATTENTION: some HR-factors are <0.\
               In the following their absolute value is used.')
      fsort=np.abs(fsort)
   uniHR=[]
   uniF=[]
   loggingwrite=logging[1].write

   loggingwrite(u'HR-fact           freq\n')
   #print(u'HR-fact           freq\n')
   for j in xrange(len(sortHR)):
      #select all 'big' HR-factors 
      if sortHR[-j]>=HRthresh:
         uniHR.append(sortHR[-j])
         uniF.append(fsort[-j])
         loggingwrite(u"%f   %f\n"%(sortHR[-j], fsort[-j]*Hartree2cm_1))
   uniHRall.append(uniHR)
   uniFall.append(uniF)
   return uniHRall, uniFall

def quantity(logging, dim, num_of_files):
   """ Here some frequencies are defined; it is just for clearer code.
   This function is called by ReadLog only.
   """
   F=np.zeros((num_of_files, dim, dim)) 
   CartCoord=np.zeros((num_of_files, 3, dim/3))
   P=np.zeros((num_of_files, dim,dim))
   Energy=np.zeros(num_of_files)
   return F, CartCoord, P, Energy

def ReadHR(logging, HRfile):
   """ This function reads the HR-factors and electronic transition energy from a given file and brings them into a 
   similar structure as they are used in the 'smallscript'.

   **PARAMETERS**
   logging:     This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
                and 4- only main information) and logging[1] is the file, already opened, to write the information in.
   HRfile:      the file where the information is found

   **RETURNS**
   initial:     a dummy-array that originally contains information about the inital states. Here at the moment only one
                is allowed and only its length is relevant in the further programme.
   HRm:         a 2-dimensional array containing all Huang-Rhys-factors
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

version=1.5
# End of functions_smsc.py
