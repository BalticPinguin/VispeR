#!/usr/bin/python2
import Spect

class FC-spect(Spect):
    """ All models that are based on FC-scheme are based on this class. 
    Since it doesn't need any more data than the DR-scheme, everything 
    initialised for Spect is initialized here already and hence no extra
    init()-function is needed.
    """

    def CalculationHR(self, initial, final, opt, HRthresh):
       """ This function gathers most essential parts for calculation of 
       HR-factors from g09-files. That is: read neccecary information from the  
       g09-files and calculate HR-factors as well as the  Duschinsky-rotation 
       matrix and the shift between minima (needed later if the option Duschinsky 
       is specified)
   
       **PARAMETERS**
       initial: name of the file containing the initial state's geometry
       final:   name of the file containing the initial state's geometry
       opt:     options that are given to this calculation; especially it is of 
                interest, whether there should be frequencies and/or normal
                modes to be read from the g09-files.
       HRthresh threshold for HR-factors (only those larger than this are taken 
                into account)
   
       **RETURNS**
       HR:      Huang-Rhys factors, sorted by size
       funi:    vibrational frequencies sorted same as HR 
       Energy:  Energy difference between the two states
       J:       Duschinsky-rotation matrix
       K:       shift between the states (in normal coordinates)
       f:       frequencies of the vibrational modes, sorted by size of the 
                frequencies
   
       """                                                                                             
       # if the input-files are G09-formcheck-files (recognised by ending):
       if (( ".fchk" in final) and (".fchk" in initial))\
             or (( ".FChk" in final) and (".FChk" in initial)): # this is also a valid ending
          dim, Coord, mass, A, E=rl.ReadGO9_fchk(logging, initial)
          F, CartCoord, P, Energy=quantity(logging, dim, 2) #creates respective quantities (empty)
          if logging[0]==0:
             logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2)+"\n")
          F[0],Energy[0]=A, E
          CartCoord[0]=Coord
          dim, CartCoord[1], mass, F[1], Energy[1]=rl.ReadGO9_fchk(logging, final) 
          #F[0]=F[1]
       #else, test what kind of file was given: G09, GAMESS or NWChem
       else:
          with open(initial, "r+b") as f: #open file as mmap
             mapping = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
             for line in iter(mapping.readline, ""): #go through every line and test characteristic part
                if "GAMESS" in line: #if it is found: read important quantities from file
                   dim, Coord, mass, A, E=rl.ReadGAMESS(logging, initial)
                   F, CartCoord, P, Energy=quantity(logging, dim, 2) #creates respective quantities (empty)
                   if logging[0]==0:
                      logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2)+"\n")
                   F[0],Energy[0]=A, E
                   CartCoord[0]=Coord
                   dim, Coord, mass, A, E=rl.ReadGAMESS(logging, final) 
                   F[1], Energy[1]=A, E
                   CartCoord[1]=Coord
                   break
                elif "Northwest Computational Chemistry Package (NWChem)" in line:
                   dim, Coord, mass, A, E=rl.ReadNWChem(logging, initial)
                   F, CartCoord, P, Energy=quantity(logging, dim, 2) #creates respective quantities (empty)
                   if logging[0]==0:
                      logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2)+"\n")
                   F[0],Energy[0]=A, E
                   CartCoord[0]=Coord
                   dim, CartCoord[1], mass, F[1], Energy[1]=rl.ReadNWChem(logging, final) 
                   break
                elif "Gaussian(R)" in line:
                   dim, Coord, mass, A, E=rl.ReadG09(logging, initial)
                   F, CartCoord, P, Energy=quantity(logging, dim, 2) #creates respective quantities (empty)
                   if logging[0]==0:
                      logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2)+"\n")
                   F[0],Energy[0]=A, E
                   CartCoord[0]=Coord
                   dim, Coord, mass, A, E=rl.ReadG09(logging, final) 
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
       J, K=Duschinsky(logging, Lmassw, mass, dim, CartCoord)
       #J, K=Duschinsky(logging, Lsorted, mass, dim, CartCoord)
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
       print opt
       if (re.search(r"makeLog", opt, re.I) is not None) is True:  
          rl.replace(logging, "S0.log", f[0], Lmassw[0])
       #Hartree2cm_1=219474.63 
       return HR, funi, Energy, J, K, f

    def calcspect(self, HR, freq, E, E0, N, M, T):
       """This is used to calculate the line spectrum assuming no mode mixing 
       (shift only) 
       and coinciding frequencies in both electronic states.
    
       **PARAMETERS:**
       logging:This variable consists of two parts: logging[0] specifies the level
               of print-out (which is between 0- very detailed
               and 4- only main information) and logging[1] is the file, already 
               opened, to write the information in.
       HR:     Huang-Rhys factors
       n:      number of modes that are considered here (with largest HR)
       freq:   frequencies (have to be in the same order as HR
       E:      energy difference of energy surfaces
       N,M:    are the numbers of vibrational quanta can be in the modes
    
       **RETURNS:**
       nothing (output into /tmp/linspect)
       """
       # M,N are maximal numbers of vibrational modes 
       #     (+1, since they should be arrived really; count from 0)
       N+=1
       M+=1
       def FCeqf( Deltag, M, N):
          """Calculate Franck-Condon factors under assumption of equal frequencies 
          for only one vibrational mode
    
          *PARAMETERS:*
          Deltag: HR-factor of respective state
          N:      excitation number of initial state
          M:      excitation number of final state
    
          *RETURNS:*
          Franck-Condon factor of the respective transition
          """
          fact=math.factorial
          faktNM=fact(M)*fact(N)
          FC=0
          for x in range(int(min(N,M))+1):# +1: to go to really reach min(M,N).
             FC+=math.pow(-1,N-x)*math.pow(Deltag,(M+N)*0.5-x)/(fact(M-x)*fact(N-x))*\
                   math.sqrt(faktNM)/fact(x)
          return FC*FC
       
       def unifSpect(intens, freqs, E, FC00):
          """ Calculation of the line spectrum respecting only shift of minima 
          (no Duschinsky rotation, no change if frequency) 
          and assuming coinciding frequencies for initial and final state
    
          **PARAMETERS:**
          intens: matrix of intensities of transitions
          freqs:  matrix of respective energies
    
          **RETURNS**
          a 3-dimensional array with energies of transition (1. column), their 
          rate (2. column) and the number of changing mode.
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
       setM=False
       if M==1: 
          # there was not specified, how many vibr. states in ground-state 
          #     should be taken into account
          setM=True
          M=max(3,int(-1.1*HR[0]*HR[0]+6.4*HR[0]+9.))
          #i.e. take M following an empirical value as function of the HR-file
       assert n>0, "There is no Huang-Rhys factor larger than the respective"+\
                                        "threshold. No mode to be calculated."
       #if setM: the size of these arrays will be overestimated.
       FC=np.zeros((n,M*N-1))
       uency=np.zeros((n,M*N-1)) #frequency
       #calculate 0->0 transition
       FC00=1.00
       #print 0,0,0, 10
       uency00=E*Hartree2cm_1 #zero-zero transition
       loggingwrite=logging[1].write #avoid dots!
       npexp=np.exp #avoiding dots accelerates python quite a lot
       #here a goes over all modes
       sgnE=np.sign(E)
       if sgnE==0:
          # sign(0) =0 ...
          sgnE=1;
       for a in xrange(n):
          #print a, HR[a]
          for j in range(N):  # initial state
             if setM:
                # set M to fit best to the value at that moment.
                M=max(3,int(-1.1*HR[a]*HR[a]+6.4*HR[a]+9.))
             for i in range(M):  #final states
                if i==0 and j==0:
                   ##skip 0-0 transitions
                   continue
                tmp=FCeqf(HR[a], i, j)
                try:
                   FC[a][j*M+i-1]=tmp*FC00*npexp(-(E0+freq[a]*j)/T)
                   uency[a][j*M+i-1]=(E+sgnE*freq[a]*(j-i))*Hartree2cm_1
                except IndexError:
                   logging[1].write("truncated spectrum for mode nr. %d"%(a))
                   break
       FC00*=npexp(-E0/T)
       spect=unifSpect(FC, uency, E*Hartree2cm_1, FC00)
       return spect

    def gs(A):
       """This function does row-wise Gram-Schmidt orthonormalization of 
       matrices.
       code for Gram-Schmidt adapted from iizukak, 
       see https://gist.github.com/iizukak/1287876
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
             proj_vec = map(lambda x : x *(npdot(X[i],inY) / 
                 npdot(inY, inY)) , inY)
             temp_vec = map(lambda x, y : x - y, temp_vec, proj_vec)
          Y.append( temp_vec/np.linalg.norm(temp_vec)) # normalise vectors
       return np.matrix(Y).T # undo transposition in the beginning
    
    def HuangR(self, K, f, HRthresh): #what is with different frequencies???
       """ Function that calculates the Huang-Rhys factors for all 
       vibrational states
    
       **PARAMETERS**
       logging: This variable consists of two parts: logging[0] specifies the 
                level of print-out (which is between 0- very detailed
                and 4- only main information) and logging[1] is the file, already 
                opened, to write the information in.
       K        The displacements of minima in internal coordinates (n-vector)
       f        The frequencies of respective modes (2xn-matrix)
       HRthresh threshold for lowest HR-factor (number)
    
       **RETURNS**
       uniHRall   HR-factors for coinciding frequencies sorted by size (decreasing)
       uniFall    respective frequencies (same order)
       """
       lenK=len(K)
       sortHR=np.zeros(lenK)
       HR=np.zeros(lenK)
       fsort=np.zeros(lenK)
       uniHRall=[]
       uniFall=[]
       #print K
       for j in range(lenK):
          HR[j]=K[j]*K[j]*f[0][j]*.5 
       index=np.argsort(HR, kind='heapsort')
       sortHR=HR[index]
       fsort0=f[0][index]
       fsort1=f[1][index]
       if np.any(fsort)<0:
          logging[1].write('ATTENTION: some HR-factors are <0.\
                   In the following their absolute value is used.')
          fsort1=np.abs(fsort1)
          fsort0=np.abs(fsort0)
       uniHR=[]
       uniF1=[]
       uniF0=[]
       loggingwrite=logging[1].write
       if sortHR[-1]>=10: 
          #if larges HR-factor is too large
          loggingwrite(u'\n!! ATTENTION!! THE HUANG-RHYS FACTOR SEEMS TO BE'+\
                                                            ' TOO LARGE !!\n')
          loggingwrite(u'the spectrum will be calculated, but most probably'+\
                                         ' the input-stat is inconsistent.\n')
       loggingwrite(u'HR-fact           freq     delta\n')
       #print(u'HR-fact           freq\n')
       for j in range(len(sortHR)-1,-1,-1):
          #select all 'big' HR-factors 
          if sortHR[j]>=HRthresh:
             uniHR.append(sortHR[j])
             uniF1.append(fsort1[j])
             uniF0.append(fsort0[j])
             loggingwrite(u"%f   %f   %f\n"%(sortHR[j], fsort1[j]*Hartree2cm_1, 
                                                np.sqrt(fsort1[j]/fsort0[j]) ))
          else:
             # there will come only smaller ones.
             break
       uniHRall.append(uniHR)
       #keep order: First one is initial, second is final state.
       uniFall.append(uniF0)
       uniFall.append(uniF1)
       return uniHRall, uniFall
    
    def quantity(self, dim, num_of_files):
       """ Here some frequencies are defined; it is just for clearer code.
       This function is called by CalculationHR.
       
       **PARAMETERS**
       logging: This variable consists of two parts: logging[0] specifies the 
                level of print-out (which is between 0- very detailed
                and 4- only main information) and logging[1] is the file, already 
                opened, to write the information in.
       dim      dimension of matrices/vectors
       num_of_files  second dimension; this is always 1
       """
       F=np.zeros((num_of_files, dim, dim)) 
       CartCoord=np.zeros((num_of_files, 3, dim//3))
       P=np.zeros((num_of_files, dim,dim))
       Energy=np.zeros(num_of_files)
       return F, CartCoord, P, Energy
    
class CFC-spect(FC-spect):
    def calcspect(self, HR, freq, E, E0, N, M, T):
       """This is used to calculate the line spectrum assuming no mode mixing 
       (shift only)  and coinciding frequencies in both electronic states.
   
       **PARAMETERS:**
       logging:This variable consists of two parts: logging[0] specifies the 
           level of print-out (which is between 0- very detailed
           and 4- only main information) and logging[1] is the file, 
           already opened, to write the information in.
       HR:     Huang-Rhys factors
       n:      number of modes that are considered here (with biggest HR)
       freq:   frequencies (have to be in the same order as HR
       E:      energy difference of energy surfaces
       N,M:    are the numbers of vibrational quanta can be in the modes
       All arguments are neccesary.
   
       **RETURNS:**
       nothing (output into /tmp/linspect)
       """
       # M,N are maximal numbers of vibrational modes (+1, since they should 
       #      be arrived really; count from 0)
       N+=1
       M+=1
   
       def FCchf(HR,N,M,freq):
          npsqrt=np.sqrt
          D=npsqrt(2.*HR) # define this to become consistent with given formula
          delta=npsqrt(freq[1]/freq[0])
          deltasquare=delta*delta
          #R00=sqrt(2*delta/(1+delta*delta))*np.exp(-.5*D*D/(1+delta*delta))
          R=np.zeros( (M,N) )
          # R_00 is normalisation -> actually calculate FC_ij/FC_00 --> 
          #                     than I can chose FC_00 outside.
          R[0][0]=1.00
          R[0][1]=-npsqrt(2.)*delta*D/(1.+deltasquare) 
          R[1][0]=npsqrt(2.)*D/(1.+deltasquare)
          R[1][1]=(D*R[0][1]+delta*npsqrt(2.)*R[0][0])*npsqrt(2.)/(1+deltasquare)
          for j in range(2,N):
             R[0][j]=(-2.*D*delta*R[0][j-1]-
                      (deltasquare-1)*npsqrt(2.*(j-1))*R[0][j-2])\
                      /(npsqrt(2.*j)*(1+deltasquare))
             R[1][j]=(2.*delta*npsqrt(2)*R[0][j-1]-2*D*delta*R[1][j-1]-\
                      (deltasquare-1)*npsqrt(2*(j-1))*R[1][j-2])\
                      /(npsqrt(2.*j)*(1.+deltasquare))
          for j in range(2,M):
             R[j][0]=(npsqrt(2.*(j-1))*(deltasquare-1)*R[j-2][0]+2.*D*R[j-1][0])\
                      /(npsqrt(2.*j)*(1.+deltasquare))
          for i in range(2,M):
             for j in range(1,N):
                R[i][j]=((deltasquare-1)*npsqrt(2.*(i-1))*R[i-2][j]+\
                      2.*D*R[i-1][j]+2.*delta*npsqrt(2.*j)*R[i-1][j-1] )\
                      /(npsqrt(2.*i)*(1.+deltasquare))
          return R # this is matrix with square roots of transition probabilities
   
       def unifSpect(intens, freqs, E, FC00):
          """ Calculation of the line spectrum respecting only shift of minima 
          (no Duschinsky rotation) 
          and assuming coinciding frequencies for initial and final state
   
          **PARAMETERS:**
          intens: matrix of intensities of transitions
          freqs:  matrix of respective energies
   
          **RETURNS**
          a 2-dimensional array with energies of transition (1. column) and their 
          rate (2. column)
          """
          J=len(intens[0])
          spect=np.zeros((3,len(intens)*J+1))
          #first transition: 0-0 transition. 
          spect[1][0]=FC00 #0->0 transition
          spect[0][0]=E
          spect[2][0]=0
          for i in range(len(intens)):
             for j in xrange(J):
                spect[2][i*J+j+1]=i+1
                spect[1][i*J+j+1]=intens[i][j]
                spect[0][i*J+j+1]=freqs[i][j]
          return spect
   
       n=len(HR) #=len(freq)
       setM=False
   
       # correct for vibrational groundstates:
       #E+=(sum(freq[0])-sum(freq[1]))*.5
       #print E, freq[0]
       if M==1: 
          # there was not specified, how many vibr. states in ground-state 
          #           should be taken into account
          setM=True
          M=max(3,int(-1.1*HR[0]*HR[0]+6.4*HR[0]+9.))
       assert n>0, "There is no Huang-Rhys factor larger than the respective"+\
                    " threshold. No mode to be calculated."
       #if setM: the size of these arrays will be overestimated.
       #calculate 0->0 transition
       FC00=1
       #print 0,0,0, 10
       uency00=E*Hartree2cm_1 #zero-zero transition
       loggingwrite=logging[1].write #avoid dots!
       npexp=np.exp                  #to accelerates python quite a lot
       FC=np.zeros((n,M*N))
       uency=np.zeros((n,M*N)) #frequency
       #here a goes over all modes
       sgnE=np.sign(E)
       if np.sign(E)==0:
          sgnE=1
       for a in xrange(n):
          if setM:
             # set M to fit best to the value at that moment.
             M=max(3,int(-1.1*HR[a]*HR[a]+6.4*HR[a]+9.))
          #print a, HR[a]
          R=FCchf(HR[a],N,M,[freq[0][a], freq[1][a]])
          for j in range(N): 
             for i in range(M):
                if i==0 and j==0:
                   ##skip 0-0 transitions
                   continue
                tmp=R[i][j]*R[i][j]
              #  if tmp*npexp(-(E0+freq[0][a]*j)/T)>0.01:
              #     print a,i,j, tmp*npexp(-(E0+freq[0][a]*j)/T),\
              #           (freq[0][a]*j-freq[1][a]*i)*Hartree2cm_1
                try:
                   FC[a][j*M+i-1]=tmp*FC00*npexp(-(E0+freq[0][a]*j)/T)
                   uency[a][j*M+i-1]=(E+sgnE*(freq[0][a]*j-freq[1][a]*i))*Hartree2cm_1
                except IndexError:
                   logging[1].write("truncated spectrum for mode nr. %d"%(a))
                   break
       FC00*=npexp(-E0/T) # rescale 0-0 transition due to thermal population
       spect=unifSpect(FC, uency, E*Hartree2cm_1, FC00)
       return spect

class GFC-spect(FC-spect):
   def gradientHR(self, initial, final, opt, HRthresh):
      """ This function gathers most essential parts for calculation of 
          HR-factors from g09-files
   
      **PARAMETERS**
      logging: This variable consists of two parts: logging[0] specifies the 
               level of print-out (which is between 0- very detailed
               and 4- only main information) and logging[1] is the file, already 
               opened, to write the information in.
      initial: name of the file containing the initial state's geometry
      final:   name of the file containing the initial state's geometry
      opt:     options that are given to this calculation; especially it is of 
               interest, whether there should be frequencies and/or normal
               modes to be read from the g09-files.
      HRthresh threshold for HR-factors (only those larger than this are taken 
               into account)
   
      **RETURNS**
      HR:      Huang-Rhys factors, sorted by size
      funi:    vibrational frequencies sorted same as HR 
      Energy:  Energy difference between the two states
      J:       Duschinsky-rotation matrix
      K:       shift between the states (in normal coordinates)
      f:       frequencies of the vibrational modes, sorted by size of the 
               frequencies
   
      """
      assert len(initial)>0, 'no initial state found!'
      assert len(final)>0, 'no final state found!'
      initial=initial[0]
      final=final[0]
      assert os.path.isfile(initial) and os.access(initial, os.R_OK),\
               initial+' is not a valid file name or not readable.'
      assert os.path.isfile(final) and os.access(final, os.R_OK),\
               final+' is not a valid file name or not readable.'
      # if the input-files are G09-formcheck-files (recognised by ending):
      if (( ".fchk" in final) and (".fchk" in initial))\
            or (( ".FChk" in final) and (".FChk" in initial)): # this is also a valid ending
         dim, Coord, mass, A, E=rl.ReadGO9_fchk(logging, initial)
         F, CartCoord, P, Energy=quantity(logging, dim, 2) #creates respective quantities (empty)
         if logging[0]==0:
            logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2)+"\n")
         F[0],Energy[0]=A, E
         CartCoord[0]=Coord
         F[1]=F[0] #force constant matrix in both states coincides
         Grad, E=rl.ReadG09_fchk2(logging, final) 
         Energy[1]=E
      #else, test what kind of file was given: G09, GAMESS or NWChem
      else:
         with open(initial, "r+b") as f: #open file as mmap
            mapping = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
            for line in iter(mapping.readline, ""): #go through every line and test characteristic part
               if "GAMESS" in line: #if it is found: read important quantities from file
                  assert 1==0, "gradient method for GAMESS-files not "+\
                                  "implemented yet."
                  dim, Coord, mass, A, E=rl.ReadGAMESS(logging, initial)
                  F, CartCoord, P, Energy=quantity(logging, dim, 2) #creates respective quantities (empty)
                  if logging[0]==0:
                     logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2)+"\n")
                  F[0],Energy[0]=A, E
                  CartCoord[0]=Coord
                  # --> here read gradient!!!
                  break
               elif "Northwest Computational Chemistry Package (NWChem)" in line:
                  assert 1==0, "gradient method for nwchem-files not "+\
                                  "implemented yet."
                  dim, Coord, mass, A, E=rl.ReadNWChem(logging, initial)
                  F, CartCoord, P, Energy=quantity(logging, dim, 2) #creates respective quantities (empty)
                  if logging[0]==0:
                     logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2)+"\n")
                  F[0],Energy[0]=A, E
                  CartCoord[0]=Coord
                  # --> here read gradient!!!
                  break
               elif "Gaussian(R)" in line:
                  dim, Coord, mass, A, E=rl.ReadG09(logging, initial)
                  F, CartCoord, P, Energy=quantity(logging, dim, 2) #creates respective quantities (empty)
                  if logging[0]==0:
                     logging[1].write("Dimensions: "+ str(dim)+ '\n Masses: '+ str(mass**2)+"\n")
                  F[0],Energy[0]=A, E
                  CartCoord[0]=Coord
                  dim, Coord, mass, A, E=rl.ReadG09(logging, final) 
                  F[1]=F[0] #force constant matrix in both states coincides
                  Grad, E=rl.ReadG092(logging, final) 
                  Energy[1]=E
                  break
            else: # is there some error??  probably this message is printed also if files were recognised
               print "file type not recognised"
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
   
   def GradientShift(self, L, mass, Grad, Freq):
      """ This function calculates the 'shift' between excited state and ground 
          state from the gradient of the excited state  at ground state geometry 
          assuming coinciding frequencies and harmonic potentials.
      **PARAMETERS**
      logging: This variable consists of two parts: logging[0] specifies the 
               level of print-out (which is between 0- very detailed
               and 4- only main information) and logging[1] is the file, already 
               opened, to write the information in.
      L        Eigen-vectors of force-constant matrix (matrix)
      mass     vector of masses (in a.u.)
      Grad     gradient vector of final states PES at ground state equilibr.
      Freq     (2xn -matrix) containing frequencies in initial/final state
      """
      #get dimensionality of the problem and initialise quantitios
      dim=len(mass)*3
      J=np.zeros(( len(L)-1,dim-6,dim-6 ))
      M=np.zeros((dim,dim)) 
      for j in range(0,dim):
         M[j,j]=1/mass[j//3]
      #K becomes now gradient in massweighted internal coordinates
      K=Grad.T.dot(M).dot(L[0])[0]
      # scale consistently: Now it is really the shift in terms of normal modes
      K/=Freq*Freq*np.sqrt(2)  ##
      
      #K*=np.sqrt(np.pi)/2. #correction factor due to some magic reason
      #calculate Duschinsky-matrix
      J=np.dot(L[0].T, L[1])
   
      if logging[0]<2:
         logging[1].write('Duschinsky rotation matrix:\n')
         k=range(0,dim-6)
         s=0
         t=min(s+5,dim-6)
         while s<dim-6:
            for temp in range(s,t):
               logging[1].write("               %d "%(k[temp]+1))
            logging[1].write("\n")
            for j in range(len(J)):
               logging[1].write(" %03d"%(j+1))
               for temp in range(s,t):
                  logging[1].write("   %+.5e"%(J[j][k[temp]]))
               logging[1].write("\n")
            s=t
            t=min(s+5,dim-6)
         logging[1].write('\nDuschinsky displacement vector:\n')
   
         for j in range(len(K)):
            logging[1].write("  %d    %e\n"%(j+1, K[j]))
      return K, J
   
   
class HR-spect(FC-spect):
      """ First, I will leave this class empty. Does one need it? 
      From its structure, this is not as the other spectra and hence there is not much to
      inherit from; maybe later this can be added to some structure.
      """
   def ReadHR(self, HRfile):
      """ This function reads the HR-factors and electronic transition energy 
      from a given file and brings them into a 
      similar structure as they are used in the 'smallscript'.
   
      **PARAMETERS**
      logging:     This variable consists of two parts: logging[0] specifies 
                   the level of print-out (which is between 0- very detailed
                   and 4- only main information) and logging[1] is the file, 
                   already opened, to write the information in.
      HRfile:      the file where the information is found
   
      **RETURNS**
      initial:     a dummy-array that originally contains information about 
                   the inital states. Here at the moment only one
                   is allowed and only its length is relevant in the further 
                   programme.
      HRm:         a 2-dimensional array containing all Huang-Rhys-factors
      freqm:       analogously to HRm, containing the respective frequencies
      Energy:      the energy-difference of the electronic states. 
                   (in atomic units)
   
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
      HRfreq=re.findall(r"HR-fact[\s]*freq[\s]*\n[\n\d\.\se\+\-]*", f, re.I)
      assert len(HRfreq)==1, "The file-format could not be read. exit now"
      HRf=re.findall(r"(?<=\n)[\d\.]*[\s]+[\d\.]*", HRfreq[0], re.I)
      for i in range(len(HRf)):
         line=re.findall(r"[\d.]+",HRf[i], re.I)
         HR.append(float(line[0]))
         funi.append(float(line[1])/Hartree2cm_1)
      initial=['excited']
      #the following is just to be consistent with structure of 
      #                         HR calculated in first part
      HRm=np.zeros((1,len(HR)))
      HRm[0]=HR
      freqm=np.zeros((2,len(HR)))
      freqm[0]=funi
      freqm[1]=funi # no changing frequencies
      return initial, HRm, freqm, Energy

class HR-factors():
      """ This is not a spectrum. Hence it doesn't inherit from Spect.
      First, I will leave this class empty. Does one need it?
      """
