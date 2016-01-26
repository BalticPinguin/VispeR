#!/usr/bin/python2
# filename: Spect.py
import numpy as np
from scipy.linalg.lapack import dsyev
import re, mmap, os.path, math, sys

AMU2au=1822.88839                                          
Angs2Bohr=1/0.52917721092                                  
Hartree2GHz=6.579684e6                                     
Hartree2cm_1=219474.63 

class  Spect:
   """ This is the general class that all spectra belong to. It contains the fundamental 
       quantities that are required and most of the functions to use. 
       Some of the functions are redefined in the respective sub-classes.
   """
   # Below are the conversion factors and fundamental constant
   AMU2au=1822.88839                                          
   Angs2Bohr=1/0.52917721092                                  
   Hartree2GHz=6.579684e6                                     
   Hartree2cm_1=219474.63 
   
   # BEGIN OF DATA-DEF.
   Energy
   dim
   CartCoord
   logging
   # END OF DATA-DEF.

   ## BEGIN OF FUNCTION DEFINITIONS
   def __init__(self, f):
      """ This function initialises the Spect-object and creates its first ... .
      The INPUT, 'f', is the input-file to smallscript. In this function,
      only the information for output is read and the logging object,
      which has the output-file and the level of output-information is defined.
      All variables/functions that are common to all spectral tyes are initialised here."""
      
      def invokeLogging(logfile, mode="important"):
         """ initialises the logging-functionality
         ==PARAMETERS==
         logfile      name of file to be used as log-file. It is expected to be an array of 
                      length 0 or one.
         mode:        5 different values are possible (see below); the lowest means: print 
                      much, higher values mean 
                      less printing
         ==RETURNS==
         logging:     number of respective mode
         log:         opened file to write in
         """
         if logfile==[]:
            log=open("calculation.log", "a")
         else:
            s=logfile[0].strip()
            log=open(s, "a")
  
         if mode in ['all', 'ALL', 'All', '0']:
            logging=0
            log.write('use log-level all\n')
         elif mode in ['detailed', 'DETAILED', 'Detailed', "1"]:
            logging=1
            log.write('use log-level detailed\n')
         elif mode in ['medium', 'MEDIUM','Medium', '2']:
            logging=2
            log.write('use log-level medium\n')
         elif mode in ['important', 'IMPORTANT','Important', '3']:
            logging=3
         elif mode in ['short', 'SHORT','Short', '4']:
            logging=4
         else:
            logging=3
            log.write("logging-mode not recognized. Using 'important' instead\n")
         return logging, log
      
      def ReadData(self, initial, final, opt):
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
          """                                                                                             
         def quantity(self, dim):
            """ Here some frequencies are defined; it is just for clearer code.
            This function is called by CalculationHR.
            
            **PARAMETERS**
            logging: This variable consists of two parts: logging[0] specifies the 
                     level of print-out (which is between 0- very detailed
                     and 4- only main information) and logging[1] is the file, already 
                     opened, to write the information in.
            dim      dimension of matrices/vectors
            """
            F=np.zeros((2, dim, dim)) 
            self.CartCoord=np.zeros((2, 3, dim//3))
            P=np.zeros((2, dim,dim))
            return F, P

          # if the input-files are G09-formcheck-files (recognised by ending):
          if (( ".fchk" in final) and (".fchk" in initial))\
                or (( ".FChk" in final) and (".FChk" in initial)): # this is also a valid ending
             read = "G09_fchk"
          else:
             with open(initial, "r+b") as f: #open file as mmap
                mapping = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
                for line in iter(mapping.readline, ""): #go through every line and test characteristic part
                   if "GAMESS" in line: #if it is found: read important quantities from file
                      read ="GAMESS"
                      break
                   elif "Northwest Computational Chemistry Package (NWChem)" in line:
                      read ="NWChem"
                      break
                   elif "Gaussian(R)" in line:
                      read ="G09"
                      break
                #  There is some error: File not recognised or an other programme was used.
                else: 
                   print "file type not recognised"
          
          #read coordinates, force constant, binding energies from log-files and 
          # from the file, using the type of file that is known now...
          self.dim, Coord, mass, A, self.Energy[0]=Read(initial, read)
          F, CartCoord, P, Energy=quantity(dim) #creates respective quantities (empty)
          if self.logging[0]==0:
             self.logging[1].write("Dimensions: "+ str(self.dim)+ '\n Masses: '+ str(mass**2)+"\n")
          F[0]=A
          self.CartCoord[0]=Coord
          dim, self.CartCoord[1], mass, self.F[1], self.Energy[1]=Read(final,read) 

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
    

      #START LOG-FILE DEFINITION     
      #invoke logging (write results into file specified by 'out: ' or into 'calculation.log')
      logfile=re.findall(r"(?<=out:)[\.\-_ \w]+",f, re.I)[-1]
      loglevel=re.findall(r"(?<=print:)[\w]+",opt, re.I)[-1]
      self.logging = invokeLogging(logfile, loglevel )
      
      # now,  write the header to the output-file.
      self.logging[1].write("\n==================================================================\n"
               "===================  output of smallscript  ======================\n"
               "==================================================================\n\n")
      self.logging[1].write("   INPUT-FILE:\n")
      self.logging[1].write(f)
      self.logging[1].write(" \n   END OF INPUT-FILE \n\n")
      self.logging[1].write("calculations to be done: %s\n"%(todo))
      #END LOG-FILE DEFINITION     
      
      Energy=np.matrix(0,0)
      #START READING DATA FROM FILE
      # get files with electronic states:
      final=re.findall(r"(?<=final: )[\w.\-]+",f, re.I)
      initial=re.findall(r"(?<=initial: )[\w.\-]+",f, re.I)
      assert len(initial)==1, 'there must be one initial state'
      assert len(final)==1, 'there must be one final state'
      initial=initial[0]
      final=final[0]
      #check, if they are valid files and through an error if not.
      assert os.path.isfile(initial) and os.access(initial, os.R_OK),\
               initial+' is not a valid file name or not readable.'
      assert os.path.isfile(final) and os.access(final, os.R_OK),\
               final+' is not a valid file name or not readable.'
      #read options from input-file:
      opt=re.findall(r"(?<=opt:)[\w.,\(\) \=;:-]", f, re.I)
      ReadData(initial, final, opt)
      #END READING DATA FROM FILE
 
   def GetL(self, dim, mass, F):
      """ Function that calculates the frequencies and normal modes from force 
      constant matrix with and without projection onto internal degrees of freedom
   
      **argumets**
      logging: This variable consists of two parts: logging[0] specifies the 
               level of print-out (which is between 0- very detailed
               and 4- only main information) and logging[1] is the file, already 
               opened, to write the information in.
      dim      The dimensions of force-constant matrix
      mass     square-root of masses dim/3-dimensional array
      F        force-constant matrix
   
      **return**
      return f, Lsorted, Lcart
      f        two vertors (f[0] and f[1]) with the frequencies
      Lsorted  matrix of vibr. normal modes (transformation matrix)
      Lcart    massweighted L for comparison with the respective matrix from 
               the g09 log-file
      """
      # Defining arrays
      lenF=len(F[0])
      L=np.zeros(( len(F), lenF, lenF-6 )) 
      Lmass=np.zeros(( len(F), lenF, lenF-6 ))
      f=np.zeros(( len(F), lenF-6 ))
      
      def SortL0(J,L,f):
         np.set_printoptions(threshold=np.nan,linewidth=200)
         print "J\n", J
         resort=np.zeros(np.shape(J))
         roundJ=np.abs(np.around(J)) # don't correct for sign-flips
        # for i in range(len(J)):
        #    j=np.argmax(J[i])
        #    k=np.argmin(J[i])
        #    if J[i][j]>-J[i][k]:
        #       resort[i][j]=1
        #    else:
        #       resort[i][k]=-1
         for i in range(len(J)):
            print i,'\n', roundJ
            if np.sum(roundJ[i])==0: # there is no large element in the row...
               # insert an element such that it does not coincide with elements from before...
               gotit=False
               index=range(len(J))
               while not gotit:
                  j=np.argmax(np.abs(J[i][index]))
                  gotit=True
                  for k in range(i):
                     if np.sum(roundJ[i]*roundJ[k])>0: # i.e. there was a 1 before
                        index.delete(k)
                        gotit=False
                        continue # try next element
               #I found a working index
               assert delete!=[], "This is a nasty matrix. Your system has a very bad configuration and I have no clue how to solve it!"
               roundJ[i][index[0]]=1 # set the missing element
   
            elif np.sum(roundJ[i])>=2: # there are multiple large elements in the row... 
               # remove all except the largest one
               index=np.where(roundJ[i]==1)
               print index
               j=np.argmax(np.abs(J[i][index]))
               roundJ[i,index]=0
               roundJ[i,j]=1
            #if np.sum(roundJ[i])==1:  -> this is always true now.
            assert np.sum(roundJ[i])==1, "dumb Hubert. Don't do such stupid things!"
            #else: # there is exactly one element in the row
            j=np.argmax(roundJ[i]) # J(i,j) is one.
            index=np.where(roundJ[:,j]==1)
            if len(index)==1: # everything is fine
               continue
            else: # there are several elements with j-th line being 1
               #get the largest elemnt in this line
               print "index:" , index
               k=np.argmax(np.abs(J[index,j]))
               roundJ[index,j]=0
               roundJ[k,j]=1
               if k>i: # I deleted the row I am looking at
                  # insert an element such that it does not conflict with those inserted before
                  gotit=False
                  index=range(len(J))
                  while not gotit:
                     l=np.argmax(np.abs(J[index,j]))
                     gotit=True
                     for k in range(i):
                        if np.sum(roundJ[i]*roundJ[k])>0: # i.e. there was a 1 before
                           index.delete(k)
                           gotit=False
                           continue # try next element
               if k<i:
                  assert 1==0, "I made a mistake. This should never be true."
   
         resort=roundJ      
         #print "{0}".format(resort)
         invresort=np.linalg.inv(resort)
         #print "J\n", invresort.dot(J).dot(resort)
         return np.abs(resort.dot(f[1])), invresort.dot(L[1]).dot(resort)
      
      def SortL(J,L,f):
         resort=np.zeros(np.shape(J))
         #chose largest elements in lines
         for i in range(len(J)):
            j=np.argmax(J[i])
            k=np.argmin(J[i])
            if J[i][j]>-J[i][k]:
               resort[i][j]=1
            else:
               resort[i][k]=1
         #now, go through rows and check if they are ok:
         #print "resortJ\n",resort
         resort=resort.T
         Nos=[]
         freePlaces=[]
         for i in range(len(J)):
            if sum(resort[i])==1:
               continue
            elif sum(resort[i])==0:
               Nos.append(i)
            else:
               index=np.where(resort[i]==1)
               x=np.argmax(np.abs(J[index,i]))
               index=np.delete(index,x) # does it work the way I want it to work?
               resort[i][index]=0 #only x remains
               freePlaces=np.append(freePlaces,index)
         #print "resortJ\n",resort.T
         assert len(Nos)==len(freePlaces), "dodododo!"
         freePlaces=np.array(freePlaces,dtype=int)
         #print(freePlaces,"\n",Nos)
         #fill remaining lines.
         for i in range(len(Nos)):
               x=np.argmax(np.abs(J[freePlaces,Nos[i]]))
               resort[Nos[i],freePlaces[x]]=1 #only x remains
               freePlaces=np.delete(freePlaces,x) # does it work the way I want it to work?
         #print freePlaces
         assert len(freePlaces)==0, "the matrix is not readily processed."
         resort=resort.T
         #print "resortJ\n",resort
         invresort=np.linalg.inv(resort)
         assert np.all(invresort==resort.T), "The matrix is total bullshit!"
         #print "J\n", invresort.dot(J).dot(resort)
         return f.dot(invresort), L.dot(invresort)
         #  END OF SortL
   
      for i in range(len(F)):
         # here one can choose between the methods: result is more or less 
         #  independent
         #ftemp,Ltemp=np.linalg.eig(F[i])
         #ftemp,Ltemp=np.linalg.eigh(F[i])
         ftemp,Ltemp,info=dsyev(F[i]) #this seems to be the best function
         
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
         M=np.eye(dim)
         for j in range(dim):
            M[j,j]/=mass[j//3]
         Lmass[i]=M.dot(L[i])
         #if logging[0]<2:
          #  logging[1].write("Frequencies (cm-1)\n"+\
           #       repr(f[i]*Hartree2cm_1)+"\nL-matrix \n"+ repr(Lmass[i])+"\n")
      #print "unity\n", L[0].T.dot(L[0]), "\n", np.linalg.norm(L[0].T.dot(L[0]))
      #print "unity\n", L[1].T.dot(L[1]), "\n", np.linalg.norm(L[1].T.dot(L[1]))
      np.set_printoptions(threshold=np.nan,linewidth=500, suppress=True)
      #print "L_1\n",L[1]
      #print "test:\n", L[1].dot((np.linalg.pinv(L[0]).dot(L[1])).T)-L[0]
      f[1],L[1]=SortL(np.linalg.pinv(L[0]).dot(L[1]),L[1],f[1])
      #print "L_1\n",L[1]
      #print "test:\n", np.linalg.pinv(L[0]).dot(L[1])
      Lmass[1]=M.dot(L[1]) # recalculate Lmass!
      return f, L, Lmass
   
   def Duschinsky(self, L, mass, dim, x):
      """
      This function calculates the shift between two electronic states 
      (whose geometry is known, see x) as well as the
      Duschinsky-rotation matrix.
      **PARAMETERS:**
      logging: This variable consists of two parts: logging[0] specifies the 
               level of print-out (which is between 0- very detailed
               and 4- only main information) and logging[1] is the file, already 
               opened, to write the information in.
      L:       Matrix having mass-weighted normal modes as column-vectors,
               L[0] refers to initial state, L[1] to final one.
      mass:    array of square-roots of nuclear masses (length: N)
      dim:     dimensionality of the problem: 3*N
      x:       cartesian coordinates of the states of interest
   
      **RETURN:**
      J:    Duschinsky-rotation matrix
      K:    displacement-vector of energy-minima in normal coordinates
      """
      J=np.zeros((dim-6,dim-6))
      K=np.zeros(dim-6)
      M=np.zeros((dim,dim))
      DeltaX=np.zeros(dim)
   
      for i in range(dim):
         M[i][i]=mass[i//3] #square root of inverse masses
      #J=np.dot(L[0].T, L[1])  # for Lsorted
      J=np.linalg.pinv(L[0]).dot(L[1]) # for Lmassw
   
      #print "J\n", J
      DeltaX=np.array(x[1]-x[0]).flatten('F')  # need initial - final here.
      if logging[0] <1:
         logging[1].write('changes of Cartesian coordinates:\n'\
               +repr(DeltaX)+'\n')
    #  K=(DeltaX.dot(M)).dot(L[0])
    #  print K
      #K=M.dot(L[0]).T.dot(DeltaX)  # with Lsorted
      K=np.linalg.pinv(L[0]).dot(DeltaX)  # w p Lmassw
      #print K
      #K=L[0].T.dot(M).dot(M).dot(DeltaX)  # with Lmassw
   
      #K*=np.sqrt(np.pi)/2. #correction factor due to some magic reason  --> want to get rid of this!!!
   
      if logging[0]<2:
         # print the Duschinsky matrix in a nice format
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
       # for debugging: Check that q'=Jq+K as requested (q': normal modes of initial state)
       #  print "Zero:"
       #  print np.linalg.pinv(L[0]).dot(np.array(x[0]).flatten("F")) - J.dot(np.linalg.pinv(L[1]).dot(np.array(x[1]).flatten("F"))) +K
      return J, K 
   
   def calcspect(self, HR, freq, E, E0, N, M, T)
        #do nothing, it is a virtual function.
   
   def outspect(self, T, opt, linspect, E=0):
      """This function calculates the broadened spectrum given the line spectrum, 
      frequency-rage and output-file whose name is first argument. 
      As basis-function a Lorentzian is assumed with a common width.
      
      **PARAMETERS:**
      logging: This variable consists of two parts: logging[0] specifies the 
               level of print-out (which is between 0- very detailed
               and 4- only main information) and logging[1] is the file, already 
               opened, to write the information in.
      T:       temperature of the system
      opt:     a string that contains all options that were given for this part 
               in the input-file. See documentation 
               for details of it's allowed/used content
      linspect:The line-spectrum that has to be broadened: A array/matrix 
               with 3(!!) rows: 
               Frequency, intentensity and mode number (last one is 
               important for making multiple-particle spectra 
      E:       energy-shift of the 0-0 transition. Important if the excited 
               state is not the lowest and
               thermal equilibration with the lower states should be considered
   
      **RETURNS:**
      nothing; the key values (broadened spectra/ many-particle-app. 
               linespectra) are printed into log-files.
      
      """
      minint=0
      logging[1].write("\n STARTING TO CALCULATE BROADENED SPECTRUM.\n")
   
      omega, spectfile, gamma, gridpt, minfreq, maxfreq, shape, stick=handel_input(opt)
      #read file in format of linspect
      #sort spectrum with respect to size of elements
      index=np.argsort(linspect[1], kind='heapsort')
      linspect[0]=linspect[0][index] #frequency
      linspect[1]=linspect[1][index] #intensity
      linspect[2]=linspect[2][index] #mode
      #find transition with minimum intensity to be respected
   
      #truncate all transitions having less than 0.0001% of
      for i in range(len(linspect[1])):
         if linspect[1][i]>=1e-6*linspect[1][-1]:
            minint=i
            break
      if logging[0]<3:
         logging[1].write('neglect '+repr(minint)+' transitions, use only '+
                                repr(len(linspect[1])-minint)+" instead.\n")
   
         if logging[0]<2:
            logging[1].write('minimal and maximal intensities:\n'+
              repr(linspect[1][minint])+' '+repr(linspect[1][-1])+"\n")
      
      #important for later loops: avoiding '.'s speeds python-codes up!!
      logwrite=logging[1].write  
      #make nPA from OPA:
      if (re.search(r"to [\d]PA", opt, re.I) is not None) is True:
         n=re.findall(r"(?<=to )[\d](?=PA)", opt, re.I)
         if n[0]=='2':
            ind=linspect[2].argmin()
                           #  spectral frequency   0-0 transition   intensities
                           #      0-0 intensit.          modes        
            TPAf, TPAi=OPA2TPA(logwrite, linspect[0][minint:],linspect[0][ind] ,
                               linspect[1][minint:], linspect[1][ind], 
                               linspect[2][minint:], stick)
            index=np.argsort(TPAi,kind='heapsort')
            TPAi=TPAi[index] #resort by intensity
            TPAf=TPAf[index]
            minint=0
            # save time: look only on every second value. Could be reduced to log(n) here...
            for i in range(1,len(index),2):
               if TPAi[i]>=0.0001*TPAi[-1]:
                  minint=i
                  break
            TPAintens=TPAi[minint:]
            TPAfreq=TPAf[minint:]
            if logging[0]<3:
               logging[1].write('for TPA: again neglect '+repr(minint)+
                        ' transitions, use only '+repr(len(TPAi)-minint-1)+" instead.\n")
         elif n[0]=='3':
            ind=linspect[2].argmin()
            TPAfreq, TPAintens=OPA23PA(logwrite, linspect[0][minint:],
                                 linspect[0][ind] ,linspect[1][minint:], 
                                 linspect[1][ind], linspect[2][minint:], stick)
            minint=0
            index=np.argsort(TPAintens,kind='heapsort')
            TPAintens=TPAintens[index] #resort by intensity
            TPAfreq=TPAfreq[index]
            for i in range(len(TPAintens)):
               if TPAintens[i]>=0.0001*TPAintens[-1]:
                  minint=i
                  break
            TPAintens=TPAintens[minint:] #resort by intensity
            TPAfreq=TPAfreq[minint:]
            if logging[0]<3:
               logging[1].write('for 3PA: again neglect '+repr(minint)+
                        ' transitions, use only '+repr(len(TPAintens)-minint)+" instead.\n")
         else:
            if n[0]!='1':
               # there should be only 1,2 or 3 given!
               logging[1].write("to <n>PA was given but not recognised.\n")
            TPAfreq=linspect[0][minint:]
            TPAintens=linspect[1][minint:]
            if stick:
               logwrite=logging[1].write
               logwrite(u"\nSTICK-SPECTRUM IN ONE-PARTICLE APPROXIMATION "+
                                             "\n intensity   frequency\n")
               for k in range(len(TPAfreq)):
                  logwrite(u" %s  %s\n" %(TPAintens[k], TPAfreq[k]))
      else:
         n=re.findall(r"(?<=to nPA:)[ \d]*", opt, re.I)
         if n==[]:
            TPAfreq=linspect[0][minint:]
            TPAintens=linspect[1][minint:]
            if stick:
               logwrite=logging[1].write
               logwrite(u"\nSTICK-SPECTRUM IN ONE-PARTICLE APPROXIMATION \n"+
                                                 " intensity   frequency\n")
               for k in range(len(TPAfreq)):
                  logwrite(u" %s  %s\n" %(TPAintens[k], TPAfreq[k]))
         else:
            n=int(n[0])
            ind=linspect[2].argmin() #get position of 0-0 transition
            if re.search(r"parallel", opt, re.I) is not None:
               logging[1].write("\n REACHING OPA TO nPA-PART. DO IT IN PARALELL."+
                                       " MANY FILES WILL BE CREATED. \n")
               logging[1].write(" --------------------------------------------"+
                                              "--------------------------- \n")
               gotit=parallelOPA2nPA(logwrite, linspect[0][minint:],
                        linspect[0][ind], linspect[1][minint:], 
                        linspect[1][ind], linspect[2][minint:], n, stick, logging)
               if gotit:
                  logging[1].write("\n SUCCESSFULLY CALCULATED FULL-nPA SPECTRUM. \n")
                  return 0
               else:
                  logging[1].write("\n AN ERROR OCCURED. THE nPA SPECTRUM OR"+
                                            " BROADENING DID NOT SUCCEED. \n")
                  return 1
            else:
               logging[1].write("\n REACHING OPA TO NPA-PART. DOING IT NOT"+
                                                  " PARALELL. \n")
               logging[1].write(" ----------------------------------------"+
                                                       "-------- \n")
               TPAfreq, TPAintens=OPA2nPA(logwrite, linspect[0][minint:],
                        linspect[0][ind], linspect[1][minint:], 
                        linspect[1][ind], linspect[2][minint:], n, stick)
            index=np.argsort(TPAintens,kind='heapsort')
            TPAintens=TPAintens[index] #resort by intensity
            TPAfreq=TPAfreq[index]
            minint=0
            for i in range(len(TPAintens)):
               if TPAintens[i]>=1e-6*TPAintens[-1]:
                  minint=i
                  break
            if logging[0]<3:
               logging[1].write('for {0}PA: again neglect {1} \n'.format(n, minint)+
                        ' transitions, use only '+repr(len(TPAintens)-minint-1)+" instead.\n")
            TPAintens=TPAintens[minint:] #resort by intensity
            TPAfreq=TPAfreq[minint:]
   
   
      #find transition with minimum intensity to be respected
      #the range of frequency ( should be greater than the transition-frequencies)
      if omega==None:
         if minfreq==0:
            minfreq=np.min(TPAfreq)-20-gamma*15
         if maxfreq==0:
            maxfreq=np.max(TPAfreq)+20+gamma*15
      else:
         minfreq=omega[0]
         maxfreq=omega[-1]
      if logging[0]<3:
         logging[1].write('maximal and minimal frequencies:\n {0} {1}\n'.format(maxfreq, minfreq))
      #if no other grid is defined: use linspace in range
      if omega==None:
         omega=np.linspace(minfreq,maxfreq,gridpt)
         if logging[0]<2:
            logging[1].write("omega is equally spaced\n")
   
      sigma=gamma*2/2.355 #if gaussian used: same FWHM
      
      if gamma*1.1<=(maxfreq-minfreq)/gridpt:
         logging[1].write("\n !WARNING!\n THE GRID SPACING IS LARGE COMPARED TO THE WIDTH OF THE PEAKS.\n"
              "THIS CAN ALTER THE RATIO BETWEEN PEAKS IN THE BROADENED SPECTRUM!")
   
      index=np.argsort(TPAfreq,kind='heapsort') #sort by freq
      freq=TPAfreq[index]
      intens=TPAintens[index]
   
      mini=0
      if spectfile==None:
         out=logging[1]
      else:
         out = open(spectfile, "w")
   
      if spectfile==None: #that means spectrum is printed into log-file
         logwrite("broadened spectrum:\n frequency      intensity\n")
      outwrite=out.write
      #this shrinks the size of the spectral lines; hopefully accelerates the script.
      #intens, freq=concise(intens,freq, sigma)
      lenfreq=len(freq)
      maxi=lenfreq-1 #just in case Gamma is too big or frequency-range too low
      for i in range(0,lenfreq-1):
         if freq[i]>=10*gamma+freq[0]:
            maxi=i
            break
      if shape=='g':
         sigmasigma=2.*sigma*sigma # these two lines are to avoid multiple calculations of the same
         npexp=np.exp
         intens/=sigma # scale all transitions according to their width.
         for i in xrange(len(omega)): 
            omegai=omega[i]
            for j in range(maxi,lenfreq):
               if freq[j]>=10*gamma+omegai:
                  maxi=j
                  break
            for j in range(mini,maxi):
               if freq[j]>=omegai-10*gamma:
                  # else it becomes -1 and hence the spectrum is wrong
                  mini=max(j-1,0) 
                  break
            spect=0
            for k in range(mini,maxi+1):
               spect+=intens[k]*npexp(-(omegai-freq[k])*(omegai-freq[k])/(sigmasigma))
            outwrite(u" %f  %e\n" %(omegai, spect))
      else:  #shape=='l':
         gammagamma=gamma*gamma
         for i in xrange(len(omega)): 
            omegai=omega[i]
            for j in range(maxi,lenfreq):
               if freq[j]>=10*gamma+omegai:
                  maxi=j
                  break
            for j in range(mini,maxi):
               if freq[j]>=omegai-10*gamma:
                  # else it becomes -1 and hence the spectrum is wrong
                  mini=max(j-1,0) 
                  break
            omegai=omega[i]
            spect=0
            for k in range(mini,maxi+1):
               spect+=intens[k]*gamma/((omegai-freq[k])*(omegai-freq[k])+gammagamma)
            outwrite(u" %f   %e\n" %(omegai, spect))
      if spectfile!=None:
         #only close file if it was opened here
         out.close()

   def ReadG09(self, fileN):
      """ This function reads the required quantities from the gaussian-files

      **PARAMETERS**
      logging:     This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
                   and 4- only main information) and logging[1] is the file, already opened, to write the information in.
      fileN:       specifies the name of the g09-file

      **RETURNS**
      dim:         dimension of the problem (3*number of atoms)
      Coord:       Cartesian coordinates of the atoms in this state (as 1x, 1y, 1z, 2x,...)
      mass:        square root of masses of the atoms in atomi units
      F:           force constant matrix (Hessian of the PES); used to calculate the normal mode's motion and the frequencies
      E:           binding energy of the respective state/geometry in Hartree
      """
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
                           "modes: {1} \n Sqrt of masses in a.u. as read from log file\n{2}\n".format(dim/3,dim,mass))
      # Reading Cartesian coordinates
      Coord=np.zeros((3, dim/3))
      for j in range(len(tmp)):
         Coord[j%3][j/3]=tmp[j]

      Coord*=Angs2Bohr
      if logging[0]<1:
         logging[1].write("Cartesian coordinates (a.u.) \n{0}\n".format(Coord.T))

      # Reading of Cartesian force constant matrix  
      f=re.findall(r"Force constants in Cartesian coordinates: [\n\d .+-D]+", log, re.M)
      if f==[]:
         #if Freq was not specified in Gaussian-file:
         f=re.findall(r"The second derivative matrix:[XYZ\n\d .-]+", log, re.M)
         #try to find matrix from option "Force"
         assert f!=[], 'The input-file does not contain information on the force-constants!'
         #else: error message. The matrix is really needed...
         f_str=str([f[-1]])#
         lines=f_str.strip().split("\\n")
         print lines
         F=np.zeros((dim,dim))
         n=0
         k=0
         line=0
         assert len(lines)>2, "Can't find any forces. Force constants given in unknown format."
         for i in xrange(2,len(lines)):
            if i == dim+k-5*n+2: 
               #these are those lines where no forces are written to
               k=i-1
               n+=1
               line=n*5 
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
         for i in xrange(2,len(lines)):
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
            F[i][j]/= (mass[i//3]*mass[j//3]) 
      Etemp=re.findall(r'HF=-[\d.\n ]+', log, re.M)
      assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
      if re.search(r'\n ', Etemp[-1]) is not None:
         Etemp[-1]=Etemp[-1].replace("\n ", "") 
      if logging[0]<=1:
         logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
      E=-float(re.findall(r'[\d.]+', Etemp[-1])[0]) #energy is negative (bound state)
      return dim, Coord, mass, F, E

   def ReadG092(self, final):
      """ This function reads the required quantities from the gaussian-files

      **PARAMETERS**
      logging:     This variable consists of two parts: logging[0] specifies 
                   the level of print-out (which is between 0- very detailed
                   and 4- only main information) and logging[1] is the file, 
                   already opened, to write the information in.
      final:       specifies the name of the g09-file

      **RETURNS**
      grad:        gradient of the PES of excited state at ground state 
                   equilibrium-geometry
      E:           binding energy of the respective state/geometry in Hartree
      """
      files=open(final, "r")
      log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close
      Etemp=re.findall(r'HF=-[\d.\n ]+', log, re.M)
      assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
      if re.search(r'\n ', Etemp[-1]) is not None:
         Etemp[-1]=Etemp[-1].replace("\n ", "") 
      if logging[0]<=1:
         logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
      E=-float(re.findall(r'[\d.]+', Etemp[-1])[0]) #energy is negative (bound state)
     
      grad=re.findall(r"Number     Number              X              Y              Z\n [\-\.\d \n]+",log)
      Grad=re.findall(r"[\-\d\.]+[\d]{9}", grad[0])
      grad=np.zeros((len(Grad),1))
      for i in xrange(len(Grad)):
         grad[i]=float(Grad[i])
      return grad, E

   def ReadGO9_fchk(self, fileN):
      """ This function reads the required quantities from the gaussian-files

      **PARAMETERS**
      logging:     This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
                   and 4- only main information) and logging[1] is the file, already opened, to write the information in.
      fileN:       specifies the name of the g09-file

      **RETURNS**
      dim:         dimension of the problem (3*number of atoms)
      Coord:       Cartesian coordinates of the atoms in this state (as 1x, 1y, 1z, 2x,...)
      mass:        square root of masses of the atoms in atomi units
      F:           force constant matrix (Hessian of the PES); used to calculate the normal mode's motion and the frequencies
      E:           binding energy of the respective state/geometry in Hartree
      """
      #####
      #####
      ##### get normal modes: Vib-Modes 
      # Mapping the log file
      files=open(fileN, "r")
      log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close

      #get geometry of the molecule
      temp=[]
      temp=re.findall(r'Current cartesian coordinates  [RN\=\n \-\+.\d E]+', log)
      tmp=re.findall(r'[ -][\d.]+E[+-][\d]+', temp[0])

      # Determine atomic masses in a.u. Note mass contains sqrt of mass!!!
      atmwgt=re.findall(r"(?<=Real atomic weights     )[RN\=\-\+\d \.E\n]+",log)
      mtemp=[]
      foonum=0
      mtemp.append(re.findall(r'[\d.]+E[+-][\d]+',atmwgt[0]))
      dim=int(re.findall("(?<=N\= )[\d ]+", atmwgt[0])[0])
      dim*=3
      #if something is wrong with the dimensionality:
      assert len(tmp)==dim, 'Not all atoms were found! Something went wrong...{0}, {1}'.format(len(tmp),dim)

      mass=np.zeros(dim/3) # this is an integer since dim=3*N with N=atomicity
      for j in range(len(mtemp)):
         for k in range(len(mtemp[j])):
            mass[k+foonum]=np.sqrt(float(mtemp[j][k])*AMU2au) #elements in m are sqrt(m_i) where m_i is the i-th atoms mass
         foonum+=len(mtemp[j])
      assert not np.any(mass==0) , "some atomic masses are zero. Please check the input-file! {0}".format(mass)
      if logging[0]<2:
         logging[1].write("Number of atoms: {0}\nNumber of vibrational "
                           "modes: {1} \n Sqrt of masses in a.u. as read from log file\n{2}\n".format(dim/3,dim,mass))
      # Reading Cartesian coordinates
      Coord=np.zeros((3, dim/3))
      for j in range(len(tmp)):
         Coord[j%3][j/3]=float(tmp[j].replace('E','e')) # need to convert by some factor!!

      #Coord*=Angs2Bohr
      if logging[0]<1:
         logging[1].write("Cartesian coordinates (a.u.) \n{0}\n".format(Coord.T))

      # Reading of Cartesian force constant matrix  
      f=re.findall(r"(?<=Cartesian Force Constants   )[RN\=\d\+\-E \.\n]+", log, re.M)
      try:
         f_str=str([f[-1]])
      except IndexError:
         assert 1==2, "the intput-file %s has no information about Force-constant\nsorry. I can't work like this!" %(fileN)
      lines=f_str.strip().split("\\n")
      F=np.zeros((dim,dim))
      k=0
      m=0
      for i in xrange(1,len(lines)): # the first line is not part of the matrix
         elements=lines[i].replace('E','e').split()
         if m==len(F): # break this outer loop if inner is broken
            break
         for j in range(len(elements)):
            F[k][m]=F[m][k]=float(elements[j])
            if k==m:
               m+=1
               if m==len(F): # matrix is full; there will be only waste
                  break
               k=0
            else:
               k+=1
      if logging[0]<1:
         logging[1].write('F matrix as read from formcheck-file\n{0} \n'.format(F))
      for i in range(0,dim):
         for j in range(0,dim):
            F[i][j]/= (mass[i//3]*mass[j//3]) 

      Etemp=re.findall(r'(?<=Total Energy                               R  )[ \-\d.E\+]+', log, re.M)
      assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
      Etemp=float(Etemp[0].replace('E','e'))
      if logging[0]<=1:
         logging[1].write('temporary energy of state: {0}\n'.format(Etemp))
      return dim, Coord, mass, F, Etemp

   def ReadG09_fchk2(self, final):
      """ This function reads the required quantities from the gaussian-files

      **PARAMETERS**
      logging:     This variable consists of two parts: logging[0] specifies 
                   the level of print-out (which is between 0- very detailed
                   and 4- only main information) and logging[1] is the file, 
                   already opened, to write the information in.
      final:       specifies the name of the g09-file

      **RETURNS**
      grad:        gradient of the PES of excited state at ground state 
                   equilibrium-geometry
      E:           binding energy of the respective state/geometry in Hartree
      """
      files=open(final, "r")
      log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close
      #get energy
      E=re.findall(r'(?<=Total Energy                               R  )[ \-\d.E\+]+', log, re.M)
      assert len(E)>=1, 'Some error occured! The states energy can not be read.'
      E=float(E[0].replace('E','e'))
      if logging[0]<=1:
         logging[1].write('temporary energy of state: %d\n'%(E))
     
      #get gradient
      grad=re.findall(r"(?<=Cartesian Gradient)[R N=\d.\+\- E\n]+", log)
      Grad=re.findall(r"[\-\d\.\+E]+", grad[0])[1:]
      grad=np.zeros(len(Grad))
      for i in xrange(len(Grad)):
         grad[i]=float(Grad[i])
      return grad, E

   def getGaussianf(self, dim):
      files=open(final[0], "r") #open file and map it for better working
      mapedlog=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ) # i-th file containing freq calculations
      files.close
      #if file is of type 'Freq'
      if re.search(r" Freque ncies -- ", mapedlog, re.M) is not None:
         f1=re.findall(r" Frequencies -- [\d .-]+", mapedlog, re.M)# 
         f2=[re.findall(r"[- ]\d+.\d+", f1[j]) for j in range(len(f1))]
         s=0
         f=np.zeros((1,dim-6))
         for j in range(len(f2)):
            f[0][s:s+len(f2[j])]=f2[j]
            s+=len(f2[j])
      #if file is of type 'Force'
      elif re.search(r" Eigenvectors of the second derivative matrix:", mapedlog, re.M) is not None:
         f2=re.findall(r"(?<=Eigenvectors of the second derivative matrix:)[\d\n\-\. XYZ Eigenvalues]+", mapedlog, re.M)
         f1=re.findall(r" Eigenvalues -- [\d .]+", f2[0], re.M)# 
         f2=[re.findall(r"\d.\d+", f1[j]) for j in range(len(f1))]
         s=0
         f=np.zeros((1,dim-6))
         for j in range(0,len(f2)):
            for i in range(len(f2[j])):
               f[0][s+i]=np.sqrt(float(f2[j][i]))*2590.839/Hartree2cm_1*1.15
            s+=len(f2[j])
      #print "gaussian freq"
      #for i in range(len(f[0])):
      #   print f[0][i]*Hartree2cm_1 
      return f

   def getGaussianL(self, dim):
      #first, check in which format they are present
      files=open(final[0], "r") #open file and map it for better working
      mapedlog=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ) # i-th file containing freq calculations
      files.close
      b=0
      L=np.zeros((dim, dim-6))
      #if file is of type 'Freq'
      if re.search(r" Eigenvectors of the second derivative matrix:", mapedlog, re.M) is not None:
         f1=re.findall(r"(?<=Eigenvalues --)  [\d .\n XYZ\-]+", mapedlog, re.M)
         # to remove the lines with 'Eigenvalues':
         for k in range(len(f1)):
            f2=re.findall(r"[- ]\d\.[\d]+", f1[k]) 
            s=len(f2)//dim 
            #s should be 5 but not in last line
            for j in range(dim):
               for i in range(s):
                  L[j][b+i]=f2[i+s*(j+1)]
            b+=s
         #renormalise L
         for j in range(dim): 
            norm=L[j].dot(L[j].T)
            if norm>1e-12:
               L[j]/=np.sqrt(norm)
      else:
         print "there is no other method implemented until now!"
      
      return L

   def replace(self, files, freq, L): # this is G09-method
      """ This function creates a new file (determined by files, ending with 
      ".rep" and copies the log-file (files) into it, replacing the frequencies and 
      normal modes by those calculated by smallscript.
      The function is suited to test, whether these results coincide qualitatively with the gaussian's.

      ** PARAMETERS: **
      logging:This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
              and 4- only main information) and logging[1] is the file, already opened, to write the information in.
      log:    i
      files:  file taken as basis ('.rep' added to be used for replacements)
      freq:   frequencies to be inserted
      L:      normal modes to be inserted  

      no return-statements (results are written to file)

      **NOTE:**
      The code is based on a function from stevaha (http://stackoverflow.com/questions/1597649/replace-strings-in-files-by-python)
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

   def ReadGAMESS(self, fileN):
      """ This function reads the required quantities from the GAMESS-log-files

      **PARAMETERS**
      logging:     This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
                   and 4- only main information) and logging[1] is the file, already opened, to write the information in.
      fileN:       specifies the name of the g09-file

      **RETURNS**
      dim:         dimension of the problem (3*number of atoms)
      Coord:       Cartesian coordinates of the atoms in this state (as 1x, 1y, 1z, 2x,...)
      mass:        square root of masses of the atoms in atomi units
      F:           force constant matrix (Hessian of the PES); used to calculate the normal mode's motion and the frequencies
      E:           binding energy of the respective state/geometry in Hartree
      """
      # Mapping the log file
      files=open(fileN, "r")
      log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close

      # Determine atomic masses in a.u. Note mass contains sqrt of mass!!!
      #atmwgt=re.findall(r"ATOMIC WEIGHTS \(AMU\)", log)
      atmwgt=re.findall(r"(?<=ATOMIC WEIGHTS \(AMU\))[ \.\w\n]+", log)
      mtemp=[]
      foonum=0
      mtemp.append(re.findall(r'[\d.]{7,10}(?=\n)',atmwgt[-1]))
      dim=3*len(mtemp[0])

      mass=np.zeros(dim/3) # this is an integer since dim=3*N with N=atomicity
      for j in range(len(mtemp[0])):
         mass[j]=np.sqrt(float(mtemp[0][j])*AMU2au) #elements in m are sqrt(m_i) where m_i is the i-th atoms mass
      assert not np.any(mass==0) , "some atomic masses are zero. Please check the input-file! {0}".format(mass)
      if logging[0]<2:
         logging[1].write("Number of atoms: {0}\nNumber of vibrational "
                           "modes: {1} \n Sqrt of masses in a.u. as read from log file\n{2}\n".format(dim/3,dim,mass))
      
      # Reading Cartesian coordinates
      temp=re.findall(r'(?<=ATOM             X              Y              Z\n\n)[\n -.\d\w]+', log)  
                 #+ need even more in structure
      tmp=re.findall(r'[ -][\d]+.[\d]+', temp[-1])
      Coord=np.zeros((3,dim//3))
      for i in range(3):
         for j in range(dim//3):
            Coord[i][j]=float(tmp[3*i+j])
      Coord*=Angs2Bohr
      if logging[0]<1:
         logging[1].write("Cartesian coordinates (a.u.) \n{0}\n".format(Coord.T))
      # Reading of Cartesian force constant matrix  
      f=re.findall(r"(?<=CARTESIAN FORCE CONSTANT MATRIX\n)[\d\w \.\-\n]+", log, re.M)
      f_str=str([f[-1]])
      lines=f_str.strip().split("\\n")[5:]
      #lines include each three six elements.
      n=0
      k=0
      s=0
      F=np.zeros((dim,dim))
      for i in xrange(len(lines)):
         f=re.findall(r"[ -]\d\.[\d]+", lines[i+s], re.M) #this is one number each
         if len(f)!=6:
            n+=6
            k=n
            s+=3
            if k==dim:
               #this will finally finish this loop (since the last lines are no more in the
               #force constant matrix; I just don't know, how to make it easier
               break
            continue
         for j in range(6):
            #some numbers will be counted twice but this shouldn't matter
            F[k][n+j]=F[n+j][k]=f[j]
         k+=1

      if logging[0]<1:
         logging[1].write('F matrix as read from log file\n{0} \n'.format(F))
      for i in range(0,dim):
         for j in range(0,dim):
            F[i][j]/= (mass[i//3]*mass[j//3]) 

      #get energy of the state
      Etemp=re.findall(r'(?<=TOTAL ENERGY =)[\-\d. ]+', log, re.M)
      assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
      if logging[0]<=1:
         logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
      E=float(Etemp[-1]) #energy is negative (bound state)
      return dim, Coord, mass, F, E

   def ReadGAMESS2(self, final):
      """ This function reads the required quantities from the GAMESS-log-files

      **PARAMETERS**
      logging:     This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
                   and 4- only main information) and logging[1] is the file, already opened, to write the information in.
      final:       specifies the name of the g09-file

      **RETURNS**
      grad:        gradient of the PES of excited state at ground state equilibrium-geometry
      E:           binding energy of the respective state/geometry in Hartree
      """
      files=open(final, "r")
      log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close
      Etemp=re.findall(r'(?<=TOTAL ENERGY =)[\-\d. ]+', log, re.M)
      assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
      if logging[0]<=1:
         logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
      E=float(Etemp[-1]) #energy is negative (bound state)
      #grad=re.findall(r"Number     Number              X              Y              Z\n [\-\.\d \n]+",log)
      #Grad=re.findall(r"[\-\d\.]+[\d]{9}", grad[0])
      #grad=np.zeros((len(Grad),1))
      #for i in xrange(len(Grad)):
         #grad[i]=float(Grad[i])
      return grad, E

   def getGamessLf(self, dim):
      #open file and map it for better working
      files=open(final[0], "r") 
      # i-th file containing freq calculations
      mapedlog=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ) 
      files.close
      b=0
      L=np.zeros((dim, dim-6))
      # this is end of heading for the frequencies. Don't wonder about this text.
      if re.search(r"     REDUCED MASSES IN AMU.\n", mapedlog, re.M) is not None:
         f1=re.findall(r"(?<=  IR INTENSITY:   )  [\d .\-\w]+", mapedlog, re.M)
         #now the whole matrix L is within f1  (each until ':' of next 'frequency'
         for k in range(len(f1)):
            f2=re.findall(r"[- ]\d\.[\d]+", f1[k]) 
            s=len(f2)//dim 
            #s should be 5 but not in last line
            for j in range(dim):
               for i in range(s):
                  L[j][b+i]=f2[i+s*(j+1)]
            b+=s
         #renormalise L
         for j in range(dim): 
            norm=L[j].dot(L[j].T)
            if norm>1e-12:
               L[j]/=np.sqrt(norm)
      else:
         print "there is no other method implemented until now!"
      f1=re.findall(r"(?<=FREQUENCY:    )  [\d .]+", mapedlog, re.M)
      f2=[re.findall(r"\d+.\d+", f1[j]) for j in range(len(f1))]
      s=0
      f=np.zeros((1,dim))
      for j in range(len(f2)):
         f[0][s:s+len(f2[j])]=f2[j]
         s+=len(f2[j])
      return f, L

   def ReadNWChem(self, fileN):
      """ This function reads the required quantities from the NWChem-files

      **PARAMETERS**
      logging:     This variable consists of two parts: logging[0] specifies the level of 
                   print-out (which is between 0- very detailed
                   and 4- only main information) and logging[1] is the file, already opened, to 
                   write the information in.
      fileN:       specifies the name of the g09-file

      **RETURNS**
      dim:         dimension of the problem (3*number of atoms)
      Coord:       Cartesian coordinates of the atoms in this state (as 1x, 1y, 1z, 2x,...)
      mass:        square root of masses of the atoms in atomic units
      F:           force constant matrix (Hessian of the PES); used to calculate the normal mode's 
                   motion and the frequencies
      E:           binding energy of the respective state/geometry in Hartree
      """
      # Mapping the log file
      files=open(fileN, "r")
      log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close

      # Determine atomic masses in a.u. Note mass contains sqrt of mass!
      WeightCoord=re.findall(r"(?<=atom    #        X              Y              Z            mass\n)[\d\w\-\+ \.\n]+", log)
      mtemp=re.findall(r"[ -]\d\.[\dD\.\+\-]+(?=\n)", WeightCoord[-1]) #only weights remain
      dim=3*len(mtemp)

      mass=np.zeros(dim/3) # this is an integer since dim=3*N with N=atomicity
      for j in range(len(mtemp)):
         mass[j]=np.sqrt(float(mtemp[j].replace('D','e'))*AMU2au) #mass[j] is square-root of masses
      assert not np.any(mass==0) , "some atomic masses are 0. Please check the input-file! {0}".format(mass)
      if logging[0]<2:
         logging[1].write("Number of atoms: {0}\nNumber of vibrational "
                          "modes: {1} \n Sqrt of masses in a.u. as read "
                          "from log file\n{2}\n".format(dim/3,dim,mass))
      
      # Reading Cartesian coordinates
      tmp=re.findall(r'[ -]\d\.[\dD\.\+\-]+', WeightCoord[-1])
      tmp
      Coord=np.zeros((3,dim//3))
      k=0
      for j in range(dim//3):
         for i in range(3):
            Coord[i][j]=float(tmp[i+3*j+k].replace('D','e'))
         k+=1 # important to skip masses in tmp 
      Coord*=Angs2Bohr
      if logging[0]<1:
         logging[1].write("Cartesian coordinates (a.u.) \n{0}\n".format(Coord.T))

      # Reading of Cartesian force constant matrix  
      f=re.findall(r"(?<=MASS-WEIGHTED PROJECTED HESSIAN \(Hartree/Bohr/Bohr/Kamu\)\n)[\dD \.\-\+\n]+", log, re.M)
      f_str=str([f[-1]])
      lines=f_str.strip().split("\\n")
      F=np.zeros((dim,dim))
      n=0
      k=0
      for i in xrange(5,len(lines)):
         if i == dim+k-10*n+8: #this is due to the lines after the first row
            k=i-4
            n+=1
            continue
         if dim+k-10*n+5<=i<=dim+k-10*n+7:
            continue
         elements=lines[i].replace('D','e').split()
         for j in range(1,len(elements)):
            # convert to a.u.  (masses)
            F[int(elements[0])-1][j-1+10*n]=float(elements[j])/(1000*AMU2au)
            F[j-1+10*n][int(elements[0])-1]=float(elements[j])/(1000*AMU2au)
      if logging[0]<1:
         logging[1].write('F matrix as read from log file\n{0} \n'.format(F))

      #get energy of the state
      Etemp=re.findall(r'(?<=Total DFT energy =)[\-\d. ]+', log, re.M)
      assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
      if logging[0]<=1:
         logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
      E=float(Etemp[-1]) #energy is negative (bound state)
      return dim, Coord, mass, F, E

   def ReadNWChem2(self, final):
      """ This function reads the required quantities from the NWChem-files

      **PARAMETERS**
      logging:     This variable consists of two parts: logging[0] specifies the level of 
                   print-out (which is between 0- very detailed and 4- only main information) 
                   and logging[1] is the file, already opened, to write the information in.
      final:       specifies the name of the g09-file

      **RETURNS**
      grad:        gradient of the PES of excited state at ground state equilibrium-geometry
      E:           binding energy of the respective state/geometry in Hartree
      """
      files=open(final, "r")
      log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close
      Etemp=re.findall(r'(?<=Total energy =)[\-\d. ]+', log, re.M)
      assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
      if logging[0]<=1:
         logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
      E=float(Etemp[-1]) #energy is negative (bound state)
      #grad=re.findall(r"Number     Number              X              Y              Z\n [\-\.\d \n]+",log)
      #Grad=re.findall(r"[\-\d\.]+[\d]{9}", grad[0])
      #grad=np.zeros((len(Grad),1))
      #for i in xrange(len(Grad)):
         #grad[i]=float(Grad[i])
      return grad, E

   def getNwchemLf(self,final, dim):
      #open file and map it for better working
      files=open(final[0], "r") 
      # i-th file containing freq calculations
      mapedlog=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ) 
      files.close
      b=0
      L=np.zeros((dim, dim-6))
      # get L-matrix
      f1=re.findall(r"(?<=P.Frequency  )  [\d .\-\n]+", mapedlog, re.M)
      #now the whole matrix L is within f1  (each until ':' of next 'frequency'
      for k in range(1,len(f1)): #start from 1 to skip transl./rotation
         f2=re.findall(r"[- ]\d\.[\d]{5}", f1[k])  # in this scheme I throw frequencyes out
         s=len(f2)//dim 
         #s should be 6 but not in last line
         for j in range(dim):
            for i in range(s):
               L[j][b+i]=f2[i+s*(j)]
         b+=s
      #renormalise L
      for j in range(dim): 
         norm=L[j].dot(L[j].T)
         if norm>1e-12:
            L[j]/=np.sqrt(norm)
      #get frequency
      f1=re.findall(r"(?<= P.Frequency   )  [\d .\-]+", mapedlog, re.M)
      f2=[re.findall(r"\d+.\d+", f1[j]) for j in range(1,len(f1))] #start from 1 to skip first 6 modes
      s=0
      f=np.zeros((1,dim-6))
      for j in range(len(f2)):
         f[0][s:s+len(f2[j])]=f2[j]
         s+=len(f2[j])
      return f, L
   ## END OF FUNCTION DEFINITIONS

version=0.0.1
# End of Spect.py
