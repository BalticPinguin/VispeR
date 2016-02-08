#!/usr/bin/python2
# filename: Spect.py
import numpy as np
from scipy.linalg.lapack import dsyev
import re, mmap, os.path, math, sys
import Read
import MultiPart

Hartree2GHz=6.579684e6                                     
Hartree2cm_1=219474.63 

# ============ CHANGELOG =================

class  Spect:
   """ This is the general class that all spectra belong to. It contains the fundamental 
   quantities that are required and most of the functions to use. 
   Some of the functions are redefined in the respective sub-classes.  """
   # Below are the conversion factors and fundamental constant
   AMU2au=1822.88839                                          
   Angs2Bohr=1/0.52917721092                                  
   Hartree2GHz=6.579684e6                                     
   Hartree2cm_1=219474.63 
   
   # BEGIN OF DATA-DEF.
   Energy=[]
   dim=0
   CartCoord=[]
   read=[]  # this object is valid only after __init__. Maybe it doesn't belong here...
   spect=[[0],[0],[0]]
   
   # The following additional elements are members of all 
   #
   #  logging
   #  f  (frequencies in a 2-D array)
   #  Lmassw
   #  T  (temperature of the system)
   #  broadopt
   #  states1, states2

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
         log:         opened file to write in
         """
         if logfile==[]:
            log=open("calculation.log", "a")
         else:
            s=logfile[-1].strip()
            log=open(s, "a")

         #remove all white spaces
         mode= mode.strip() 
         #search, which option was set.
         if mode in ['all', "ALL", "All", "0", 0]:
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
            log.write("logging-mode "+mode+" not recognized. Using 'important' instead\n")
         return logging, log
      

      #START LOG-FILE DEFINITION     
      #invoke logging (write results into file specified by 'out: ' or into 'calculation.log')
      logfile=re.findall(r"(?<=out: )[\w.,\(\) \=;:\-_]+", f, re.I)
      try:
         loglevel=re.findall(r"(?<=print:)[\w \d]+",f, re.I)[-1]
         self.logging = invokeLogging(logfile, loglevel )
      except IndexError:
         self.logging = invokeLogging(logfile)
      
      # now,  write the header to the output-file.
      self.logging[1].write("\n==================================================================\n"
               "===================  output of smallscript  ======================\n"
               "==================================================================\n\n")
      self.logging[1].write("   INPUT-FILE:\n")
      self.logging[1].write(f)
      self.logging[1].write(" \n   END OF INPUT-FILE \n\n")
      #self.logging[1].write("calculations to be done: %s\n"%(todo))
      #END LOG-FILE DEFINITION     
      
      self.Energy=np.zeros(2)
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
      opt=re.findall(r"(?<=opt:)[\d\s\w.,\(\) \=;:-]+", f, re.M)
      if opt!=[]:
         self.opt=opt[-1]
      else:
         self.opt=" " 
         #make some empty string that the search-commands don't fail in ReadData.
      broadopt=re.findall(r"(?<=broaden:)[\d\s\w.,\(\) \=;:-]+", f, re.M)
      if broadopt!=[]:
         self.broadopt=broadopt[-1]
      else:
         self.broadopt=" " 
      self.ReadData(initial, final)
      self.T=re.findall(r"(?<=T=)[\d .]+", self.opt, re.M)
      if self.T==[]:
         self.T=300
      else:
         self.T=float(self.T[-1])
      self.T*=8.6173324e-5/27.21138386 # multiplied by k_B in hartree/K

      # get number of maximal excited/ground from opt:
      states=re.findall(r"(?<=states=)[\d ,]*", self.opt, re.I)
      if len(states)==0:
         # i.e. no states given
	 self.states1=5
         self.states2=0
      else:
	 try:
            #i.e. 1 number given (else, the 'int' will give an error)
	    self.states1=int(states[0])
	    self.states2=0
	    logging[1].write("number of states: %d and %d \n"%(self.states1, self.states2))
	 except ValueError:
            try:
               # two numbers given, try to extract them
               self.states1=int(states[0].split(",")[0])
               self.states2=int(states[0].split(",")[1])
               logging[1].write("number of states: %d and %d \n"%(states1, states2))
            except ValueError:
               #unknown format. Use default and give a respective message.
               self.states1=5
               self.states2=0
               logging[1].write("!!number of vibrational states {0} is not an integer.",
                                    " Use default instead.\n".format(self.states1, self.states2))
      #END READING DATA FROM FILE

   def GetL(self, mass, F):
      """ Function that calculates the frequencies and normal modes from force 
      constant matrix with and without projection onto internal degrees of freedom
   
      **argumets**
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
            self.logging[1].write('imaginary frequencies occured. The absolute'
                    ' values are used in the following.\n{0}\n'.format(f[i]))
            f[i]=np.abs(f[i])
         M=np.eye(self.dim)
         for j in range(self.dim):
            M[j,j]/=mass[j//3]
         Lmass[i]=M.dot(L[i])
         if self.logging[0]<2:
            self.logging[1].write("Frequencies (cm-1)\n"+\
                  repr(f[i]*Hartree2cm_1)+"\nL-matrix \n"+ repr(Lmass[i])+"\n")
      #np.set_printoptions(threshold=np.nan,linewidth=500, suppress=True)
      f[1],L[1]=SortL(np.linalg.pinv(L[0]).dot(L[1]),L[1],f[1])
      Lmass[1]=M.dot(L[1]) # recalculate Lmass!
      return f, L, Lmass
   
   def Duschinsky(self, mass):
      """
      This function calculates the shift between two electronic states 
      (whose geometry is known, see x) as well as the
      Duschinsky-rotation matrix.
      **PARAMETERS:**
      mass:    array of square-roots of nuclear masses (length: N)
      **RETURN:**
      J:    Duschinsky-rotation matrix
      K:    displacement-vector of energy-minima in normal coordinates
      """
      self.J=np.zeros((self.dim-6,self.dim-6))
      self.K=np.zeros(self.dim-6)
      M=np.zeros((self.dim,self.dim))
      DeltaX=np.zeros(self.dim)
   
      for i in range(self.dim):
         M[i][i]=mass[i//3] #square root of inverse masses
      #J=np.dot(L[0].T, L[1])  # for Lsorted
      self.J=np.linalg.pinv(self.Lmassw[0]).dot(self.Lmassw[1]) # for Lmassw
   
      #print "J\n", J
      DeltaX=np.array(self.CartCoord[1]-self.CartCoord[0]).flatten('F')  # need initial - final here.
      if self.logging[0] <1:
         self.logging[1].write('changes of Cartesian coordinates:\n'\
               +repr(DeltaX)+'\n')
    #  K=(DeltaX.dot(M)).dot(L[0])
    #  print K
      #K=M.dot(L[0]).T.dot(DeltaX)  # with Lsorted
      self.K=np.linalg.pinv(self.Lmassw[0]).dot(DeltaX)  # w p Lmassw
      #print K
      #K=L[0].T.dot(M).dot(M).dot(DeltaX)  # with Lmassw
   
      #K*=np.sqrt(np.pi)/2. #correction factor due to some magic reason  --> want to get rid of this!!!
   
      if self.logging[0]<2:
         # print the Duschinsky matrix in a nice format
         self.logging[1].write('Duschinsky rotation matrix:\n')
         k=range(0,self.dim-6)
         s=0
         t=min(s+5,self.dim-6)
         while s<self.dim-6:
            for temp in range(s,t):
               self.logging[1].write("               %d "%(k[temp]+1))
            self.logging[1].write("\n")
            for j in range(len(self.J)):
               self.logging[1].write(" %03d"%(j+1))
               for temp in range(s,t):
                  self.logging[1].write("   %+.5e"%(self.J[j][k[temp]]))
               self.logging[1].write("\n")
            s=t
            t=min(s+5,self.dim-6)
         self.logging[1].write('\nDuschinsky displacement vector:\n')
   
         for j in range(len(self.K)):
            self.logging[1].write("  %d    %e\n"%(j+1, self.K[j]))
      #return J, K   # now, everything is part of 'self'.
   
   def outspect(self, E=0):
      """This function calculates the broadened spectrum given the line spectrum, 
      frequency-rage and output-file whose name is first argument. 
      As basis-function a Lorentzian is assumed with a common width.
      
      **PARAMETERS:**
      E:       energy-shift of the 0-0 transition. Important if the excited 
               state is not the lowest and
               thermal equilibration with the lower states should be considered
   
      **RETURNS:**
      nothing; the key values (broadened spectra/ many-particle-app. 
               linespectra) are printed into log-files.
      
      """
      def handel_input(opt):
         #set default values (to have all variables set)
         gridfile=None
         gamma=1 #by default: only slight broadening
         gridpt=5000
         omega=None
         minfreq=0
         maxfreq=0
         shape='g'
         stick=False
      
         tmpgrid=re.findall(r"(?<=grid=)[ \=\s\w\.;]+", opt, re.M)
         if len(tmpgrid)==1: 
         # i.e. if grid is specified
            grid=re.findall(r"[\w\.]+", tmpgrid[0], re.M)
            if len(grid)==1:
               #that means, if either one number (# of gridpoints or a file) is given
               try:
                  gridpt=float(grid[0])
               except ValueError: # if grid is no a number
                  gridfile=grid[0]
            elif len(grid)==3:
               # that means there is the number of gridpoints, min- and max frequency given
               gridpt=float(grid[0])
               minfreq=float(grid[1])
               maxfreq=float(grid[2])
            if gridfile!=None:
               #read file in format of spect
               grid=[]
               with open(gridfile) as f:
                  lis=[line.split() for line in f]  # create a list of lists
                  for i,x in enumerate(lis):        # get the list items 
                     grid.append(float(x[0]))
               omega=np.zeros(len(grid))
               for i in range(len(grid)):
                  omega[i]=grid[i]
         #see, whether a broadening is given
         if (re.search(r"(?<=gamma=)[ \d\.,]+", opt, re.I) is not None)  is True:
            gamma=re.findall(r"(?<=gamma=)[ \d\.]+", opt, re.I)
            gamma=float(gamma[0])
      
         shape=re.findall(r"(?<=shape=)[ \w]+", opt, re.I)
         if len(shape)>0:
         # there are several options each
            if shape[0] in ["lorentzian", "Lorentzian", "L", "l"]:
               shape="l"
            elif shape[0] in ["gaussian", "Gaussian", "G", "g"]:
               shape="g"
   
         if (re.search(r"stick", opt, re.I) is not None) is True:
            stick=True
   
         spectfile=re.findall(r"(?<=spectfile=)[\w.]+", opt, re.I)
         if spectfile==[]:
            spectfile=re.findall(r"(?<=spectfile= )[\w.]+", opt, re.I)
            if spectfile==[]:
               spectfile=None
            else:
               spectfile=spectfile[-1]
         else:
            spectfile=spectfile[-1]
         return omega, spectfile, gamma, gridpt, minfreq, maxfreq, shape, stick

      minint=0
      self.logging[1].write("\n STARTING TO CALCULATE BROADENED SPECTRUM.\n")
      
      omega, spectfile, gamma, gridpt, minfreq, maxfreq, shape, stick=handel_input(self.broadopt)
      #read file in format of spect
      #sort spectrum with respect to size of elements
      index=np.argsort(self.spect[1], kind='heapsort')
      self.spect[0]=self.spect[0][index] #frequency
      self.spect[1]=self.spect[1][index] #intensity
      self.spect[2]=self.spect[2][index] #mode
      #find transition with minimum intensity to be respected
   
      #truncate all transitions having less than 0.0001% of
      for i in range(len(self.spect[1])):
         if self.spect[1][i]>=1e-6*self.spect[1][-1]:
            minint=i
            break
      if self.logging[0]<3:
         self.logging[1].write('neglect '+repr(minint)+' transitions, use only '+
                                repr(len(self.spect[1])-minint)+" instead.\n")
   
         if self.logging[0]<2:
            self.logging[1].write('minimal and maximal intensities:\n'+
              repr(self.spect[1][minint])+' '+repr(self.spect[1][-1])+"\n")
      
      #important for later loops: avoiding '.'s speeds python-codes up!!
      logwrite=self.logging[1].write  
     
      #make nPA from OPA if requested.
      n=re.findall(r"(?<=to nPA:)[ \d]*", self.broadopt, re.I)
      if n!=[]:
         MakeFull=MultiPart.OPAtoNPA(float(n[-1].strip()))
         self.logging[1].write("\n REACHING OPA TO NPA-PART. \n")
         self.logging[1].write(" ----------------------------------------"+
                                                 "-------- \n")
         MakeFull.GetSpect(self.spect, minint)
         TPAintens, TPAfreq=MakeFull.Calc()
      else: 
         TPAfreq=linspect[0][minint:]
         TPAintens=linspect[1][minint:]
   
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
      if self.logging[0]<3:
         self.logging[1].write('maximal and minimal frequencies:\n {0} {1}\n'.format(maxfreq, minfreq))
      #if no other grid is defined: use linspace in range
      if omega==None:
         omega=np.linspace(minfreq,maxfreq,gridpt)
         if self.logging[0]<2:
            self.logging[1].write("omega is equally spaced\n")
   
      sigma=gamma*2/2.355 #if gaussian used: same FWHM
      
      if gamma*1.1<=(maxfreq-minfreq)/gridpt:
         self.logging[1].write("\n !WARNING!\n THE GRID SPACING IS LARGE COMPARED TO THE WIDTH OF THE PEAKS.\n"
              "THIS CAN ALTER THE RATIO BETWEEN PEAKS IN THE BROADENED SPECTRUM!")
   
      index=np.argsort(TPAfreq,kind='heapsort') #sort by freq
      freq=TPAfreq[index]
      intens=TPAintens[index]
   
      mini=0
      if spectfile==None:
         out=self.logging[1]
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

   def ReadData(self, initial, final):
      """ This function gathers most essential parts for calculation of 
      HR-factors from g09-files. That is: read neccecary information from the  
      g09-files and calculate HR-factors as well as the  Duschinsky-rotation 
      matrix and the shift between minima (needed later if the option Duschinsky 
      is specified)
      
      """                                                                                             
      # create the member read. It is an instance of a class, made for reading data.
      self.read = Read.Read(final, initial)
         
      #read coordinates, force constant, binding energies from log-files and 
      # from the file, using the type of file that is known now...
      self.dim, Coord, mass, A, self.Energy[0]=self.read.ReadAll('i')
      F, P=self.quantity(self.dim) #creates respective quantities (empty)
      if self.logging[0]==0:
         self.logging[1].write("Dimensions: "+ str(self.dim)+ '\n Masses: '+ str(mass**2)+"\n")
      F[0]=A
      self.CartCoord[0]=Coord
      self.dim, self.CartCoord[1], mass, F[1], self.Energy[1]=self.read.ReadAll('f') 

      if self.logging[0]<3:
         self.logging[1].write('difference of minimum energy between states:'
                        ' Delta E= {0}\n'.format((self.Energy[0]-self.Energy[1])*Hartree2cm_1))
         if self.logging[0]<2:
            self.logging[1].write('Cartesion coordinates of initial state: \n{0}\n'.format( self.CartCoord[0].T/self.Angs2Bohr))
            self.logging[1].write('Cartesion coordinates of final state: \n{0}\n Forces:\n'.format( self.CartCoord[1].T/self.Angs2Bohr))
            self.logging[1].write('initial state: \n{0}\n'.format(F[0]))
            self.logging[1].write('final state: \n {0}\n'.format(F[1]))
   
   
      #Calculate Frequencies and normal modes
      self.f, Lsorted, self.Lmassw=self.GetL(mass, F)
      self.Duschinsky(mass)
   
   def quantity(self, dim):
      """ Here some frequencies are defined; it is just for clearer code.
      This function is called by CalculationHR.
   
      **PARAMETERS**
      dim      dimension of matrices/vectors
      """
      F=np.zeros((2, dim, dim)) 
      self.CartCoord=np.zeros((2, 3, dim//3))
      P=np.zeros((2, dim,dim))
      return F, P
   ## END OF FUNCTION DEFINITIONS

#version=0.0.1
# End of Spect.py
