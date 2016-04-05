#!/usr/bin/python2
#filename: FC_Spects.py

#including the class everything stems from:
#include [[Spect.py]]
import Spect
import re, mmap, math
import numpy as np

#  CHANGELOG 
# ===========
# to version 0.1.5:  
#  1) removed GFC-class  
#  2) Added term for energy-correction in case of gradient.
#  3) Changed setting self.type.
#  4) repaired CFC_spect to work with the new structure.
#
# to version 0.1.0:  
#  1) Added some sense to GFC  
#  2) Removed the HR-classes from the code  
#  3) In FC: final states frequencies in use.

class FC_spect(Spect.Spect): # import class Spect from file Spect.
   """ All models that are based on FC-scheme are based on this class. 
      Since it doesn't need any more data than the DR-scheme, everything 
      initialised for Spect is initialized here already and hence no extra
      init()-function is needed.
   """
   def __init__(self, f): 
      """Initialisation of the function.
         In addition to the main class Spect, here only the Huang-Rhys factors
         are computed and a respective threshold for smallest HR-factors to be 
         accounted for is set.
      """
      # first, initialise as done for all spectra:
      if self.type=="Spect":
         #don't set it if this is CFC-calculation.
         self.type='FC'
      Spect.Spect.__init__(self, f)
      # now, calculate HR-spect in addition.
      HRthresh=re.findall(r"(?<=HRthreshold=)[ \d.]+",self.opt,re.M)
      if HRthresh==[]:
         HRthresh=0.015
      else: 
         HRthresh=float(HRthresh[-1])
      self.HuangR(HRthresh)

   def calcspect(self):
      """The main function of this class, calculating the spectrum in one-particle 
         approximation.
      """
      # M,N are maximal numbers of vibrational modes 
      #     (+1, since they should be arrived really; count from 0)
      self.states1+=1
      self.states2+=1
      #as frequency, use the final states frequency.
      f=self.f[1]

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
          faktNM=float(fact(M)*fact(N))
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
    
      n=len(self.HR) #=len(freq)
      setM=False
      if self.states2==1: 
         # there was not specified, how many vibr. states in ground-state 
         #     should be taken into account
         setM=True
         self.states2=max(3,int(-1.1*self.HR[0]*self.HR[0]+6.4*self.HR[0]+9.))
         #i.e. take M following an empirical value as function of the HR-file
      assert n>0, "There is no Huang-Rhys factor larger than the respective"+\
                                        "threshold. No mode to be calculated."
      #if setM: the size of these arrays will be overestimated.
      FC=np.zeros((n,self.states2*self.states1-1))
      uency=np.zeros((n,self.states2*self.states1-1)) #frequency

      #avoiding dots accelerates python quite a lot
      loggingwrite=self.log.write
      npexp=np.exp 
      E=self.Energy[0]-self.Energy[1]
      #set  0->0 transition:
      sgnE=np.sign(E)
      FC00=1.00
      if sgnE==0:
         sgnE=1
      #here a goes over all modes
      for a in xrange(n):
         if setM:
            # set M to fit best to the value at that moment.
            #M=max(3,int(-1.1*self.HR[a]*self.HR[a]+6.4*self.HR[a]+9.))
            M=self.states2
         for j in range(self.states1):  # initial state
            for i in range(self.states2):  #final states
               if i==0 and j==0:
                  #skip 0-0 transitions
                  continue
               tmp=FCeqf(self.HR[a], i, j)
               try:
                  FC[a][j*M+i-1]=tmp*FC00*npexp(-(f[a]*j)/self.T)
                  uency[a][j*M+i-1]=(sgnE*E+sgnE*f[a]*(j-i))*self.Hartree2cm_1
               except IndexError:
                  loggingwrite("WARNING: truncated spectrum for mode nr. %d\n"%(a))
                  break
      self.spect=unifSpect(FC, uency, sgnE*E*self.Hartree2cm_1, FC00)

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
    
   def HuangR(self, HRthresh): #what is with different frequencies???
      """ Function that calculates the Huang-Rhys factors for all 
         vibrational states
            
         **PARAMETERS**
         HRthresh threshold for lowest HR-factor (number)
      
         **CALCULATES**
         uniHRall   HR-factors for coinciding frequencies sorted by size (decreasing)
         uniFall    respective frequencies (same order)
   
         **RETURNES**
         nothing. All results are members of class  FC_Spect.
      """
      lenK=len(self.nm.K)
      sortHR=np.zeros(lenK)
      HR=np.zeros(lenK)
      fsort=np.zeros(lenK)
      uniFall=[]
      for j in range(lenK):
         HR[j]=self.nm.K[j]*self.nm.K[j]*self.f[0][j]*.5 
      index=np.argsort(HR, kind='heapsort')
      sortHR=HR[index]
      fsort0=self.f[0][index]
      fsort1=self.f[1][index]
      if np.any(sortHR)<0:
         self.log.write('WARNING: some HR-factors are <0.\
                  In the following their absolute value is used.')
         sortHR=np.abs(sortHR)
      uniHR=[]
      uniF1=[]
      uniF0=[]
      loggingwrite=self.log.write
      # now remove modes with too low frequency;
      # they can easily cause some trouble!
      i=0
      for j in range(len(fsort1)):
         if fsort1[i]*self.Hartree2cm_1<5:
            #print fsort1[i]*self.Hartree2cm_1, sortHR[i]
            loggingwrite(u'\nWARNING: VERY SMALL FREQUENCIES OCCURED.\n'+\
                        '        the mode %d will be disabled.\n'%(i))
            loggingwrite('      freq: %f,  HR-factor: %f \n'%(fsort1[i]*self.Hartree2cm_1, sortHR[i]) )
            fsort1=np.delete(fsort1,i)
            fsort0=np.delete(fsort0,i)
            sortHR=np.delete(sortHR,i)
         else:
            i+=1

      if sortHR[-1]>=10: 
         #if larges HR-factor is too large

         loggingwrite(u'\nWARNING: THE HUANG-RHYS FACTOR SEEMS TO BE'+\
                                                         ' TOO LARGE !!\n')
         loggingwrite(u'        the spectrum will be calculated, but most probably'+\
                                       ' the input-state is inconsistent.\n')
         loggingwrite(u'        Please check the xyz-file created.\n')
         self.makeXYZ()
      # the HR-factors are printed allways.
      loggingwrite(u'HR-fact           freq     delta\n')
      #print(u'HR-fact           freq\n')
      for j in range(len(sortHR)-1,-1,-1):
         #select all 'big' HR-factors 
         if sortHR[j]>=HRthresh:
            uniHR.append(sortHR[j])
            uniF1.append(fsort1[j])
            uniF0.append(fsort0[j])
            loggingwrite(u"%f   %f   %f\n"%(sortHR[j], fsort1[j]*self.Hartree2cm_1, 
                                             np.sqrt(fsort1[j]/fsort0[j]) ))
         else:
            # there will come only smaller ones.
            break
      #keep order: First one is initial, second is final state.
      uniFall.append(uniF0)
      uniFall.append(uniF1)
      # now, make the new results being objects of that class:
      self.HR=uniHR
      self.f=uniFall
      Esign=np.sign(self.Energy[0]-self.Energy[1])
      if Esign==0:
         for i in range(len(self.HR)):
            self.Energy[1]-=self.f[1][i]*self.HR[i]
      else:
         for i in range(len(self.HR)):
            self.Energy[1]+=Esign*self.f[1][i]*self.HR[i]

class CFC_spect(FC_spect):
   """This is more general class compared to FC_spect. Here, both states need to have different
      force constant matrices and the change of the respective frequencies are taken into account
      also for the intensity of the transition.
      The only difference hence is a modified variant of computing the intensities, all other functions
      are inhericed from FC_spect or Spect respectively.
   """
   def __init__(self, f): 
      """The initialisation of this class is the same as for FC_spect
         except for the name.
      """
      self.type='CFC'
      FC_spect.__init__(self, f)

   def calcspect(self):
       """This is used to calculate the line spectrum assuming no mode mixing 
         (shift only)  and coinciding frequencies in both electronic states.
   
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
       # M,N are maximal numbers of vibrational modes (+1, since they should 
       #      be arrived really; count from 0)
       self.states1+=1
       self.states2+=1
   
       def FCchf(HR,N,M,freq):
          """This function calculates the FC-factors for a mode with HR-factor HR,
            frequency freq and the highest achievable states N,M in initial and
            final state respectively. 
            Here, for the overlap the change in frequency between the modes is taken
            into account.
          """
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
   
       n=len(self.HR) #=len(freq)
       setM=False
   
       # correct for vibrational groundstates:
       #E+=(sum(freq[0])-sum(freq[1]))*.5
       #print E, freq[0]
       if self.states2==1: 
          # there was not specified, how many vibr. states in ground-state 
          #           should be taken into account
          setM=True
          M=max(3,int(-1.1*self.HR[0]*self.HR[0]+6.4*self.HR[0]+5.))
       assert n>0, "There is no Huang-Rhys factor larger than the respective"+\
                    " threshold. No mode to be calculated."
       #if setM: the size of these arrays will be overestimated.
       #calculate 0->0 transition
       FC00=1.00
       #print 0,0,0, 10
       npexp=np.exp                  #to accelerates python quite a lot
       freq=self.f

       FC=np.zeros((n,M*self.states1))
       uency=np.zeros((n,M*self.states1)) #frequency
       E=self.Energy[0]-self.Energy[1]
       #here a goes over all modes
       sgnE=np.sign(E)
       if np.sign(E)==0:
          sgnE=1
       for a in xrange(n):
          if setM:
             # set M to fit best to the value at that moment.
             M=max(3,int(-1.1*self.HR[a]*self.HR[a]+6.4*self.HR[a]+5.))
          #print a, HR[a]
          R=FCchf(self.HR[a],self.states1,M,[freq[0][a], freq[1][a]])
          for j in range(self.states1): 
             for i in range(M):
                if i==0 and j==0:
                   #skip 0-0 transitions
                   continue
                tmp=R[i][j]*R[i][j]
                try:
                   FC[a][j*M+i-1]=tmp*FC00*npexp(-(freq[0][a]*j)/self.T)
                   uency[a][j*M+i-1]=sgnE*(E+(freq[0][a]*j-freq[1][a]*i))*self.Hartree2cm_1
                except IndexError:
                   self.log.write("WARNING: truncated spectrum for mode nr. %d\n"%(a))
                   break
       self.spect=unifSpect(FC, uency, sgnE*E*self.Hartree2cm_1, FC00)
   
version='0.1.5'  
#End of FC_Spects.py
