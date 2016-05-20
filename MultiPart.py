#!/usr/bin/python2
# filename: MultiPart.py
import numpy as np
import ctypes as C

# CHANGELOG
# =========
#to version 0.2.0:  
#   a) add ctypes library to make OPA2nPA more efficient
#
#to version 0.1.0:  
#   a) fixed error with self.function and self.n
#   b) added some documentation

class OPAtoNPA:
   """This class gets a spectrum in one-particle approximation of the format
      frequency_1, intensity_1, mode_1  
      frequency_2, intensity_2, mode_2  
           :     ,     :      ,   :  
           :     ,     :      ,   :  
      for all transitions (hence a 3xN-matrix)
      quantity for all transition. 
      The class's interface-functions are:

      init  - gets the number of considered combinations as parameter
      concise - computes an approximate stick-spectrum with fewer transitions
      Calc - computes the n-particle spectrum
      GetSpect - obtains the spectrum as described above. Second parameter:
                 index for cutoff; assumes the spectrum to be sorted by
                 increasing intensity.
   """

   #member variables:
   function = None
   frequency=[]
   intensity=[]
   mode     =[]
   n        =0
   
   # DEFINITION OF FUNCTIONS START
   def __init__(self,n):
      """ when initialising the class for one-particle to n-particle approximation only
         the n is needed and, depending on its size, the function for calculatio is set.
      """
      if n<=0:
         assert 1==2 , "You can not chose less than 1 particle for your approximation method."
      elif n>500:
         assert 1==2 , "you are a bit optimistic. Don't use more than 500 simultaneously changing modes!"
      elif n==1:
         self.function=None
         self.n=n
      elif n==2:
         self.function='OPA2TPA'
         self.n=n
      elif n==3:
         self.function='OPA23PA'
         self.n=n
      else: 
         self.function='OPA2nPA'
         self.n=n
   
   def __OPA2nPA(self,ind00):
      """This function calculates a spectrum with arbitrarily many chaneing modes per transition.
         How ever, the implementation is quite demanding in memory and time, so the number should
         not be chosen to be too high.
      """
      def allnonzero(foo):
         """ a more efficient version of np.all(foo!=0), made for arrays and integers...
         """
         try:
            foo=np.array(foo)
            for s in range(len(foo)):
               if foo[s]==0:
                  return False
         except TypeError:
            if foo==0:
               return False
         return True
      
      def putN(intens, freq, mode, n):
         """Function wrapper for the C-function that calculates the combination transition
            energies. I expect a large speedup and some spare memory from this construction.
         """
         # add the corresponding C-libraries:
         #LIBRARY_PATH = './putN.so'
         lib = C.CDLL('/home/tm162/bin/smallscript/libputN.so')
         length= C.c_int(len(intens))
         print freq
         print freq.ctypes.data_as(C.POINTER(C.c_double))[0]
         TPA=lib.putN( intens.ctypes.data_as(C.POINTER(C.c_double)), 
                      freq.ctypes.data_as(C.POINTER(C.c_double)), 
                      mode.ctypes.data_as(C.POINTER(C.c_int)),
                      C.c_int(int(n)),length)
         TPAintens=np.array(np.fromiter(TPA, dtype=np.float64, count=length))[0]
         TPAfreq=np.array(np.fromiter(TPA, dtype=np.float64, count=length))[1]
         print TPA.ctypes.data_as(C.POINTER(C.c_double))[0][0]
         print freq[0]

         return intens, freq

      length=len(self.frequency)
      freq00=self.frequency[ind00]
      intens00=self.intensity[ind00]
      for i in range(length):
         self.frequency[i]-=freq00
         self.intensity[i]/=intens00
      #newmode=np.zeros((1,len(self.mode))) #for n>1: matrix-structure needed
      #newmode[0]=self.mode
      x=self.mode.max()
      if self.n>x:
         self.n=x
         #save complexity, that it does not try to make more combinations 
         # than actually possible...
      
      TPAfreq, TPAintens=putN(self.intensity, self.frequency, self.mode, self.n)

      #TPAfreq, TPAintens=putN( byref(opa), self.n)
      for i in xrange(len(TPAfreq)):
         TPAfreq[i]+=freq00
         TPAintens[i]*=intens00
      return TPAfreq, TPAintens
   
   def __OPA2TPA(self,ind00):
      """ This function computes combination transit with up to two modes changeing at once.
      """
      length=len(self.intensity)
      TPAfreq=np.zeros((length+1)*(length+2)//2+length+1) #this is overestimation of size...
      TPAintens=np.zeros((length+1)*(length+2)//2+length+1)
      ind=0
      intens00 = self.intensity[ind00]
      freq00 = self.frequency[ind00]
      #now, make a mapping of names for better performance:
      intens=self.intensity
      freq=self.frequency
      mode=self.mode
      intens00 = self.intensity[ind00]
      freq00 = self.frequency[ind00]
      for i in range(length):
         if mode[i]==0:
            #0-0 transition is included, but only once!
            continue
         TPAintens[ind]=intens[i] #this is OPA-part
         TPAfreq[ind]=freq[i]
         ind+=1
         for j in range(i+1,length):
            #not only same mode but all modes with lower number should not be taken into account here!?
            if mode[i]==mode[j] or mode[j]==0: 
               #both have same mode... or mode[j] is 0-0 transition
               continue
            if intens[i]*intens[j]<intens00*intens00*0.0001:
               #save memory by not saving low-intensity-modes
               continue
            TPAintens[ind]=intens[i]*intens[j]/intens00
            TPAfreq[ind]=freq[i]+freq[j]-freq00
            ind+=1
      index=np.argsort(TPAfreq,kind='heapsort')
      TPAfreq=TPAfreq[index]
      TPAintens=TPAintens[index]
      return TPAfreq, TPAintens

   def __OPA23PA(self,ind00):
      """ This function computes combination transit with up to three modes changeing at once.
      """
      #FIRST initialise some quantities making later synatax 
      # more convenient
      length=len(self.frequency)
      TPAfreq=[]
      TPAintens=[]
      intens00 = self.intensity[ind00]
      freq00 = self.frequency[ind00]
      #now, make a mapping of names for better performance:
      intens=self.intensity
      freq=self.frequency
      mode=self.mode
      TPAintens.append(intens00) #this is OPA-part
      TPAfreq.append(freq00)
      #SECOND go through the whole spectrum (besides 0-0) and compute all combinations 
      #        besides self-combination
      for i in range(length):
         if mode[i]==0:
            #0-0 transition is included, but only once!
            continue
         TPAintens.append(intens[i]) #this is OPA-part
         TPAfreq.append(freq[i])
         # here the combination part starts
         for j in range(i+1,length):
            #both have same mode... or mode[j] is 0-0 transition
            if mode[i]==mode[j] or mode[j]==0: 
               continue
            if intens[i]*intens[j]<intens00*intens00*0.00001:
               #save memory by not saving low-intensity-modes
               continue
            TPAintens.append(intens[i]*intens[j]/intens00)
            TPAfreq.append(freq[i]+freq[j]-freq00)
            # look for all combinations of combinations for three-particle approx.
            for k in range(j+1,length):
               if mode[k]==mode[j] or mode[k]==0: #both have same mode...
                  continue
               if mode[k]==mode[i]:
                  continue
               if intens[i]*intens[j]*intens[k]<\
                                      (intens00*intens00)*intens00*0.0001:
                  #save memory by not saving low-intensity-modes
                  continue
               TPAintens.append(intens[i]*intens[j]*intens[k]/
                                                    (intens00*intens00))
               TPAfreq.append(freq[i]+freq[k]+freq[j]-2*freq00)
      #FINALLY save the spectrum into numpy-matrices
      freq=np.zeros(len(TPAfreq))
      intens=np.zeros(len(TPAintens))
      #this can not be done by np.matrix() due to dimensionality...
      for i in xrange(len(freq)): 
           freq[i]=TPAfreq[i]
           intens[i]=TPAintens[i]
      return freq, intens

   def concise(self,broadness):
      """ This function shrinks length of the stick-spectrum to speed-up the 
         calculation of broadened spectrum (folding with lineshape-function).
         It puts all transitions within a tenth of the Gaussian-width into one line.
      
         ==PARAMETERS==
         broadness:   gamma from the Lorentian-courve; specifying, 
                     how many lines will be put together
      
         ==RETURNS==
         intens2:     shrinked intensity-vector, sorted by increasing frequency
         freq2:       shrinked frequency-vector, sorted by increasing frequency
      """
      # both arrays are frequency-sorted 
      #initialize arrays
      intens2=[]
      freq2=[]
      mini=0
      #go through spectrum and sum everything up.
      freq=self.frequency
      intens=self.intensity
      for i in range(len(freq)):
         # index-range, put into one line is broadness/5; this should be enough
         tempintens=0
         for j in xrange(mini, len(freq)):
            tempintens+=intens[j]
            if freq[j]>=freq[mini]+broadness/5.:
               mini=j # set mini to its new value
               intens2.append(tempintens)
               freq2.append(freq[j]-broadness/10.)
               break
      return freq2, intens2

   def GetSpect(self,linspect, minint=0):
      """ This function needs to be called before calculating the combination spectrum.
         It copies the quantites needed ( frequency, intensity and the number of respective mode
         that changes) to members of the class OPAtoNPA.
      """
      self.frequency=linspect[0][minint:]
      self.intensity=linspect[1][minint:]
      self.mode     =linspect[2][minint:]
  
   def Calc(self):
      """This function interfaces the computation of combination transitions to the respective spectrum
         that is calculated.
      """
      #first, get the index of the purely electronic transition
      ind=self.mode.argmin()
      if self.function==None:
         return self.intensity, self.frequency
      elif self.function=='OPA2TPA':
         #now, calculate the full spectrum in 2-particle approx.
         TPAfreq, TPAintens=self.__OPA2TPA(ind)
      elif self.function=='OPA23PA':
         #or three-particle approx. 
         TPAfreq, TPAintens=self.__OPA23PA(ind)
      elif self.function=='OPA2nPA':
         # or any other number of particles.
         TPAfreq, TPAintens=self.__OPA2nPA(ind)
      # finally sort and truncate the full spectrum. (by low-intensity-transitions)
      index=np.argsort(TPAintens,kind='heapsort')
      TPAintens=TPAintens[index] #resort by intensity
      TPAfreq=TPAfreq[index]
      minint=0
      for i in range(len(TPAintens)):
         if TPAintens[i]>=1e-6*TPAintens[-1]:
            minint=i
            break
      return TPAintens[minint:], TPAfreq[minint:]
   # DEFINITION OF FUNCTIONS END
   # end of class definition.

version='0.2.0'  
# End of MultiPart.py
