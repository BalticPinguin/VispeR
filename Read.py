#!/usr/bin/python2
# filename: Read.py
import numpy as np
from scipy.linalg.lapack import dsyev
import re, mmap, os.path, math, sys

# CHANGELOG 
# =========
#to version=0.1.6   
# 1. Changed the Estring of G09 again; hopefully this is the best 
#    way reading it.  
#
#to version=0.1.5   
# 1. Reduced amount of redundant code; added strings to rtype.  
# 2. Gradients for NWChem work now. Gradient() returns now only
#    final states gradient.  
# 3. repaired Estring of G09  
# 4. repaired gradpolishstring of g09_fchk  
#
#to version=0.1.0    
# 1. Resorted the class: reads data-based, not logfile-based.
#    This might be slower at runtime but is better structured
#    and easier to read. Moreover, it avoids redundant code.  
# 2. Added class rtype that knows the strings to search for
#    some of the data.  
#  
#to version=0.0.1  
# 1. I deleted the functions reading L and f from the log-file.  
#    They are available only in older versions (non-class)
#    of smallscript.  
# 2. Changed functions _Grad to read only the gradient and added 
#    function ReadGrad  
# 3. Fixed error with units of coordinates in ReadNWChem  
# 4. refixed energy in Gaussian. Need to do it for NWChem as well...  

# The Class Read
# ==============
class Read: 
   """This class is made to outsource all work (and much code) related to reading
      data from log-files in one of the supported formats. 
      At the moment, for each format there are two functions: One that reads all data
      required for the models FC,CFC,DR,URDR and one function that reads the gradient
      as required for the class GFC.
   
      Description of the class
      -----------------------
      It mainly consists of the following three members:
      the name of the file with **initial** state and that of
      **final** state. The third member is **type** which
      is ' ' if the files differ in their format.
      In that case than, the two members *ftype* and *itype
      replace it.
   """
   #members:
   AMU2au=1822.88839
   Angs2Bohr=1/0.52917721092
   ev2Hartree = 0.0367493081366
   init='unknown'
   final='unknown'
   ftype=False
   itype=False
   dim=0

   class rtype:
      """ inline-class for Read. It provides the strings for the
          read-functions to make them more general and reduce the amount
          of redundant code.
      """
      def __init__(self, typestring):
         self.type=typestring
         self.mass=0
         if typestring=="NWChem":
            self.gradString=r"(?<=atom               coordinates                        gradient)[\d\w\- .\n]+"
            self.gradPolishString=r"[-\d.]+[ ]+[-\d.]+[ ]+[-\d.]+(?=\n)"
            self.Estring=r'(?<= Excited state energy = =)[\-\d ]+'
            self.ForceString=r"(?<=MASS-WEIGHTED PROJECTED HESSIAN \(Hartree/Bohr/Bohr/Kamu\)\n)[\dD \.\-\+\n]+"
            self.MassString=r"(?<=atom    #        X              Y              Z            mass\n)[\d\w\-\+ \.\n]+"
            self.MassString2=r"[ -]\d\.[\dD\.\+\-]+(?=\n)"
            self.CoordS=r"(?<=atom    #        X              Y              Z            mass\n)[\d\w\-\+ \.\n]+"

         if typestring=="G09":
            self.gradString=r"Number     Number              X              Y              Z\n [\-\.\d \n]+"
            self.gradPolishString=r"[\-\d\.]+[\d]{9}"
            self.Estring=r"(?<=SCF Done:  E)[\d.\n \(\)\w\+\-=]+"
            self.Estring2=r"(?<=SCF Done:  E)[\d.\n \(\)\w\+\-=]+"
            self.ForceString=r"Force constants in Cartesian coordinates: [\n\d .+-D]+"
            self.MassString=r"AtmWgt= [\d .]+"
            self.CoordS=r' Number     Number       Type             X           Y           Z[\n -.\d]+'

         if typestring=="GAMESS":
            assert 1==2, "GAMES-files are not supported at the moment. I am sorry."
            self.gradString=""
            self.gradPolishString=""
            self.Estring=r'(?<=TOTAL ENERGY =)[\-\d. ]+'
            self.ForceString=r"(?<=CARTESIAN FORCE CONSTANT MATRIX\n)[\d\w \.\-\n]+"
            self.MassString=r"(?<=ATOMIC WEIGHTS \(AMU\))[ \.\w\n]+"
            self.MassString2=r'[\d.]{7,10}(?=\n)'
            self.CoordS=r'(?<=ATOM             X              Y              Z\n\n)[\n -.\d\w]+'

         if typestring=="G09_fchk":
            self.gradString=r"(?<=Cartesian Gradient)[R N=\d.\+\- E\n]+"
            self.gradPolishString=r"[ -]\d\.[\d\-\+E]+"
            self.Estring=r'(?<=Total Energy                               R  )[ \-\d.E\+]+'
            self.ForceString=r"(?<=Cartesian Force Constants   )[RN\=\d\+\-E \.\n]+"
            self.MassString=r"(?<=Real atomic weights     )[RN\=\-\+\d\. E\n]+"
            self.CoordS=r'Current cartesian coordinates  [RN\=\n \-\+.\d E]+'

   #attributes:
   def __init__(self, initial, final):
      """ initialises the class Read. It will get the following properties:
         init: the name of the file with initial state information. 
               Having this as class member, we don't can access it any time we need.  
         final: name of file with final states information.  
         type/ itype+ftype : the type of file, meaning: its format. With the
               possiblitity of having itype and ftype defined, the two states
               are allowed to have different formats (which is NOT recommended).
      """
      #initialise the class members:
      self.init=initial
      self.final=final
      #check, which file type they have. Do it individually
      # for both to allow for the (not very sensible) case
      # that initial state and final state have different format.
      if (".fchk" in self.init) or (".FChk" in self.init): 
         self.itype =self.rtype("G09_fchk")
      if (".fchk" in self.final) or (".FChk" in self.final) :
         self.ftype=self.rtype('G09_fchk')
      
      #Besides the formcheck-files, all log-files are recognised by a
      # certain line/string in the following:
      if not self.itype:
         with open(self.init, "r+b") as f: #open file as mmap
            mapping = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
            for line in iter(mapping.readline, ""): #go through every line and test characteristic part
               if "GAMESS" in line: #if it is found: read important quantities from file
                  self.itype =self.rtype("GAMESS")
                  break
               elif "Northwest Computational Chemistry Package (NWChem)" in line:
                  self.itype =self.rtype("NWChem")
                  break
               elif "Gaussian(R)" in line:
                  self.itype=self.rtype("G09")
                  break
            #  There is some error: File not recognised or an other programme was used.
            else: 
               print "file type not recognised"

      if not self.ftype:
         with open(self.final, "r+b") as f: #open file as mmap
            mapping = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
            for line in iter(mapping.readline, ""): #go through every line and test characteristic part
               if "GAMESS" in line: #if it is found: read important quantities from file
                  self.ftype=self.rtype("GAMESS")
                  break
               elif "Northwest Computational Chemistry Package (NWChem)" in line:
                  self.ftype=self.rtype("NWChem")
                  break
               elif "Gaussian(R)" in line:
                  self.ftype=self.rtype("G09")
                  break
            #  There is some error: File not recognised or an other programme was used.
            else: 
               print "file type not recognised"
        # End of __init__()

   def __Read_Mass(self, logfile, rtype):
      files=open(logfile, "r")
      log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close()

      atmwgt=re.findall(rtype.MassString, log)
      if rtype.type=='G09':
         if atmwgt==[]:
            #it may happen, that the masses are not given in the above format. Then let's see,
            #if we can recognise some other format.
            atmwgt=re.findall(r"(?<= and mass  )[\d. ]+", log)
            if atmwgt==[]:
               # sometimes this format is not available as well.
               # So just leave it away and take masses of other file.
               return [0,0]
         # Determine atomic masses in a.u. Note mass contains sqrt of mass!!!
         mtemp=[]
         for j in range(len(atmwgt)/2): # because atomic masses are printed twize in log-files...
            mtemp.append(re.findall(r'[\d.]+',atmwgt[j]))
         dim=0
         for j in range(len(mtemp)):
               dim+=len(mtemp[j]) # dim will be sum over all elements of temp
         dim*=3
         
         temp=re.findall(r' Number     Number       Type             X           Y           Z[\n -.\d]+', log)
         tmp=re.findall(r'[ -][\d]+.[\d]+', temp[-1])
         if dim!=len(tmp):
            # this is necessary since they are not always printed twice...
            #atmwgt=re.findall(r"AtmWgt= [\d .]+",log) -> do not look for the same twice!!
            mtemp=[]
            for j in range(len(atmwgt)):
               mtemp.append(re.findall(r'[\d.]+',atmwgt[j]))
            dim=0
            for j in range(len(mtemp)):
                  # dim will be sum over all elements of temp
                  dim+=len(mtemp[j])
            dim*=3
         #if still something is wrong with the dimensionality:
         assert len(tmp)==dim, 'Not all atoms were found! Something went wrong...{0}, {1}'.format(len(tmp),dim)

         foonum=0
         mass=np.zeros(dim/3) # this is an integer since dim=3*N with N=atomicity
         for j in range(len(mtemp)):
            for k in range(len(mtemp[j])):
               mass[k+foonum]=np.sqrt(float(mtemp[j][k])*self.AMU2au) #elements in m are sqrt(m_i) where m_i is the i-th atoms mass
            foonum+=len(mtemp[j])
         assert not np.any(mass==0) , "some atomic masses are zero. Please check the input-file!\n {0}".format(mass)

      elif rtype.type=='G09_fchk':
         # Determine atomic masses in a.u. Note mass contains sqrt of mass!!!
         assert len(atmwgt)>0, "Formcheck-file not complete. No masses available.\n Searched for %s"%rtype.MassString
         mtemp=re.findall(r'[\d.]+E[+-][\d]+',atmwgt[0])
         dim=int(re.findall("(?<=N\= )[\d ]+", atmwgt[0])[0])
         dim*=3
   
         foonum=0
         mass=np.zeros(dim/3) # this is an integer since dim=3*N with N=atomicity
         for j in range(len(mtemp)):
            mass[j]=np.sqrt(float(mtemp[j])*self.AMU2au) #elements in m are sqrt(m_i) where m_i is the i-th atoms mass
         assert not np.any(mass==0) , "some atomic masses are zero. Please check the input-file! {0}".format(mass)

      elif rtype.type=='NWChem' or rtype.type=='GAMESS':
         # Determine atomic masses in a.u. Note mass contains sqrt of mass!
         if atmwgt==[]:
            # this is the case especially, if no frequency-calculation was performed.
            return 0
         mtemp=re.findall(rtype.MassString2, atmwgt[-1]) #only weights remain
         dim=3*len(mtemp)
   
         mass=np.zeros(dim/3) # this is an integer since dim=3*N with N=atomicity
         for j in range(len(mtemp)):
            mass[j]=np.sqrt(float(mtemp[j].replace('D','e'))*self.AMU2au) #mass[j] is square-root of masses
         assert not np.any(mass==0) , "some atomic masses are 0. Please check the input-file! {0}".format(mass)

      self.dim=dim
      rtype.mass=mass
      return mass

   def __Read_Force(self, logfile, rtype):
      files=open(logfile, "r")
      log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close()

      dim=self.dim
      mass=rtype.mass
      # Reading of Cartesian force constant matrix  
      f=re.findall( rtype.ForceString, log, re.M)
      if rtype.type=='G09':
         #in case of G09-files, test an other string:
         if f==[]:
            #try to find matrix from option "Force"
            f=re.findall(r"The second derivative matrix:[XYZ\n\d .-]+", log, re.M)
            #else: return 0; no force constant matrix given.
            second=True
         else:
            second=False

      #Now, do general preparation for the FC-matrix:
      # Reading of Cartesian force constant matrix  
      if f==[]:
         return 0
      f_str=str([f[-1]])
      lines=f_str.strip().split("\\n")
      # Reading of Cartesian force constant matrix  
      F=np.zeros((dim,dim))
      n=0
      k=0

      if rtype.type=='G09':
         line=0
         if second:
            assert len(lines)>2, "Can't find any forces. Force constants given in unknown format. \n Searched for %s"%rtype.ForcString
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
            return F
         else:
            for i in xrange(2,len(lines)):
               if i == dim+k-5*n+2: 
                  k=i-1
                  n+=1
                  continue
               elements=lines[i].replace('D','e').split()
               for j in range(1,len(elements)):
                  F[int(elements[0])-1][j-1+5*n]=float(elements[j])
                  F[j-1+5*n][int(elements[0])-1]=float(elements[j])

      elif rtype.type=='G09_fchk':
         m=0
         for i in xrange(1,len(lines)): # the first line is not part of the matrix
            elements=lines[i].replace('E','e').split()
            if n==len(F): # break this outer loop if inner is broken
               break
            for j in range(len(elements)):
               F[k][n]=F[n][k]=float(elements[j])
               if k==n:
                  n+=1
                  if n==len(F): # matrix is full; there will be only waste
                     break
                  k=0
               else:
                  k+=1
      elif rtype.type=='NWChem':
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
               F[int(elements[0])-1][j-1+10*n]=float(elements[j])/(1000*self.AMU2au)
               F[j-1+10*n][int(elements[0])-1]=float(elements[j])/(1000*self.AMU2au)
         return F

      elif rtype.type=='GAMESS':
         lines=lines[5:]
         #lines include each three six elements.
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

      #mass-weight the Force-constant matrix. This is reached for all 
      # routines where no return-statement is given in place.
      for i in range(0,dim):
         for j in range(0,dim):
            F[i][j]/= (mass[i//3]*mass[j//3]) 
      return F

   def __Read_Coords(self, logfile, rtype):
      files=open(logfile, "r")
      log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close()

      dim=self.dim
      # Reading Cartesian coordinates
      temp=re.findall(rtype.CoordS, log)
      Coord=np.zeros((3,dim//3))
      
      if rtype.type=='G09':
         tmp=re.findall(r'[ -][\d]+.[\d]+', temp[-1])
         for j in range(len(tmp)):
            Coord[j%3][j/3]=tmp[j]
         Coord*=self.Angs2Bohr
   
      elif rtype.type=='G09_fchk':
         tmp=re.findall(r'[ -][\d.]+E[+-][\d]+', temp[0])
         # Reading Cartesian coordinates
         for j in range(len(tmp)):
            Coord[j%3][j/3]=float(tmp[j].replace('E','e')) # need to convert by some factor!!

      elif rtype.type=='NWChem':
         if temp!=[]:
            tmp=re.findall(r'[ -]\d\.[\dD\.\+\-]+', temp[-1])
            k=0
            for j in range(dim//3):
               for i in range(3):
                  Coord[i][j]=float(tmp[i+3*j+k].replace('D','e'))
               k+=1 # important to skip masses in tmp 
         else:
            # this is coordinate set, if no frequencies are calculated.
            # Maybe I should use this in all cases.
            temp=re.findall(r"No.       Tag          Charge          X              Y              Z[\d\w\-\+ \.\n]+",log)
            tmp=re.findall(r'[ -]\d\.[\d]{5,10}', temp[-1])
            for j in range(dim//3):
               for i in range(3):
                  Coord[i][j]=float(tmp[i+3*j])
            # depending on the units: need to convert to a.u.!
            if re.search(r"Output coordinates in angstroms ", log) is not None:
               # need them in atomic units later:
               Coord*=self.Angs2Bohr 

      elif rtype.type=='GAMESS':
         for i in range(3):
            for j in range(dim//3):
               Coord[i][j]=float(tmp[3*i+j])
         Coord*=self.Angs2Bohr
      return Coord

   def __Read_Energy(self, logfile, rtype):
      files=open(logfile, "r")
      log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close()

      Etemp=re.findall(rtype.Estring, log, re.M)
      if rtype.type=='G09':
         if Etemp==[]:
            Etemp=re.findall(rtype.Estring2, log, re.M)
            assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.\n Searched for %s'%rtype.Estring2
            #this is only for Gaussian-case but doesn't disturb otherwise..
            if re.search(r'\n ', Etemp[-1]) is not None:
               Etemp[-1]=Etemp[-1].replace("\n ", "") 
            E=-float(re.findall(r'[\d.]+', Etemp[-1])[0]) #energy is negative (bound state)
         else:
            #this number is already given in Hartree.
            if len(Etemp)>1:
               E=float(Etemp[-2].split()[2])
            else:
               E=float(Etemp[0].split()[2])
      elif rtype.type=='G09_fchk' or rtype=='GAMESS':
         assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.\n Searched for %s'%rtype.Estring
         # replacement has only effect in case of G09_fchk but doesn't disturb for GAMESS.
         E=float(Etemp[-1].replace('E','e'))
      elif rtype.type=='NWChem':
         #Etemp is energy of the (excited) state
         if Etemp==[]:
            #if the ground state is calculated
            Etemp=re.findall(r'(?<=Total DFT energy =)[\-\d. ]+', log, re.M)
         assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.\n Searched for %s'%rtype.Estring
         #FIXME: need to fix it for tddft-case.
         assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.\n Searched for %s'%rtype.Estring
         E=float(Etemp[-1]) 
      return E

   def __Read_Grad(self, logfile, rtype):
      """ This function reads the required quantities from logfiles
         of variable format and returns the gradient, if given.
         **PARAMETERS**
         logfile:       specifies the name of the log-file
         rtype :        specifies the file-type.

         **RETURNS**
         grad:        gradient of the PES of excited state at ground state equilibrium-geometry
      """
      files=open(logfile, "r")
      log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close()

      grad=re.findall(rtype.gradString, log)
      assert len(grad)>0, "ERROR: No gradient given. Searched for: \n      %s"%rtyp.gradString
      if rtype.type=='NWChem':
         # in this case, I need to do it differently:
         Grad=re.findall(rtype.gradPolishString, grad[0])
         grad=np.zeros(len(Grad)*3)
         for i in xrange(len(Grad)):
            tmp=Grad[i].split()
            for j in [0,1,2]:
               grad[3*i+j]=float(tmp[j])
         return grad
      #if self.type!='NWChem': 
      Grad=re.findall(rtype.gradPolishString, grad[-1])
      grad=np.zeros(len(Grad))
      for i in xrange(len(Grad)):
         grad[i]=float(Grad[i])
      return grad

   # END: FUNCTIONS TO READ THE NEEDED DATA FROM OUTPUT-FILES IN DIFFERENT FORMATS.

   #USER-FUNCTIONS ASKING FOR CERTAIN DATA:
   def Gradient(self):
      # the gradient makes only sense when given for the final state.
      return self.__Read_Grad(self.final, self.ftype)
   
   def Energy(self):
      return [self.__Read_Energy(self.init, self.itype) ,self.__Read_Energy(self.final, self.ftype)]

   def mass(self):
      return [self.__Read_Mass(self.init, self.itype), self.__Read_Mass(self.final, self.ftype)]
   
   def Coordinates(self):
      return [self.__Read_Coords(self.init, self.itype), self.__Read_Coords(self.final, self.ftype)]

   def Force(self):
      return [self.__Read_Force(self.init, self.itype), self.__Read_Force(self.final, self.ftype)]
   #END USER-FUNCTIONS ASKING FOR CERTAIN DATA:

version='0.1.6'
# End of Read.py
