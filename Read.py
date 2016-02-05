#!/usr/bin/python2
# filename: Read.py
import numpy as np
from scipy.linalg.lapack import dsyev
import re, mmap, os.path, math, sys

# ============ CHANGELOG =================
# I deleted the functions reading L and f from the log-file.
# They are available only in older versions (non-class)
# of smallscript.


#  FUNCTIONS TO READ THE NEEDED DATA FROM OUTPUT-FILES IN DIFFERENT FORMATS.
class Read: 
   #members:
     init='unknown'
     final='unknown'
     type='unknown'
   #attributes:
   def __init__(self, initial,final)
      """ initialises the class Read. It will get the following properties:
       init:  the name of the file with initial state information. 
              Having this as class member, we don't can access it any time we need.
       final: name of file with final states information.
       type, itype+ftype : the type of file, meaning: its format. With the 
                           possiblitity of having itype and ftype defined, the two states
                           are allowed to have different formats (which is NOT recommended).
      """
      self.init=initial
      self.final=final
      if (".fchk" in self.init))\ or (".FChk" in self.init)): 
         self.type = "G09_fchk"
      if (( ".fchk" in self.final) or (( ".FChk" in self.final) :
         if self.type!="G09_fchk":
            self.type=""
            self.ftype='G09_fchk'

      if self.type!='G09_fchk':
         with open(self.initial, "r+b") as f: #open file as mmap
            mapping = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
            for line in iter(mapping.readline, ""): #go through every line and test characteristic part
               if "GAMESS" in line: #if it is found: read important quantities from file
                  self.type ="GAMESS"
                  break
               elif "Northwest Computational Chemistry Package (NWChem)" in line:
                  self.type ="NWChem"
                  break
               elif "Gaussian(R)" in line:
                  self.type="G09"
                  break
            #  There is some error: File not recognised or an other programme was used.
            else: 
               print "file type not recognised"

         with open(self.final, "r+b") as f: #open file as mmap
            mapping = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
            for line in iter(mapping.readline, ""): #go through every line and test characteristic part
               if "GAMESS" in line: #if it is found: read important quantities from file
                  if self.type !="GAMESS":
                     self.itype=self.type
                     self.type=""
                     self.ftype="GAMESS"
                  break
               elif "Northwest Computational Chemistry Package (NWChem)" in line:
                  if self.type !="NWChem":
                     self.itype=self.type
                     self.type=""
                     self.ftype="NWChem"
                  break
               elif "Gaussian(R)" in line:
                  if self.type!="G09":
                     self.itype=self.type
                     self.type=""
                     self.ftype="G09"
                  break
            #  There is some error: File not recognised or an other programme was used.
            else: 
               print "file type not recognised"
   
   def ReadAll(state):
      """This function is supposed to be the main interface
      for accessing this class. """
      #first chose the correct file and its format
      if state=='i':
         if self.type=='':
            kind=self.itype
         else:
            kind=self.type
         inputfile=self.init
      else:
         if self.type=='':
            kind=self.ftype
         else:
            kind=self.type
         inputfile=self.init

      # than perform the correct calculation with it and return all data.
      if kind == "G09_fchk":
         return self.__ReadG09_fchk(inputfile)
      if kind =="GAMESS":
         return self.__ReadGAMESS(inputfile)
      if kind =="NWChem":
         return self.__ReadNWChem(inputfile)
      if kind =="G09":
         return self.__ReadG09(inputfile)
      else:
         print "errer in read-options."
         return 2

   def __ReadG09(self, fileN):
      """ This function reads the required quantities from the gaussian-files

      **PARAMETERS**
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
      if self.logging[0]<2:
         self.logging[1].write("Number of atoms: {0}\nNumber of vibrational "
                           "modes: {1} \n Sqrt of masses in a.u. as read from log file\n{2}\n".format(dim/3,dim,mass))
      # Reading Cartesian coordinates
      Coord=np.zeros((3, dim/3))
      for j in range(len(tmp)):
         Coord[j%3][j/3]=tmp[j]

      Coord*=self.Angs2Bohr
      if self.logging[0]<1:
         self.logging[1].write("Cartesian coordinates (a.u.) \n{0}\n".format(Coord.T))

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
      if self.logging[0]<1:
         self.logging[1].write('F matrix as read from log file\n{0} \n'.format(F))
      for i in range(0,dim):
         for j in range(0,dim):
            F[i][j]/= (mass[i//3]*mass[j//3]) 
      Etemp=re.findall(r'HF=-[\d.\n ]+', log, re.M)
      assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
      if re.search(r'\n ', Etemp[-1]) is not None:
         Etemp[-1]=Etemp[-1].replace("\n ", "") 
      if self.logging[0]<=1:
         self.logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
      E=-float(re.findall(r'[\d.]+', Etemp[-1])[0]) #energy is negative (bound state)
      return dim, Coord, mass, F, E

   def __ReadG09_Grad(self, final):
      """ This function reads the required quantities from the gaussian-files

      **PARAMETERS**
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
      if self.logging[0]<=1:
         self.logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
      E=-float(re.findall(r'[\d.]+', Etemp[-1])[0]) #energy is negative (bound state)
     
      grad=re.findall(r"Number     Number              X              Y              Z\n [\-\.\d \n]+",log)
      Grad=re.findall(r"[\-\d\.]+[\d]{9}", grad[0])
      grad=np.zeros((len(Grad),1))
      for i in xrange(len(Grad)):
         grad[i]=float(Grad[i])
      return grad, E

   def __ReadGO9_fchk(self, fileN):
      """ This function reads the required quantities from the gaussian-files

      **PARAMETERS**
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
      if self.logging[0]<2:
         self.logging[1].write("Number of atoms: {0}\nNumber of vibrational "
                           "modes: {1} \n Sqrt of masses in a.u. as read from log file\n{2}\n".format(dim/3,dim,mass))
      # Reading Cartesian coordinates
      Coord=np.zeros((3, dim/3))
      for j in range(len(tmp)):
         Coord[j%3][j/3]=float(tmp[j].replace('E','e')) # need to convert by some factor!!

      #Coord*=self.Angs2Bohr
      if self.logging[0]<1:
         self.logging[1].write("Cartesian coordinates (a.u.) \n{0}\n".format(Coord.T))

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
      if self.logging[0]<1:
         self.logging[1].write('F matrix as read from formcheck-file\n{0} \n'.format(F))
      for i in range(0,dim):
         for j in range(0,dim):
            F[i][j]/= (mass[i//3]*mass[j//3]) 

      Etemp=re.findall(r'(?<=Total Energy                               R  )[ \-\d.E\+]+', log, re.M)
      assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
      Etemp=float(Etemp[0].replace('E','e'))
      if self.logging[0]<=1:
         self.logging[1].write('temporary energy of state: {0}\n'.format(Etemp))
      return dim, Coord, mass, F, Etemp

   def __ReadG09_fchk_Grad(self, final):
      """ This function reads the required quantities from the gaussian-files

      **PARAMETERS**
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
      if self.logging[0]<=1:
         self.logging[1].write('temporary energy of state: %d\n'%(E))
     
      #get gradient
      grad=re.findall(r"(?<=Cartesian Gradient)[R N=\d.\+\- E\n]+", log)
      Grad=re.findall(r"[\-\d\.\+E]+", grad[0])[1:]
      grad=np.zeros(len(Grad))
      for i in xrange(len(Grad)):
         grad[i]=float(Grad[i])
      return grad, E

   def __ReadGAMESS(self, fileN):
      """ This function reads the required quantities from the GAMESS-log-files

      **PARAMETERS**
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
      if self.logging[0]<2:
         self.logging[1].write("Number of atoms: {0}\nNumber of vibrational "
                           "modes: {1} \n Sqrt of masses in a.u. as read from log file\n{2}\n".format(dim/3,dim,mass))
      
      # Reading Cartesian coordinates
      temp=re.findall(r'(?<=ATOM             X              Y              Z\n\n)[\n -.\d\w]+', log)  
                 #+ need even more in structure
      tmp=re.findall(r'[ -][\d]+.[\d]+', temp[-1])
      Coord=np.zeros((3,dim//3))
      for i in range(3):
         for j in range(dim//3):
            Coord[i][j]=float(tmp[3*i+j])
      Coord*=self.Angs2Bohr
      if self.logging[0]<1:
         self.logging[1].write("Cartesian coordinates (a.u.) \n{0}\n".format(Coord.T))
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

      if self.logging[0]<1:
         self.logging[1].write('F matrix as read from log file\n{0} \n'.format(F))
      for i in range(0,dim):
         for j in range(0,dim):
            F[i][j]/= (mass[i//3]*mass[j//3]) 

      #get energy of the state
      Etemp=re.findall(r'(?<=TOTAL ENERGY =)[\-\d. ]+', log, re.M)
      assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
      if self.logging[0]<=1:
         self.logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
      E=float(Etemp[-1]) #energy is negative (bound state)
      return dim, Coord, mass, F, E

   def __ReadGAMESS_Grad(self, final):
      """ This function reads the required quantities from the GAMESS-log-files

      **PARAMETERS**
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
      if self.logging[0]<=1:
         self.logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
      E=float(Etemp[-1]) #energy is negative (bound state)
      #grad=re.findall(r"Number     Number              X              Y              Z\n [\-\.\d \n]+",log)
      #Grad=re.findall(r"[\-\d\.]+[\d]{9}", grad[0])
      #grad=np.zeros((len(Grad),1))
      #for i in xrange(len(Grad)):
         #grad[i]=float(Grad[i])
      return grad, E

   def __ReadNWChem(self, fileN):
      """ This function reads the required quantities from the NWChem-files

      **PARAMETERS**
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
      if self.logging[0]<2:
         self.logging[1].write("Number of atoms: {0}\nNumber of vibrational "
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
      Coord*=self.Angs2Bohr
      if self.logging[0]<1:
         self.logging[1].write("Cartesian coordinates (a.u.) \n{0}\n".format(Coord.T))

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
      if self.logging[0]<1:
         self.logging[1].write('F matrix as read from log file\n{0} \n'.format(F))

      #get energy of the state
      Etemp=re.findall(r'(?<=Total DFT energy =)[\-\d. ]+', log, re.M)
      assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
      if self.logging[0]<=1:
         self.logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
      E=float(Etemp[-1]) #energy is negative (bound state)
      return dim, Coord, mass, F, E

   def __ReadNWChem2_Grad(self, final):
      """ This function reads the required quantities from the NWChem-files

      **PARAMETERS**
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
      if self.logging[0]<=1:
         self.logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
      E=float(Etemp[-1]) #energy is negative (bound state)
      return grad, E

# END: FUNCTIONS TO READ THE NEEDED DATA FROM OUTPUT-FILES IN DIFFERENT FORMATS.
   
#version=0.0.1
# End of Read.py
