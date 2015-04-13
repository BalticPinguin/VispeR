#!/usr/bin/python
# filename: readLogs.py
import numpy as np
import re, mmap
AMU2au=1822.88839                                          
Hartree2cm_1=219474.63 

def ReadG09(logging, fileN):
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
      f_str=str([f[-1]])#[2:-2]
      lines=f_str.strip().split("\\n")
      F=np.zeros((dim,dim))
      n=0
      k=0
      line=0
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

def ReadG092(logging, final):
   """ This function reads the required quantities from the gaussian-files

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

def getGaussianf(final, dim):
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

def getGaussianL(final, dim):
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

def replace(logging, files, freq, L): # this is G09-method
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

def ReadGAMESS(logging, fileN):
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
      print lines[i+s]
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
   print F

   if logging[0]<1:
      logging[1].write('F matrix as read from log file\n{0} \n'.format(F))
##########  is this correct?? Probably F is already mass-weighted?
   for i in range(0,dim):
      for j in range(0,dim):
         F[i][j]/= (mass[i//3]*mass[j//3]) 

   #get energy of the state
#         TOTAL ENERGY =    
   Etemp=re.findall(r'(?<=TOTAL ENERGY =)[\-\d. ]+', log, re.M)
   assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
   if logging[0]<=1:
      logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
   E=float(Etemp[-1]) #energy is negative (bound state)
   return dim, Coord, mass, F, E

def ReadGAMESS2(logging, final):
   """ This function reads the required quantities from the gaussian-files

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

def ReadNWChem(logging, fileN):
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
     #atom    #        X              Y              Z            mass
   WeightCoord=re.findall(r"(?<=atom    #        X              Y              Z            mass\n)[\d\w\-\+ \.\n]+", log)
   mtemp=re.findall(r"[ -]\d\.[\dD\.\+\-]+(?=\n)", WeightCoord[-1]) #only weights remain
   dim=3*len(mtemp)

   mass=np.zeros(dim/3) # this is an integer since dim=3*N with N=atomicity
   for j in range(len(mtemp)):
      mass[j]=np.sqrt(float(mtemp[j].replace('D','e'))*AMU2au) #mass[j] is square-root of masses
   assert not np.any(mass==0) , "some atomic masses are zero. Please check the input-file! {0}".format(mass)
   if logging[0]<2:
      logging[1].write("Number of atoms: {0}\nNumber of vibrational "
                        "modes: {1} \n Sqrt of masses in a.u. as read from log file\n{2}\n".format(dim/3,dim,mass))
   
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
         F[int(elements[0])-1][j-1+10*n]=float(elements[j])
         F[j-1+10*n][int(elements[0])-1]=float(elements[j])
   if logging[0]<1:
      logging[1].write('F matrix as read from log file\n{0} \n'.format(F))

   #get energy of the state
   Etemp=re.findall(r'(?<=Total energy =)[\-\d. ]+', log, re.M)
   assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
   if logging[0]<=1:
      logging[1].write('temporary energy of state: {0}\n'.format(Etemp[-1]))
   E=float(Etemp[-1]) #energy is negative (bound state)
   return dim, Coord, mass, F, E

def ReadNWChem2(logging, final):
   """ This function reads the required quantities from the gaussian-files

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

version=0.2
# End of readLogs.py
