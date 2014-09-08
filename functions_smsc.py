#!/usr/bin/python
# filename: functions_smsc.py
import numpy as np, glob , re , shutil , mmap, os, sys, math, logging, matplotlib.pyplot as plt
from copy import deepcopy # for function sort(). Probably find a better function!!
#for python-3 compatibility
from io import open 
#include [[Btree.py]]
import Btree as bt #class of binary tree

# Below are the conversion factors and fundamental constant
AMU2au=1822.88839                                          
Angs2Bohr=1/0.52917721092                                  
Hartree2GHz=6.579684e6                                     
Hartree2cm_1=219474.63 

def quantity(dim):
   F=np.zeros((2,dim,dim)) 
   CartCoord=np.zeros((2,3, dim/3))
   X=np.zeros((2,3,3))
   P=np.zeros((2,dim,dim))
   Energy=np.zeros(2)
   return F, CartCoord, X, P, Energy

def sort(f):
   """This function is a self-writen sort-function for floats-arrays 
   (increasing arguments by absolute value). Its input is the array whose elements 
   should not exceed 3e300 (otherwise the sorting will fails) and returns the number 
   of indices sorted by the size of respective elements. Hence, sorting an array A
   by the size of its elemens (largest first) can be done by 
   index=sort(A) 
   B=A[index]
   where B will be the sorted array.
   For not absolute values see resort()."""
   index=np.zeros(len(f), dtype=int)
   tmp=deepcopy(f) #for avoiding side effects
   for i in range(len(f)):
      index[i]=np.argmin(np.abs(tmp)) # there can be frequencies < 0 as well...
      tmp[index[i]]=5e+300 # this can be considered as smaller than all elements...
   return index

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
      a=np.sqrt(sum(temp_vec[j]*temp_vec[j] for j in range(len(temp_vec)) )) #calculate norm
      Y.append( temp_vec/a) #normalize all vectors
   return np.matrix(Y).T # undo transposition in the beginning

def replace(files, freq, L):
   """ This function creates a new file (determined by files, ending with 
   ".rep" and copies the log-file (files) into it, replacing the frequencies and 
   normal modes by those calculated by smallscript.
   The function is suited to test, whether these results coincide qualitatively with the gaussian's.

   ** PARAMETERS: **
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
	    logging.debug('frequencies not yet written to file:'+ repr(len(freq[s:].T))+ repr(freq[s:].T))
	    if len(freq[s:].T)> 2: # there are at least three more frequencies
	       out.write(re.sub(r'Frequencies -- [\d .-]+',
		     'Frequencies --'+'    '+repr(freq[s])+'          '\
		     +repr(freq[s+1])+'          '+repr(freq[s+2]), line))
	    elif len(freq[s:].T)== 2: # there are only two frequencies left
	       out.write(re.sub(r'Frequencies -- [\d .-]+',
		     'Frequencies --'+'     '+repr(freq[s])+'          '\
		     +repr(freq[s+1]), line))
	    elif len(freq[s:].T)== 1: # there is just one additional freq
	       out.write(re.sub(r'Frequencies -- [\d .-]+',
		     'Frequencies --'+'    '+repr(freq[s]), line))
	    s+=3
	 elif re.search(r'[ ]+\d+[ ]+\d+[ -]+\d.\d\d[ -]+\d.\d\d+[ \d.-]+', line) is not None:
	    if len(L[t][u:].T)> 2: # there are at least three more frequencies
	       out.write(re.sub(r'[\d .-]+', '    '+repr(t/3)+'    '+repr(s)+'   '+
		  '   '+str("%.2f" % L[t+0][u+0])+'  '+str("%.2f" % L[t+1][u+0])+' '+
			str("%.2f" % L[t+2][u+0])+
		  '   '+str("%.2f" % L[t+0][u+1])+'  '+str("%.2f" % L[t+1][u+1])+' '+
			str("%.2f" % L[t+2][u+1])+
		  '   '+str("%.2f" % L[t+0][u+2])+'  '+str("%.2f" % L[t+1][u+2])+' '+
			str("%.2f" % L[t+2][u+2]), line))
	    elif len(L[t][u:].T)== 2:
	       out.write(re.sub(r'[\d .-]+', '    '+repr(t/3)+'    '+repr(s)+'   '+
		  '   '+str("%.2f" % L[t+0][u+0])+'  '+str("%.2f" % L[t+1][u+0])+' '+
			str("%.2f" % L[t+2][u+0])+
		  '   '+str("%.2f" % L[t+0][u+1])+'  '+str("%.2f" % L[t+1][u+1])+' '+
			str("%.2f" % L[t+2][u+1]), line))
	    elif len(L[t][u:].T)== 1:
	       out.write(re.sub(r'[\d .-]+', '    '+repr(t/3)+'    '+repr(s)+'   '+
		  '   '+str("%.2f" % L[t+0][u+0])+'  '+str("%.2f" % L[t+1][u+0])+' '+
			str("%.2f" % L[t+2][u+0]), line))
	    t+=3 # 
	 else: 
	    out.write(re.sub('replace nothing','by nothing', line)) #just write line as it is
      out.close()

def gaussianfreq(ContntInfo, dim):
   """Extraction of the frequencies from g09-log file"""
   f=np.zeros((len(ContntInfo), dim-6))
   for i in range(len(ContntInfo)): 
      files=open(ContntInfo[i][0], "r") #open file and map it for better working
      mapedlog=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ) # i-th file containing freq calculations
      files.close
      f1=re.findall(r" Frequencies -- [\d .-]+", mapedlog, re.M)# 
      f2=[re.findall(r"[- ]\d+.\d+", f1[j]) for j in range(len(f1))]
      s=0
      for j in range(len(f2)):
	 f[i][s:s+len(f2[j])]=f2[j]
	 s+=len(f2[j])
   return f

# extracting the L-matrix from the log file for comparing with the one calculated above
def extractL(ContntInfo, dim):
   for i in range(len(ContntInfo)): 
      files=open(ContntInfo[i][0], "r") #open file and map it for better working
      mapedlog=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ) # i-th file containing freq calculations
      files.close
      f1=re.findall(r" Atom  AN [\n\d XYZ.-]+", mapedlog, re.M)
      mapedlog.close() 
      Ltemp=[re.findall(r"[- ]\d.\d\d", f1[j]) for j in range(len(f1))]
      assert len(Ltemp)>=0, 'PROGRAMM ERROR: file'+repr(ContntInfo[i][0])+ \
	       ' does not contain frequency informations.\
	       This should have been avoided above already.'
      if i==0:
	 L2=np.zeros((len(ContntInfo),dim,dim-6)) # similar to L
      for k in range(len(Ltemp)): # 
	 for j in range(len(Ltemp[k])): 
	    L2[i][j%3+(j/9)*3][(j/3)%3+3*k]=Ltemp[k][j]
      for k in range(len(L2[i][0])): #renormalization
	 L2[i][:].T[k]/=np.sqrt(np.sum( L2[i][:].T[k]*L2[i][:].T[k]))
   return L2

def TrafoCoord(F, P, Coord, dim):
   #assuming two files !! (more are impossible here...)
   Y=np.zeros((dim,dim)) 
   one=np.eye(dim)
   S=np.zeros((3,3))
   for i in range(3):
      for j in range(3):
	 S[i][j]=np.dot(Coord[0][i],Coord[1][j])
      S[i]/=np.linalg.norm(S[i])
   print 'Overlapmatrix of coordinates:\n', S
   for i in range(3):
      for j in range(3):
	 S[i][j]=round(S[i][j])
   print 'rounded overlapmatrix:\n', S
   Coord[1]=np.dot(S,Coord[1])
   for j in range(dim/3):
      Y[3*j:3*j+3].T[3*j:3*j+3]=S
   # Returning the read values
   logging.debug( 'new Coords\n'+ repr(S)+'\n'+repr(Y))
   F[1]=np.dot(np.dot(Y.T,F[1]),Y)
   #print P[1]
   #the following seems to have problems!!
   P[1]=one-np.dot(np.dot(Y.T,one-P[1]),Y) #for second file: D -> YD (Y.T: redirect of axes)
   #print P[1]
   return F, P, Coord

def Contentcheck(inputs): #look through it in detail! In the way it is, it is senseless!!
   """
   This function opens the given g09-log files and checks, whether frequency-calculations are done and if they converged. Otherwise these states are not taken into account.
   The function is valid for any number of files.
   
   Arguments:
   1.

   returns:
   A list containing the file-names and calculation-names (of g09) of successfull frequency-calculations
   """
   ContntInfo=[] #create array: name(char), Stoichiometry(char), root(int), freq(bool), 
   #for files in glob.glob('*.log'):
   for files in inputs:
      name= files
      log=open(name, "r") #open file and map it for better working
      mapedlog=mmap.mmap(log.fileno(), 0, prot=mmap.PROT_READ) # improve by finding a better size
      log.close
      stoi=re.findall(r"[\w]+" ,re.findall(r"-+\n [\w \d]+\n -+\n Charge =", mapedlog, re.I)[0])[0]
      freq=re.search(r"-+\n \#[\w\d\- \)\(\/=\n]+\
	    freq [\w\d \)\(\/=\-\n]+\n -+",
	    mapedlog, re.M) is not None #I hope, this is reasonable
      freq=re.search(r" freq ", mapedlog, re.I) is not None #I hope, this is reasonable
      if freq is True:
	 if (re.search(r" Error termination", mapedlog) is None) is True: #if successful calculations were done
	    ContntInfo.extend([[name, stoi]]) #write into array all information about files
	 else:
	    logging.error("Some error occured. In file "+\
		     name + ". It seems as if the calculation did not end correctly. Please check the log-file.")
   Contnt=ContntInfo
   return Contnt

def ReadLog(fileN):
   # Mapping the log file
   files=open(fileN, "r")
   log=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   files.close

   # Determine atomic masses in a.u. Note mass contains sqrt of mass!!!
   atmwgt=re.findall(r"AtmWgt= [\d .]+",log)
   mtemp=[]
   foonum=0
   for j in range(len(atmwgt)/2): # because atomic masses are printed twize in log-files...
      mtemp.append(re.findall(r'[\d.]+',atmwgt[j]))
   dim=0
   for j in range(len(mtemp)):
	 dim+=len(mtemp[j]) # dim will be sum over all elements of temp
   dim*=3
   mass=np.zeros(dim/3) # this is an integer since dim=3*N with N=atomicity
   for j in range(len(mtemp)):
      for k in range(len(mtemp[j])):
	 mass[k+foonum]=np.sqrt(float(mtemp[j][k])*AMU2au) #elements in m are sqrt(m_i) where m_i is the i-th atoms mass
      foonum+=len(mtemp[j])
   logging.debug("atmwgt: "+ repr(atmwgt)+'\nmtemp: '+ repr(mtemp))
   logging.info("Number of atoms: "+ repr(dim/3)+'\nDimensions of a problem: '+\
	 repr(dim)+'Sqrt of masses in a.u. as read from log file\n'+ repr(mass))

   # Reading Cartesian coordinates
   temp=[]
   temp=re.findall(r' Number     Number       Type             X           Y           Z[\n -.\d]+', log)
   tmp=re.findall(r'[ -][\d]+.[\d]+', temp[-1])
   assert len(tmp)==dim, 'Not all atoms were found! Something went wrong...'
   Coord=np.zeros((3, dim/3))
   MassCenter=np.zeros(3)
   for j in range(len(tmp)):
      Coord[j%3][j/3]=tmp[j]
   for j in range(3):
      MassCenter[j]=np.sum(Coord[j]*mass)
      MassCenter[j]/=np.sum(mass) #now it is cartesian center of mass
   logging.debug("Cartesian (Angstrom) coordinates before alignment to center of mass\n"+ repr(Coord.T)+\
	 "Center of mass coordinates (Angstrom)\n"+ repr(MassCenter))
   for j in range(3):#displacement of molecule into center of mass:
      Coord[j]-=MassCenter[j] # if commented we get rotational constants in agreement with Gaussian log
   logging.debug("Cartesian coords with respect to center of mass\n"+ repr(Coord.T))
   Coord*=Angs2Bohr
   logging.info("Cartesian coordinates (a.u.) in center of mass system\n"+repr(Coord.T))

   # Getting tensor of inertia, transforming to principlas axes
   moi=np.zeros((3,3))# this is Moment Of Inertia
   # print Coord[1]
   for j in range(3):
      for k in range(3):
	 if k is j:
	    moi[j][k]=np.sum(mass*mass*(Coord[0]*Coord[0]+\
		     Coord[1]*Coord[1]+Coord[2]*Coord[2]-\
		     Coord[j]*Coord[k]))
	 else:
	    moi[j][k]=np.sum(mass*mass*(Coord[j]*Coord[k]))
   logging.debug("Moments of intertia as read from log file\n"+repr(moi))
   diagI,X=np.linalg.eig(moi) # this can be shortened of course!
   index=sort(diagI)
   #X=np.matrix(X[index]) #sorting by eigenvalues
   X=np.matrix(X) #sorting by eigenvalues
   diagI=diagI[index]
   logging.debug("Moments of inertia (a.u.) in principle axes\n"+repr(diagI.T)+\
	 '\nRotational constants (GHz) in principle axes\n'+ repr(1/(2*diagI.T)*Hartree2GHz)+\
      "Rotation matrix\n"+repr(X))

   # Reading of Cartesian force constant matrix  
   f=re.findall(r"Force constants in Cartesian coordinates: [\n\d .+-D]+", log, re.M)
   f_str=str(f)[2:-2]
   lines=f_str.strip().split("\\n")
   F=np.zeros((dim,dim))
   n=0
   k=0
   for i in range(2,len(lines)):
      if i == dim+k-5*n+2: # is 'is' ok as well?
         k=i-1
         n+=1
         continue
      elements=lines[i].replace('D','e').split()
      for j in range(1,len(elements)):
         F[int(elements[0])-1][j-1+5*n]=float(elements[j])
         F[j-1+5*n][int(elements[0])-1]=float(elements[j])
   logging.debug('F matrix as read from log file\n'+ repr(F) +'\n0:9x0:9\n'+ repr(F[:9].T[:9].T))
   for i in range(0,dim):
      for j in range(0,dim):
         F[i][j]/= (mass[i/3]*mass[j/3]) 

   Etemp=re.findall(r'HF=-[\d.\n ]+', log, re.M)
   assert len(Etemp)>=1, 'Some error occured! The states energy can not be read.'
   if re.search(r'\n ', Etemp[-1]) is not None:
      Etemp[-1]=Etemp[-1].replace("\n ", "") 
   logging.info('temporary energy of state:'+repr(Etemp[-1]))
   E=-float(re.findall(r'[\d.]+', Etemp[-1])[0])# energy is negative (bound state)
   return dim, Coord, mass, X, F, E

def GetProjector(X, dim, m, Coord):
   D=np.zeros((dim,6))
   for k in range(3):# first three rows in D: The translational vectors
      for j in range(dim/3):
	 D[3*j+k][k]=m[j]
   for k in range(dim):# next three rows in D: The rotational vectors
      D[k][3:6]=(np.cross(np.dot(X,Coord)[:].T[k/3],X[:].T[k%3]))*m[k/3]

   logging.debug("Original translational and rotational displacement vectors"+repr(D[3:13].T))
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

def GetL(dim, mass, F, D):
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
   N=np.zeros(( len(F), len(F[0]), len(F[0]) )) 
   Lsorted=np.zeros(( len(F), len(F[0]), len(F[0])-6 )) 
   f=np.zeros(( len(F), len(F[0])-6 ))
   Ltemp=np.zeros(( len(F[0]), len(F[0])-6 ))
   ftemp=np.zeros(len(F[0]-6))
   for i in range(len(F)): #temporary for understanding only (this loop)
      ftemp,Ltemp=np.linalg.eig(F[i])
      index=sort(np.real(ftemp)) # ascending sorting f
      N[i]=Ltemp[index]
      foo=ftemp[index]
      logging.debug("Before projecting onto internal coords subspace:\n"+ 'Forces:\n'+\
	    repr(F[i])+'\nFrequencies (cm-1) \n'+ repr(np.sqrt(np.abs(ftemp[index]))*Hartree2cm_1)+\
            "\nL-matrix \n"+ repr(N[i]))

   for i in range(len(F)):
      ftemp,Ltemp=np.linalg.eig(np.dot(np.dot(D[i].T,F[i]),D[i]))
      assert np.any(ftemp<0) or np.imag(ftemp)!=0, 'Frequencies smaller than 0 occured. Please check the input-file!!'
      index=sort(np.real(ftemp)) # ascending sorting f
      f[i]=np.real(ftemp[index]).T[:].T[6:].T
      L[i]=np.real(Ltemp[index]).T[:].T[6:].T
      logging.debug("Frequencies (cm-1) \n"+ repr(np.sqrt(np.abs(ftemp[index]))*Hartree2cm_1))
      N[i]=np.real(Ltemp[index]) # or =Ltemp[index].T (see above)
      M=np.zeros((dim,dim))
      for j in range(0,dim):
         M[j,j]=1/mass[j/3]
      Lcart=np.dot(M,np.dot(D[i],np.real(Ltemp)))
      for j in range(0,dim):
         norm=np.sum(Lcart.T[j]*Lcart.T[j])
	 if np.abs(norm)>1e-12:
	    Lcart.T[j]/=np.sqrt(norm)
      Lsorted[i]=(Lcart.T[index].T)[:].T[6:].T
      logging.debug("Normalized Lcart\n"+ repr(Lcart)+"\nNormalized, sorted and truncated Lcart\n"+ repr(Lsorted[i]))

      for j in range(len(f[i])):
     	 f[i][j]=np.sign(f[i][j])*np.sqrt(np.abs(f[i][j]))
      logging.info("After projecting onto internal coords subspace\n"+"Frequencies (cm-1)\n"+\
	    repr(f[i]*Hartree2cm_1)+"L-matrix \n"+ repr(L[i]))
      
   #np.set_printoptions(precision=4, linewidth=122, suppress=True)
   return N, L, f, Lsorted

def GetLstab(dim, mass, F, G09f, P):
   """ Function that calculates the frequencies and normal modes from force constant matrix 
   with projection onto internal degrees of freedom and uses inverse iteration to get 
   numerically better results for normal modes. Therefore the g09-frequencies are used

   **argumets**
   1. The dimensions of force-constant matrix
   2. square-root of masses dim/3-dimensional array
   3. force-constant matrix
   4. dim-6-dimensional array of frequencies calculated by g09
   5. projection-matrix onto vibrations

   **return**
   1. matrix of vibr. normal modes (including global modes) 
   2. frequencies 
   3. massweighted L for comparison with the respective matrix from the g09 log-file
   """
   np.set_printoptions(suppress=True)
   np.set_printoptions(precision=4, linewidth=122)

   def Newton(F,L,Lambda):			   # inverse iteration to find better modes
      ones=np.eye(len(F))
      Tmp=np.zeros((len(F)+1,len(F)+1))
      tmp=np.zeros(len(F)+1)
      for i in range(len(L)):
	 for j in range(len(ones)):			   # there is no built-in for scalar*matrix
	    ones[j][j]*=Lambda[i]
	 A=np.matrix(F-ones)
      Tmp[:-1].T[:-1].T=A
      for i in range(len(L)):
	 Tmp[:].T[-1]=L[i]
	 Tmp[-1]=L[i].T
	 tmp[:-1].T=L[i]
	 tmp[-1]=Lambda[i]
	 resid=Tmp.dot(tmp)
	 while np.linalg.norm(resid, ord=2)<Lambda[i]*5e-5:
	    tmp[:-1]+=resid[:-1]
	    tmp[-1]-=resid[-1]#probably leaf this to fix frequencies...
      return L
   L=np.zeros(( len(F), len(F[0]), len(F[0])-6 )) 
   Lsorted=np.zeros(( len(F), len(F[0]), len(F[0])-6 )) 
   f=np.zeros(( len(F), len(F[0])-6 ))
   Ltemp=np.zeros(( len(F[0]), len(F[0])-6 ))
   ftemp=np.zeros(len(F[0]-6))
   for i in range(len(F)):
      ftemp,Ltemp=np.linalg.eig(np.dot(np.dot(P[i].T,F[i]),P[i]))
      #assert np.any(ftemp<0) or np.imag(ftemp)!=0, 'Frequencies smaller than 0 occured. Please check the input-file!!'
      #'repair' L-matrix using G09-frequencies
      index=sort(np.real(ftemp)) # ascending sorting f
      f[i]=np.real(ftemp[index]).T[:].T[6:].T
      L[i]=np.real(Ltemp[index]).T[:].T[6:].T # or =Ltemp[index].T (see above)
      logging.debug("Frequencies (cm-1) \n"+ repr(np.sqrt(np.abs(ftemp[index]))*Hartree2cm_1))
      L[i]=Newton(np.dot(np.dot(P[i].T,F[i]),P[i]), L[i].T, G09f[i]).T #'repair' L-matrix (truncated)
      M=np.zeros((dim,dim))
      for j in range(0,dim):
         M[j,j]=1/mass[j/3]
      Lcart=np.dot(M,np.dot(P[i],np.real(Ltemp)))
      for j in range(0,dim):
         norm=np.sum(Lcart.T[j]*Lcart.T[j])
	 if np.abs(norm)>1e-12:
	    Lcart.T[j]/=np.sqrt(norm)
      Lsorted[i]=(Lcart.T[index].T)[:].T[6:].T
      index=resort(f[i]) # only of temporary interest!! .......
      bar=L[i] #...............................................
      L[i]=bar[:].T[index].T #.................................
      bar=Lsorted[i] #.........................................
      Lsorted[i]=bar[:].T[index].T #...........................
      bar=f[i] #...............................................
      f[i]=bar[index] #........................................
      logging.debug("Normalized Lcart\n"+ repr(Lcart)+"\nNormalized, sorted and truncated Lcart\n"+ repr(Lsorted[i]))

      for j in range(len(f[i])):
     	 f[i][j]=np.sign(f[i][j])*np.sqrt(np.abs(f[i][j]))
      logging.info("After projecting onto internal coords subspace\n"+"Frequencies (cm-1)\n"+\
	    repr(f[i]*Hartree2cm_1)+"L-matrix \n"+ repr(L[i]))
   return L, f, Lsorted  #attention! The frequencies are NOT stabilized...

def Geometries(ContntInfo, problems):
   name=[]
   geometry=[] # will contain names of the files refering to different states
   coordName=[]
   coordNumb=[]
   geometry.append(ContntInfo[0][0])
   for i in range(len(ContntInfo)):  #loop over all states investigated
      files=open(ContntInfo[i][0], "r")#open file and map it for better working
      mapedlog=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
      files.close
      coords=re.findall(
	    r" ! Name  Definition[ ]+Value[]+Derivative Info.[ ]+![!,()RADEXcalutenyi\d /=\n .-]+",
	    mapedlog, re.I)
      if len(coords)==0:
	    print("in file "+ContntInfo[i][0]+' are no coordinate-informations')
	    continue
      mapedlog.close()
      coordName.append(re.findall(r" [RAD]{1}\([\d,]+\)", coords[-1]))
      coordNumb.append(re.findall(r" \d+.\d+", coords[-1]))
      check=np.zeros(len(coordName)-1)
      if i==0:
         for k in range(len(coordName[0])):
            coordNumb[0][k]=float(coordNumb[0][k])
      for j in range(len(coordName)-1):  #loop over all previously investigated files
         if len(coordName[j]) != len(coordName[i]): #if number of coordinates doesn't mach
	    check[j]+=1
	    continue
	 for k in range(len(coordName[j])):
	    if coordName[j][k] != coordName[i][k]: ##here an error occurs!!!!
	       check[j]+=1
	       continue
	    coordNumb[-1][k]=float(coordNumb[i][k])
	    if np.abs(coordNumb[-1][k])<0.1: #avoid dividing by 0
	       if np.abs(coordNumb[j][k]-coordNumb[-1][k])>3:
		  check[j]+=1
		  continue
	    elif np.abs(
		  (coordNumb[j][k]-coordNumb[-1][k])
		  /coordNumb[-1][k])>0.04: #difference more than 4%
	       check[j]+=1
	       continue
      if np.all(check > 0):
	 geometry.append(ContntInfo[i][0])
   return geometry

def Duschinsky(L, mass, dim, x):
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
   J=np.zeros((dim-6,dim-6))
   K=np.zeros(dim-6)
   I=np.zeros((len(L),dim))
   M=np.zeros((dim,dim))
   DeltaX=np.zeros((len(L),dim))

   for i in range(dim):
      M[i][i]=mass[i/3] #square root of masses
   for i in range(len(DeltaX[0])/3):
      DeltaX[0][3*i:3*i+3]=x[0].T[i]
   J=np.dot(L[0].T, np.linalg.pinv(L[1].T)) #always relaxation into ground state

   DeltaX=np.array(x[0]-x[1]).flatten('F')
   logging.debug('Flatted\n'+repr(DeltaX))
   #K=np.dot(L[1].T,DeltaX.T)
   K=np.dot(L[1].T.dot(M),DeltaX.T) #at the moment: mass-weighted

   np.set_printoptions(suppress=True)
   np.set_printoptions(precision=5, linewidth=138)
   #print('Duschinsky-rotation')
   #for i in range(len(J)):
      #for j in range(len(J[0])):
	 #print len(J)*2-2*i-1,'  ',2*j,'   ',J[i][j]
	 #print len(J)*2-2*i-1,'  ',2*j+1,'   ',J[i][j]
      #print
      #for j in range(len(J[0])):
	 #print len(J)*2-2*i-2,'  ',2*j,'   ',J[i][j]
	 #print len(J)*2-2*i-2,'  ',2*j+1,'   ',J[i][j]
      #print

   logging.info('Duschinsky rotation matrix, '+\
	 repr(np.linalg.norm(J[1]-np.eye(dim-6))/(dim-6))+ '  :\n'+ repr(J)+\
	 '  :\n'+ repr(J[:4].T[11:25])+\
	 '\nDuschinsky displacement vector:\n'+ repr(K))
   return J, K #J[0] is unity (if not everything goes wrong), second is of interenst...

def HuangR(K, f): #what is with different frequencies???
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
   #unif=np.zeros(len(K))
   logging.info('Delta Q:'+repr( K))
   unif=K*K*f[1]/(2)
   index=sort(unif)
   sortuni=unif[index]
   funi=f[1][index]
   if any(unif)<0:
      logging.warning('ATTENTION: some HR-factors are <0.\
	    In the following their absolute value is used.')
   uniHR=[]
   uniF=[]
   print 'HR-fact           freq'
   for j in range(len(unif)):
      #select all 'big' HR-factors 
      if sortuni[-j]>0.2: #can be changed to 'while' ore ommit sorting!?
	 uniHR.append(sortuni[-j])
	 uniF.append(funi[-j])
	 print  sortuni [-j],'  ',funi[-j]*Hartree2cm_1
   return uniHR, uniF

def calcspect(HR, freq, E, E0, N, M, T):
   """This is used to calculate the line spectrum assuming no mode mixing (shift only) and coinciding frequencies in both electronic states.

   **PARAMETERS:**
   HR:   Huang-Rhys factors
   n:    number of modes that are considered here (with biggest HR)
   freq: frequencies (have to be in the same order as HR
   E:    energy difference of energy surfaces
   N,M:  are the numbers of vibrational quanta can be in the modes
   All arguments are neccesary.

   **RETURNS:**
   nothing (output into /tmp/linspect)
   """

   def FCeqf( Deltag, M, N):
      """Calculate Franck-Condon factors under assumption of equal frequencies 
      for only one vibrational mode

      PARAMETERS:
      Deltag: HR-factor of respective state

      RETURNS:
      """
      exg=np.exp(-np.abs(Deltag)/2) #actually Deltag should be >0, but is not always due to negative frequencies
      faktNM=math.factorial(M)*math.factorial(N)
      FC=0
      for x in range(int(min(N,M))+1):
	 FC+=exg*math.pow(-1,N-x)*math.pow(np.abs(Deltag),(M+N-2*x)//2)/(math.factorial(M-x)*math.factorial(N-x))*\
	       math.sqrt(faktNM/(math.factorial(x)*math.factorial(x)))
      return FC
   
   def unifSpect(intens, freqs, E, FC00):
      """ Calculation of the line spectrum respecting only shift of minima (no Duschinsky rotation) 
      and assuming coinciding frequencies for initial and final state

      **PARAMETERS:**
      intens: matrix of intensities of transitions
      freqs:  matrix of respective energies

      **RETURNS**
      a 2-dimensional array with energies of transition (1. column) and their rate (2. column)
      """
      logging.debug('Spectrum\n'+ repr(N)+' '+repr(M)+'  '+repr(len(intens))+'  '+repr(len(intens[0])))
      J=len(intens[0])
      spect=np.zeros((2,len(intens)*len(intens[0])+1))
      spect[1][0]=FC00 #0->0 transition
      spect[0][0]=E
      #print "intensities"
      #print '0', spect[1][0], spect[0][0]
      for i in range(len(intens)):#make this easier: reshapeing of intens, freqs
	 for j in range(J):
	    spect[1][i*J+j+1]=intens[i][j]
	    spect[0][i*J+j+1]=freqs[i][j]
	    #print i*J+j+1, spect[1][i*J+j+1], spect[0][i*J+j+1]
      return spect
   
   n=len(HR) #=len(freq)
   FC=np.zeros((n,(M+1)*(N+1)-1)) 
   uency=np.zeros((n,(M+1)*(N+1)-1)) #freqUENCY
   #calculate 0->0 transition
   tmp=1
   for i in range(len(HR)): #this is actually not needed, right?
      tmp*=FCeqf(HR[i], 0,0)
   #scale whole spectrum
   FC00=tmp*tmp*1e4
   uency00=E*Hartree2cm_1 #zero-zero transition
   print "?   ?        0    0    ",E*Hartree2cm_1,"  ", FC00
   for a in range(n):
      temp=FCeqf(HR[a],0,0)
      for j in range(N+1):
	 for i in range(M+1):
	    if i==0 and j==0: 
	       continue #skip 0-0 transitions
	    tmp=FCeqf(HR[a], i, j)/temp
	    FC[a][j*(M+1)+i-1]=tmp*tmp*FC00*np.exp(-(E0+freq[a]*i)/T)
	    uency[a][j*(M+1)+i-1]=(E+freq[a]*(i-j))*Hartree2cm_1
	    print a,"  ",a ,'     ', i,'  ',j,'  ', (E+freq[a]*(i-j))*Hartree2cm_1,"  ",FC[a][j*M+i-1]
   FC00*=np.exp(-E0/T)
   return unifSpect(FC, uency,E*Hartree2cm_1, FC00)

def FCf(J, K, f, Energy, N):
   """Calculates the FC-factors for given Duschinsky-effect. No restriction to OPA
 
   
   *PARAMETERS:*
   J:      Duschisky-matrix
   K:      Displacement-Vector
   f:      frequency: two-dim array (freq_initial, freq_final)
   Energy: Energy-difference of minima
   N.      Max. number of excitation quanta state considered
     
   All parameters are obligatory.

   *RETURNS:*
   linespectrum 
   """
   def CalcI00(J, K, Gamma, Gammap, E):
      """This function calculates the overlap-integral for zero vibrations """
      invJ=np.linalg.inv(J)
      pref=math.pow(2,len(Gamma))*np.linalg.det(Gamma)
      TMP=J.dot(invJ.dot(Gammap).dot(J)+Gamma)
      pref/=np.linalg.det(TMP)

      #error!!!! This has to be positive but is not always!!
      pref=np.sqrt(pref)
      TMP=invJ.dot(Gammap).dot(J)+Gamma
      TMP=Gammap.dot(J).dot(np.linalg.inv(TMP)).dot(invJ)-np.eye(len(J))
      exp=np.exp(0.5*K.T.dot(TMP).dot(Gammap).dot(K))

      Tree=bt.Tree(2*len(K))
      Tree.fill(0)
      Zero=np.zeros(2*len(K))
      Tree.insert(Zero, [pref*exp, (E+sum(sum(Gammap-Gamma))/2)*Hartree2cm_1] ) #sum(sum()) due to matrix
      #I_00 transition-probability [[Btree.py#extract]]
      linspect.append(Tree.extract()) 
      return Tree
   
   def iterate(L1, L2, Energy, i, f, J, K):
      """ Calculates the Franck-Condon factors of an eletronic transition using the lower levels L1 and L2
   
      *PARAMETERS:*
      L1:     binary tree where i-2 quanta are excited (structure: [[Btree.py]]
      L2:     binary tree where i-1 quanta are excited
      Energy: Energy-difference between the states (minimal energy)
      i:      number of excitation-quanta
      f:      (2xN) frequencies of both states
      J:	   Duschisky-rotation matrix
      K:	   Displacement-vector

      *RETURNS:*
      L2:     input-parameter (needed for next iteration)
      L3:     new binary tree 
      """
   
      #quantities for the iterative spectrum-calculation
      Gamma=np.diag(f[0])              	# in atomic units. It is equivalent to 4pi^2/h f_i
      Gammap=np.diag(f[1])             	# for final state
      sqGamma=np.diag(np.sqrt(f[0]))   
      sqGammap=np.diag(np.sqrt(f[1]))  
      unity=np.eye(len(Gamma))
      invJ=np.linalg.inv(J)
   
      C=np.linalg.inv(invJ.dot(Gammap).dot(J)+Gamma) #C is only temporary matrix here
      A=np.dot(J,np.dot(C,invJ)) 
      A=2*np.dot(Gammap,np.dot(J,A))-unity
      b=unity-J.dot(C).dot(invJ).dot(Gammap)
      b=2*sqGammap.dot(unity-b).dot(K)
      E=4*sqGamma.dot(C).dot(invJ).dot(sqGammap)
      d=-2*sqGamma.dot(C).dot(invJ).dot(Gammap).dot(K)
      C=2*Gamma.dot(C)-unity 		#this is 'real' C-matrix
   
      #initialize new tree
      alpha=2*len(b)
      L3=bt.Tree(i)    	       		# initialize root-node
      L3.fill(alpha)         		# initialize tree
      States=states(alpha, i) 		# States are all possible

      def freq(E, Gamma, Gammap):
         """Calculates the frequency of respective transition including vibrational frequency

	 *PARAMETERS:*
	 E:	 energy-difference of states
	 Gamma:	 vector of vibrational frequencies in inital state (in atomic units)
	 Gammap: vector of vibrational frequencies in final state (in atomic units)

	 *RETURNS;*
	 frequency of respective transition
	 """
	 return (E+sum(Gammap-Gamma))*Hartree2cm_1 
   
      def FirstNonzero(n): 
	 """Find first non-zero elements in first and second half of array n """
	 ni=n[len(n)//2:] #interger division (python3-compatible)
	 nf=n[:len(n)//2]
	 m=len(ni)+1 #this means there is no excitation in this state
	 mp=len(nf)+1
	 for j in range(len(ni)):
	    if ni[j]>0:
	       m=j
	       break
	 for j in range(len(nf)):
	    if nf[j]>0:
	       mp=j
	       break
	 return m, mp

      for n in States: #for each possible state, described by n(vector)
	 m, mp= FirstNonzero(n)# index of first-non-zero element of (initial, final) state
	 # if the 'first' excited state is in initial state: need first iteration formula
	 I_nn=0
	 if m<=mp:
	    n_m=n[m]
	    ntemp=deepcopy(n)
	    ntemp[m]-=1 #n[m] is at least 1
	    Ps=L2.getState(ntemp)[0]
	    if not math.isnan(Ps) and abs(Ps)>1e-8:
	       I_nn=b[m]*Ps					# first term 
	    if ntemp[m]>0:
	       ntemp[m]-=1
	       Ps=L1.getState(ntemp)[0]
	       if not math.isnan(Ps) and abs(Ps)>1e-8:
		  I_nn+=np.sqrt(2*(n_m-1))*A[m][m]*Ps		# second term
	    for i in range(m+1, len(n)/2):
	       if n[i]>0:
		  ntemp=deepcopy(n)
		  ntemp[m]-=1
		  ntemp[i]-=1
		  Ps=L1.getState(ntemp)[0]
		  if not math.isnan(Ps) and abs(Ps)>1e-8:
		     I_nn+=np.sqrt(n[i]/2)*(A[m][i]+A[i][m])*Ps	# second term

	    for i in range(mp+len(n)//2, len(n)): 			# sum over respective final states
	       if mp>len(n)//2:					# that means: there are no excited vibrations
		  break
	       if n[i]>0:
		  ntemp=deepcopy(n)
		  ntemp[m]-=1
		  ntemp[i]-=1
		  Ps=L1.getState(ntemp)[0]
		  if not math.isnan(Ps) and abs(Ps)>1e-8:
		     I_nn+=np.sqrt(n[i]/2)*(E[i-len(n)//2][m])*Ps		# second term
	 #else: need the other iteration-formula
	 else: 
	    n_m=n[mp]
	    ntemp=deepcopy(n)
	    ntemp[mp]-=1
	    Ps=L2.getState(ntemp)[0]
	    if not math.isnan(Ps) and abs(Ps)>1e-8:
	       I_nn=d[mp]*Ps					# first term 
	    if ntemp[mp]>0:
	       ntemp[mp]-=1
	       Ps=L1.getState(ntemp)[0]
	       if not math.isnan(Ps) and abs(Ps)>1e-8:
		  I_nn+=np.sqrt(2*(n_m-1))*C[mp][mp]*Ps          	# second term
	    for i in range(mp+1, len(n)):
	       if n[i]>0:
		  ntemp=deepcopy(n)
		  ntemp[mp]-=1
		  ntemp[i]-=1
		  Ps=L1.getState(ntemp)[0]
		  if not math.isnan(Ps) and abs(Ps)>1e-8:
		     I_nn+=np.sqrt(n[i]/2)*(C[mp][i-len(n)//2]+    # second term
			      C[i-len(n)//2][mp])*Ps	
	    for i in range(m, len(n)): 				#sum over respective final states
	       if m>len(n)//2:					# that means: there are no excited vibrations
		  break
   	       if n[i]>0:
   		  ntemp=deepcopy(n)
   		  ntemp[mp]-=1
   		  ntemp[i]-=1
   		  Ps=L1.getState(ntemp)[0]
   		  if not math.isnan(Ps) and abs(Ps)>1e-8:
   		     I_nn+=np.sqrt(n[i]/2)*(E[mp][i-len(n)//2])*Ps 		# second term
     	 I_nn/=np.sqrt(2*n_m)
	 #threshold for insertion: saves memory, since int insead of float is used
	 if I_nn>1e-8:
	    try:
	       L3.insert(n, [I_nn, freq(Energy, f[0]*n[:len(n)//2], f[1]*n[len(n)//2:]) ])
	    except MemoryError: 
	       print('memory-error by inserting data. Finishing calculation.')
	       linspect=open('/tmp/linspect', "a")
	       linspect.writelines("%s\n" % item  for item in L2.extract())
	       linspect.close()
	       return 0,0
      return L2, L3

   def states(alpha, n): 
      """This function creates all possible states having a total number of n excitations in alpha different states
   
      *PARAMETERS:*
      alpha: number of degrees of freedom
   
      *RETURNS:*
      """

      def unlabeled_balls_in_labeled_boxes(balls, box_sizes): #needed for 'states'
	 """
	 These functions are part of python-package: 'combinatorics' 
	 (download from https://pypi.python.org/pypi/Combinatorics)
	 unlabeled_balls_in_labeled_boxes(balls, box_sizes): This function 
	 returns a generator that produces all distinct distributions of indistinguishable balls
	 among labeled boxes with specified box sizes (capacities). This is 
	 a generalization of the most common formulation of the problem, where each box is
	 sufficiently large to accommodate all of the balls, and is an important 
	 example of a class of combinatorics problems called 'weak composition' problems.
      
       	 OVERVIEW
      
       	 This function returns a generator that produces all distinct distributions of
	 indistinguishable balls among labeled boxes with specified box sizes
	 (capacities).  This is a generalization of the most common formulation of the
	 problem, where each box is sufficiently large to accommodate all of the
	 balls, and is an important example of a class of combinatorics problems
	 called 'weak composition' problems.
   
   
	 CONSTRUCTOR INPUTS
      
       	 n: the number of balls
	 
	 box_sizes: This argument is a list of length 1 or greater.  The length of
	 the list corresponds to the number of boxes.  `box_sizes[i]` is a positive
	 integer that specifies the maximum capacity of the ith box.  If
	 `box_sizes[i]` equals `n` (or greater), the ith box can accommodate all `n`
	 balls and thus effectively has unlimited capacity.
   
   
	 ACKNOWLEDGMENT
   
	 I'd like to thank Chris Rebert for helping me to convert my prototype
	 class-based code into a generator function.
	 """
	 def _unlabeled_balls_in_labeled_boxes(balls, box_sizes): #needed for 'unlabeled_balls_in_labeled_boxes' needed for 'states'
	    """
	    This recursive generator function was designed to be returned by
	    `unlabeled_balls_in_labeled_boxes`.
	    """
      
	    # If there are no balls, all boxes must be empty:
	    if not balls:
	       yield len(box_sizes) * (0,)
	 
	    elif len(box_sizes) == 1:
      
	       # If the single available box has sufficient capacity to store the balls,
	       # there is only one possible distribution, and we return it to the caller
	       # via `yield`.  Otherwise, the flow of control will pass to the end of the
	       # function, triggering a `StopIteration` exception.
	       if box_sizes[0] >= balls:
   		  yield (balls,)
      
	    else:
	       # Iterate over the number of balls in the first box (from the maximum
	       # possible down to zero), recursively invoking the generator to distribute
	       # the remaining balls among the remaining boxes.
	       for balls_in_first_box in xrange( min(balls, box_sizes[0]), -1, -1 ):
		  balls_in_other_boxes= balls - balls_in_first_box
		  for distribution_other in _unlabeled_balls_in_labeled_boxes(
		  balls_in_other_boxes, box_sizes[1:]):
		     yield (balls_in_first_box,) + distribution_other
	    # end three alternative blocks
   
	 if not isinstance(balls, int):
	       raise TypeError("balls must be a non-negative integer.")
	 if balls < 0:
	    raise ValueError("balls must be a non-negative integer.")
      
       	 if not isinstance(box_sizes,list):
	    raise ValueError("box_sizes must be a non-empty list.")
      
       	 capacity= 0
	 for size in box_sizes:
	    if not isinstance(size, int):
	       raise TypeError("box_sizes must contain only positive integers.")
	    if size < 1:
	       raise ValueError("box_sizes must contain only positive integers.")
	    capacity+= size
      
       	 if capacity < balls:
	    raise ValueError("The total capacity of the boxes is less than the "
	    "number of balls to be distributed.")
   
	 return _unlabeled_balls_in_labeled_boxes(balls, box_sizes)
	 # end def _unlabeled_balls_in_labeled_boxes(balls, box_sizes)
   
      #States=np.zeros((math.factorial(n+alpha-1)/(math.factorial(n)*math.factorial(alpha-1)),alpha))
      States2=[]
      a=np.ones(alpha).tolist()
      for i in range(len(a)):
	 a[i]=n*int(a[i]) #create the needed list
      i=0
      for distributions in unlabeled_balls_in_labeled_boxes(n,a):
	 #States[i]=np.matrix(distributions)
	 try:
	    States2.append(np.array(distributions, dtype=np.int8)) #save memory!
	 except MemoryError: 
	    #if the memory is not enough: don't add further states and use what is availibleA
	    print('Memory is full for state',n,'having only', len(States2),'objects in it. Use these states only.')
	    break #in principle: do more stuff here!
	 i+=1
      return States2


   Gamma=np.diag(f[0]) #in atomic units. It is equivalent to 4pi^2/h f_i
   Gammap=np.diag(f[1]) # for final state

   linspect=[]
   L2=CalcI00(J, K, Gamma, Gammap, Energy)
   #both trees can be expected to coincide for first state. 
   L1=L2 
   for i in range(1,N+1):
      L1, L2=iterate(L1, L2, Energy, i, f, J,K)
      #only by assert: MemoryError
      if L1==0 and L2==0:
	 linspect.append(L2.extract) #for this probably there is no more space
	 break #finish calculation
      linspect.append(L2.extract)
   return linspect #2-dimensional array

def outspect(gridpt, linspect, gamma, spectfile):
   """This function calculates the broadened spectrum given the line spectrum, frequency-rage and output-file whose name is first argument. 
   As basis-function a Lorentzian is assumed with a common width.
   
   **PARAMETERS:**
   spectfile: file, the result is written in (ascii-table). In addition a graph is created and shown on the fly. This graph is not saved.
   gridpt:    number of grid-points to be used for the calculation
   linspect:  line-spectrum list (frequency, intensity) 
   gamma:     broadening constant for the Lorentzians. It is the same for all peaks
   
   All parameters are obligatory."""
   out = open(spectfile, "w")
   #sort spectrum with respect to size of elements
   index=sort(linspect[1])
   linspect[1]=linspect[1][index]
   linspect[0]=linspect[0][index]
   #find transition with minimum intensity to be respected
   minint=0
   #for i in range(len(linspect[1])):
      #if linspect[1][i]>=0.001*linspect[1][-1]:
	 #minint=i
	 #break
   print('minimal and maximal intensities:\n', linspect[1][minint], linspect[1][-1])
   minfreq=linspect[0][np.argmin(linspect[0][minint:].T)] # min-freq   of fluorescence
   maxfreq=linspect[0][np.argmax(linspect[0][minint:].T)] # max freq
   print('maximal and minimal frequencies:\n', maxfreq, minfreq)
   minfreq-=300 #the range should be greater than the transition-frequencies
   maxfreq+=300 
   omega=np.linspace(minfreq,maxfreq,gridpt)
   spect=np.zeros(len(omega))
   #only those with high-enough intensities are respected
   for i in range(len(omega)): 
      intens=sum(linspect[1][j]/np.pi*gamma/((omega[i]-linspect[0][j])*(omega[i]-linspect[0][j])+ gamma*gamma)
	    for j in range(minint,len(linspect[1])))
      out.write(u" '{0}'  '{1}'\n".format(omega[i] ,intens))
      spect[i]=intens
   plt.plot(omega, spect)
   #plt.axis([2000,22000,0,5e-7])
   plt.title('Broadened spectrum of Ir-PS')
   plt.xlabel('Frequency [$cm^{-1}$]')
   plt.ylabel('Intensity (arb. units)')
   plt.show()
   out.close()

version=2.8
# End of functions_smsc.py
