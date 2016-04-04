#!/usr/bin/python2
# filename: normal_modes.py
import math, sys
import numpy as np
#from scipy.linalg.lapack import dsyev

# CHANGELOG
# =========
#in version 0.1.0:  
#   1) initialised class
#
class NormalMode():
   """Frament the functions below as well.
      This class handles objects related to normal mode analysis.
      Hence, the calculation of normal modes and frequencies is performed
      here (GetL) and their shift-vector/couplings are calculated (Duschinsky).
      The other functions (SortL and GetProjector) are functions to manipulate the 
      data: SortL assignes the normal modes of final state to fit the initial states
      coordinates best (which might become problematic when almost degenerate frequencies
      occur) and GetProjector projects out the rotations and vibrations from the force
      constant matrix. 
      This is also mostly for numerical purposes and turned out to be not that
      easy in practice.
   """
   Hartree2cm_1=219474.63 
   L=[]
   K=[]
   J=[]
   grad=[]
   F=[]
   Grad=[]
   mass=[]
   Coord=[]
   dim=0

   def __init__(self, parent):
      """Here, I just need to initialise the respective data.
         The rest will be done at the respective places in Spect.
      """
      self.log=parent.log
      self.F=parent.F
      self.mass=parent.mass
      self.Coord=parent.CartCoord
      self.dim=len(self.F[0])
      self.Grad=parent.Grad
      self.parent=parent

   def GetL(self):
      """Function that calculates the frequencies and normal modes from force constant matrix.It
         mainly calls a hermitian eigenvalue-solver and computes the frequencies as 
         srqrt(eigenvalue) and the normal modes by cutting the translation/rotation off.
         A projection method once has been in use but turned out to be pointless.
   
         This function is a member of Spect but syntactically not bound to it at the moment.
         Hence, it has access to member-variables.
      
         **PARAMETERS** 
         self -> it is member of class
         F    -> two matrices, F[0] and F[1] ; the respective nuclear energy Hessians.
   
         **RETURN**
         f  -> frequencies of initial (f[0]) and final (f[1]) state. The dimesion is (2, dim-6)
         L  -> unitary matrices (L[0], L[1]) that diagonalise M*F (massweighted Hessian). Its columns are normal modes 
               in Cartesian Coordinats The dimension is (2, dim, dim-6)
         Lmass -> L*M where M_ij=1/sqrt(m_i m_j). This is mass-weighted L and will be used later for most systems.
        
      """
      # Defining arrays
      lenF=len(self.F[0])
      L=np.zeros(( 2, lenF, lenF-6 )) 
      Lmass=np.zeros(( 2, lenF, lenF-6 ))
      f=np.zeros(( 2, lenF-6 ))

      # do the following for both states:
      for i in [0,1]:
         # solve the eigenvalue-equation for F:

         #REMOVE ROTATIONS AND TRANSLATIONS FROM HESSIAN.
         project=True
         #don't do this, if the force-constant matrix is just a copy
         # of the initial state to keep the results consistent.
         if self.parent.sameF==True:
            if i==1:
               project=False
         #project out the rotations and vibrations and not just throw away the smallest 6 eigen values
         # and respective eigen modes.
         if project:
            D=self.GetProjector(i)
            self.F[i]=np.dot(np.dot(D.T,self.F[i]),D)
            #Why can I not construct D such that I don't need to throw away anything?
         else:
            #A copy of F was used before, so I need to recopy it to have the projected
            # version here as well.
            self.F[i]=self.F[0]

         ftemp,Ltemp=np.linalg.eigh(self.F[i])
         #ftemp,Ltemp,info=dsyev(self.F[i])

         #sort the results with increasing frequency (to make sure it really is)
         # use absolute values for this.
         index=np.argsort(np.abs(ftemp),kind='heapsort') # ascending sorting f
         # and then cut off the 6 smallest values: rotations and vibrations.
         f[i]=np.real(ftemp[index]).T[:].T[6:].T
         L[i]=(Ltemp.T[index].T)[:].T[6:].T
         #END REMOVING ROTATIONS AND TRANSLATIONS.
   
         #the frequencies are square root of the eigen values of F
         # Here, we need to take care for values <0 due to bad geometry
         for j in range(len(f[i])):
            f[i][j]=np.sign(f[i][j])*np.sqrt(np.abs(f[i][j]))

         #through a warning if imaginary frequencies occured and take their pos. values for that.
         if np.any(f[i]<0):
            self.log.write('WARNING: imaginary frequencies occured. The absolute'
                    ' values are used in the following.\n{0}\n'.format(f[i]))
            f[i]=np.abs(f[i])

         #calculate the mass matrix, for Lmass
         M=np.eye(self.dim)
         for j in range(self.dim):
            M[j,j]/=self.mass[j//3]
         Lmass[i]=np.dot(M,L[i])
         #if log level=0 or 1:
         if i==0:
            self.log.write("Initial states ",2)
         else:
            self.log.write("Final states ",2)
         self.log.write("frequencies (cm-1)\n",2)
         if self.log.level<2:
            self.log.printVec(f[i]*self.Hartree2cm_1)
      
      # to account for the effect that the frequencies independently change between the states 
      #  and hence may change their order, the final states L and f are resorted via the
      #  largest overlap. When strong changes occur, there is no fair method. This one tries
      #  to be as fair as possible.
      f[1],L[1]=self.SortL(np.dot(np.linalg.pinv(L[0]),L[1]),L[1],f[1])

      #recalculate Lmass for final state.
      Lmass[1]=np.dot(M, L[1]) # recalculate Lmass!
      self.Lmassw=Lmass
      self.f=f
      self.L=L
      #finally, print the L-matrices as well:
      if self.log.level<2:
         self.log.write("initial state L-matrix \n")
         self.log.printMat(Lmass[0])
         self.log.write("final state L-matrix \n")
         self.log.printMat(Lmass[1])
   
   def __gs(self, A):
      """This function does row-wise Gram-Schmidt orthonormalization of matrices. 
         code for Gram-Schmidt adapted from iizukak, see https://gist.github.com/iizukak/1287876
      """
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

   def GetProjector(self, i):
      """ This function calculates a projection-matrix that is used to project the mass-weighted
         Hessian onto the space of vibrations. Therefore, we first construct a projector D onto
         translations and rotations and than apply 1-D onto F.
      """
      # Getting tensor of inertia, transforming to principlas axes
      moi=np.zeros((3,3))# this is Moment Of Inertia
      # print Coord[1]
      for j in [0,1,2]:
         for k in [0,1,2]:
            if k == j:
               moi[j][j]=np.sum(self.mass*self.mass*(self.Coord[i][0]*self.Coord[i][0]+\
                        self.Coord[i][1]*self.Coord[i][1]+self.Coord[i][2]*self.Coord[i][2]-\
                        self.Coord[i][j]*self.Coord[i][k]))
            else:
               moi[j][k]=np.sum(self.mass*self.mass*(self.Coord[i][j]*self.Coord[i][k]))
      diagI,Moi_trafo=np.linalg.eig(moi) # this can be shortened of course!
      index=np.argsort(diagI,kind='heapsort')
      #X=np.matrix(X[index]) #sorting by eigenvalues
      Moi_trafo=np.matrix(Moi_trafo) #notsorting by eigenvalues
   
      #now, construct the projector onto frame of rotation and translation using Sayvetz conditions
      D=np.zeros((self.dim,6))
      for k in [0,1,2]:# first three rows in D: The translational vectors
         for j in range(self.dim//3):
            #translations in mass-weighted coordinates
            D[3*j+k][k]=self.mass[j]
      for k in range(self.dim):# next three rows in D: The rotational vectors
         #rotations in mass-weighted coordinates
         D[k][3:6]=(np.cross(np.dot(Moi_trafo,self.Coord[i])[:].T[k//3],Moi_trafo[:].T[k%3]))*self.mass[k//3]
      D_orthog=self.__gs(np.array(D)) #orhogonalize it
      ones=np.identity(self.dim)
      one_P=ones-np.dot(D_orthog,D_orthog.T)
      prob_vec=(D_orthog.T[1]+D_orthog.T[4]+D_orthog.T[0]+D_orthog.T[5]).T #what is this actually??
      assert not np.any(np.abs(prob_vec-np.dot(np.dot(D_orthog,D_orthog.T),prob_vec))>0.00001), \
               'Translations and rotations are affected by projection operator.'+\
               repr(np.abs(prob_vec-np.dot(np.dot(D_orthog,D_orthog.T),prob_vec)))
      assert not  np.any(np.abs(np.dot(one_P,prob_vec))>0.00001), \
               "Projecting out translations and rotations from probe vector"
      return one_P
   
   def SortL(self, J,L,f):
      """This functions resorts the normal modes (L, f) such that the Duschinsky-Rotation
         matrix J becomes most close to unity (as possible just by sorting).
         In many cases, here chosing max(J[i]) does not help since there will be rows/columns occur
         more often. 
         Since I don't know any closed theory for calculating this, it is done by cosidering all possible cases.
      """
      
      #initialize the matrix that resorts the states:
      resort=np.zeros(np.shape(J))
      #FIRST, DO SOME GUESS HOW THE MATRIX COULD LOOK LIKE
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

      # NOW LOOK FOR ERRORS IN THE GUESS
      for i in range(len(J)):
         if sum(resort[i])==1:
            #this is the normal case: the order of
            # states did not change.
            continue
         elif sum(resort[i])==0:
            Nos.append(i)
         else:
            index=np.where(resort[i]==1)
            x=np.argmax(np.abs(J[index,i]))
            index=np.delete(index,x)
            resort[i][index]=0 #only x remains
            freePlaces=np.append(freePlaces,index)
      # By construction, this should always be true!
      assert len(Nos)==len(freePlaces), "dodododo!"
      freePlaces=np.array(freePlaces,dtype=int)
      
      #FIXING THE ERRORS NOW:
      #fill remaining lines. Therefore, set that element to one
      # whose value is largest under all remaining ones.
      # This method is not fair since the first has most choise but should
      # be fair enough in most cases
      for i in range(len(Nos)):
            x=np.argmax(np.abs(J[freePlaces,Nos[i]]))
            resort[Nos[i],freePlaces[x]]=1 #only x remains
            freePlaces=np.delete(freePlaces,x) # does it work the way I want it to work?
      assert len(freePlaces)==0, "the matrix is not readily processed."
      #FIX DONE.
      
      #since resort is a permutation matrix, it is unitary. Using this:
      return np.dot(f,resort), np.dot(L,resort)
      #  END OF SortL
   
   def Duschinsky(self):
      """This function calculates the shift between two electronic states 
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
   
      for i in range(self.dim):
         M[i][i]=self.mass[i//3] #square root of inverse masses
      #J=np.dot(L[0].T, L[1])  # for Lsorted
      self.J=np.dot(np.linalg.pinv(self.Lmassw[0]),self.Lmassw[1]) # for Lmassw
   
      #print "J\n", J
      if any(self.Grad[i]>0 for i in range(len(self.Grad))):
         self.K=np.dot(self.Grad.T,self.Lmassw[0])
         #self.K=(np.linalg.pinv(self.Lmassw[0]).dot(self.Grad)).T  # w p Lmassw
         # scale consistently: Now it is really the shift in terms of normal modes
         self.K/=self.f[0]*self.f[0]#*np.sqrt(2)  #
         #self.K=self.K[0] # it is matrix and needs to be a vector...
         #FIXME: this seems to be inconsistent! may be matrix or vector...

      else:
         DeltaX=np.array(self.Coord[1]-self.Coord[0]).flatten('F')  # need initial - final here.
         if self.log.level <1:
            self.log.write('changes of Cartesian coordinates:\n')
            self.log.printVec(DeltaX)
         self.K=np.dot(np.linalg.pinv(self.Lmassw[0]), DeltaX)  # w p Lmassw
      
      if self.log.level<2:
         # print the Duschinsky matrix in a nice format
         self.log.write('Duschinsky rotation matrix:\n')
         self.log.printMat(self.J)
         self.log.write('\nDuschinsky displacement vector:\n')
         self.log.printVec(self.K)

version='0.1.0'
#End of normal_modes.py
