#!/usr/bin/python
# filename: Btree.py
import numpy as np, math
#contstant is important for [[Btree.py#extract]]
DATATHRESHOLD=1e-6

class Tree:
   """the class Tree is a binary tree having a certain structure depending on 'alpha' 
   and (n+apha-1)!/(n!(alpha-1)!) data-points.
   For better structurisation, it has  4 different data-types:

   1. n: has two child-trees attached ('left' and *right*)
   2. l: has one child-tree ('left') and one data-point
   3. r: has one child-tree ('right') and one data-point
   4. u: has two data-points attached

   in addition, for the 0-th order tree there is a type '_' having only one data-point attached.
   """
   #all possible attributes; it saves memory especially for huge trees, (saves about 1/3)
   #see http://tech.oyster.com/save-ram-with-python-slots/
   __slots__=['data', 'data2', 'left','right','alpha'] #why is 'type' not neccesary?

   def __init__(self,alph):
      """ initializes the tree root of a (sub-) tree"""
      #self.alpha is not the same for different ns in the tree.
      # For the main-tree self.alpha =2*alpha-1 where alpha is the number of vibrational modes
      self.alpha=alph

   def fill(self, n):
      """ fills the Tree with the specific structure required for the level-fixed binary tree"""
      assert self.alpha!=0, 'There must be at least 1 vibrational mode!!'
      assert n>=0, 'The dimensionality of a tree can not be smaller 0'
      if n==0:
         self.data2=0 #this is extra-tree
         self.type='_'
      elif n==1 and self.alpha==1:
         self.data=np.array([0, 0],dtype=np.int8) #saves memory
         self.data2=np.array([0, 0],dtype=np.int8)
         self.type='u'
      elif self.alpha==1:
         self.data=np.array([0, 0],dtype=np.int8)
         self.right=Tree(self.alpha)
         self.right.fill(n-1)
         self.type='r'
      elif n==1:
         self.data=np.array([0, 0],dtype=np.int8)
         self.left=Tree(self.alpha-1)
         self.left.fill(1)
         self.type='l'
      else:
         self.right=Tree(self.alpha)
         self.right.fill(n-1)
         self.left=Tree(self.alpha-1)
         self.left.fill(n)
         self.type='n'

   def insert(self, N, FC):
      """ function that inserts some data into a specific placed denoted by the 2*alpha-dimensional array
      (note: 2*alpha-1 is 'self.alpha' of the main n)

      **PARAMETERS:** 
      N:  A 2*alpha-dimensional array 'attention: The size is never checked.  
      FC  The data to be filled into the respective n (2x1)-array

      **RETURNS:**
      nothing

      **NOTE:**
      Wrong-sized arrays will can lead to unexpected beaviour.
      Additionally the sum over all elemens has to coincide with 'n' used for 
      the function 'fill()'. This is not checked as well and will NOT lead 
      to asserts but will fill the element into wrong ns'
      """
      def summ(array):
         total=0
         for i,n in enumerate(N):
            total+=n
         return total

      n=M=0
      test=1
      #m=np.sum(N[i] for i in range(len(N)))
      m=summ(N)
      for i in range(self.alpha+1): #self.alpha==alpha since this is root-n
         if test==0:
            break
         M+=N[i]
         if self.type=='n':
            for j in range(0,int(N[i])):
               self=self.right
               n+=1
               if self.type=='l': 
                  if M==m:
                     #additional decreas of memory by factor 2
                     # float16 no more feasible due to iteration and strong rounding effects
                     self.data=np.array(FC, dtype=np.float32) 
                     test=0
                     break
            self=self.left
         elif self.type=='r': 
            for j in range(0,int(N[i])):
               self=self.right
               n+=1
               if self.type == 'u': #important for last elements
                  break #inner for-loop
            if M!=m and self.type=='u':
               self.data2=np.array(FC, dtype=np.float32)
            else:
               self.data=np.array(FC, dtype=np.float32)
            break
         elif self.type=='l':
            if M==m:
               self.data=np.array(FC, dtype=np.float32)
               break
            else:
               self=self.left
         else: #self.type=='u' or 0-th tree
            if N[i]>0:
               self.data=np.array(FC, dtype=np.float32)
            else:
               self.data2=np.array(FC, dtype=np.float32)
            break

   def extract(self): #extract all elements 
      """
      **extract**
      This function is for ordered extraction of all elements in the tree and additionally can be used 
      instead of print (by 'print instance.extract()').
      the returned statement is a vector containing all values created by the inner function

      *RETURN:*
      two-dimensional array containing intensities([0]) and frequencies([1])
      """
      result=[]

      def extra(self,result):
         """creates the vector to be returned in a recursive way 
            **PARAMETER**
            self:       object
            result:     2-dimensional array containing data already read from the tree, data will be added to it using 
                        side effects
            """
         if self.type=='n':
            extra(self.right,result)
            extra(self.left,result)
         elif self.type=='l':
            if self.data[0]>DATATHRESHOLD and not math.isnan(self.data[0]):
               result.append(self.data)
            extra(self.left,result)
         elif self.type=='r':
            extra(self.right,result)
            if self.data[0]>DATATHRESHOLD and not math.isnan(self.data[0]):
               result.append(self.data)
         elif self.type=='u':
            if self.data[0]>DATATHRESHOLD and not math.isnan(self.data[0]):
               result.append(self.data)
            if self.data2[0]>DATATHRESHOLD and not math.isnan(self.data2[0]):
               result.append(self.data2)
         else:
            if self.data2[0]>DATATHRESHOLD and not math.isnan(self.data2[0]):
               result.append(self.data2)

      extra(self,result)
      #return np.matrix(result) #
      if len(result)==0: #if no intensity is higher than DATATHRESHOLD
         return [[0, 0]]
      return result

   def getState(self, N): 
      """This function returns a particular state denoted by the array N (comprising thi initial and final state).

      **Input-Arguments:**
      Vector of length 2*alpha

      **returns:**
      the Franck-Condon factor for this state
      """
      def summ(array):
         total=0
         for i,n in enumerate(N):
            total+=n
         return total

      if any(N)<0: 
         return 0 ##from iteration it will happen
      n=M=0
      #m=np.sum([N[i] for i in range(len(N))])
      m=summ(N)
      for i in range(self.alpha+1): #self.alpha==alpha since this is root-n
         M+=N[i]
         if self.type=='n':
            for j in range(0,int(N[i])):
               self=self.right
               n+=1
               if self.type=='l': #exception for 'diagonal' elements
                  if M==m:
                     return self.data
            self=self.left
         elif self.type=='r': #allways has success
            for j in range(0,int(N[i])):
               self=self.right
               n+=1
               if self.type == 'u': #important for last elements
                  break #inner for-loop
            if M!=m and self.type=='u':
               return self.data2
            else:
               return self.data
         elif self.type=='l':
            if M==m:
               return self.data
            else:
               self=self.left
         else: #self.type=='u':
            if N[i]>0:
               return self.data
            else:
               return self.data2

   def size(self):
      """ function that returns the size of the tree (number of data-points). 
      It is needed for temporary needs and debugging only since the size can be calculated
      by the input-argumens of '__init__' and 'fill'.
      """
      if self.type=='u':
         return 2
      if self.type=='r':
         return self.right.size() + 1
      if self.type=='l':
         return self.left.size() + 1
      if self.type=='n':
         return self.left.size() + self.right.size()

version=1.3
# End of Btree.py
