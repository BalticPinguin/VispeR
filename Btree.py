#!/usr/bin/python
# filename: Btree.py
import numpy as np, math
#contstant is important for [[Btree.py#extract]]
DATATHRESHOLD=1e-6

class Tree:
   """the class Tree is a binary tree having a certain structure depending on 'alpha' 
   and (n+apha-1)!/(n!(alpha-1)!) data-points.
   For better structurisation, it has  4 different data-types:

   1. node: has two child-trees attached ('left' and *right*)
   2. leaf: has one child-tree ('left') and one data-point
   2. reaf: has one child-tree ('right') and one data-point
   2. lleaf: has two data-points attached

   in addition, for the 0-th order tree there is a type _ having only one data-point attached.
   """
   #all possible attributes; it saves memory especially for huge trees,
   #see http://tech.oyster.com/save-ram-with-python-slots/
   __slots__=['data', 'data2', 'left','right','alpha'] 

   def __init__(self,alph): #should it be __new__?
      """ initializes the tree root of a (sub-) tree"""
      #self.alpha is not the same for different nodes in the tree.
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
	 self.data=int([0, 0])
	 self.data2=int([0, 0])
	 self.type='lleaf'
      elif self.alpha==1:
  	 self.data=int([0, 0])
  	 self.right=Tree(self.alpha)
  	 self.right.fill(n-1)
	 self.type='reaf'
      elif n==1:
  	 self.data=int([0, 0])
  	 self.left=Tree(self.alpha-1)
  	 self.left.fill(1)
	 self.type='leaf'
      else:
  	 self.right=Tree(self.alpha)
  	 self.right.fill(n-1)
	 self.left=Tree(self.alpha-1)
	 self.left.fill(n)
	 self.type='node'

   def insert(self, N, FC): #data=[........] alpha-dim array...
      """ function that inserts some data into a specific placed denoted by the 2*alpha-dimensional array
      (note: 2*alpha-1 is 'self.alpha' of the main node)
      **Argumens**
      1. A 2*alpha-dimensional array 'attention: The size is never checked. 
      Wrong-sized arrays will can lead to unexpected beaviour.
      Additionally the sum over all elemens has to coincide with 'n' used for 
      the function 'fill()'. This is not checked as well and will NOT lead 
      to asserts but will fill the element into wrong nodes'
      2. The data to be filled into the respective node
      """
      n=M=0
      test=1
      m=np.sum(N[i] for i in range(len(N)))
      for i in range(self.alpha+1): #self.alpha==alpha since this is root-node
	 if test==0:
	    break
	 M+=N[i]
	 if self.type=='node':
	    for j in range(0,int(N[i])):
	       self=self.right
	       n+=1
	       if self.type=='leaf': #exception for 'diagonal' elements
		  if M==m:
		     self.data=FC
		     test=0
		     break
	    self=self.left
	 elif self.type=='reaf': #allways has success
	    for j in range(0,int(N[i])):
	       self=self.right
	       n+=1
	       if self.type == 'lleaf': #important for last elements
		  break #inner for-loop
	    if M!=m and self.type=='lleaf':
  	       self.data2=FC
	    else:
	       self.data=FC
	    break
	 elif self.type=='leaf':
	    if M==m:
	       self.data=FC
	       break
	    else:
	       self=self.left
	 else: #self.type=='lleaf' or 0-th tree
	    if N[i]>0:
	       self.data=FC
	    else:
	       self.data2=FC
	    break

   def extract(self): #extract all elements 
      """
      ===extract===
      This function is for ordered extraction of all elements in the tree and additionally can be used 
      instead of print (by 'print instance.extract()').
      the returned statement is a vector containing all values created by the inner function

      *RETURN:*
      two-dimensional array containing intensities([0]) and frequencies([1])
      """
      result=[]

      def extra(self,result):
   	 """creates the vector to be returned in a recursive way """
	 if self.type=='node':
	    extra(self.right,result)
	    extra(self.left,result)
	 elif self.type=='leaf':
	    if self.data[0]>DATATHRESHOLD and not math.isnan(self.data[0]):
   	       result.append(self.data)
	    extra(self.left,result)
	 elif self.type=='reaf':
	    extra(self.right,result)
	    if self.data[0]>DATATHRESHOLD and not math.isnan(self.data[0]):
   	       result.append(self.data)
	 elif self.type=='lleaf':
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
	 return 0
      return result

   def getState(self, N): 
      """This function returns a particular state denoted by the array N (comprising thi initial and final state).

      **Input-Arguments:**
      Vector of length 2*alpha

      **returns:**
      the Franck-Condon factor for this state
      """
      if any(N)<0: 
	 return 0 ##from iteration it will happen
      n=M=0
      m=np.sum(N[i] for i in range(len(N)))
      for i in range(self.alpha+1): #self.alpha==alpha since this is root-node
	 M+=N[i]
	 if self.type=='node':
	    for j in range(0,int(N[i])):
	       self=self.right
	       n+=1
	       if self.type=='leaf': #exception for 'diagonal' elements
		  if M==m:
		     print  self.data
		     return self.data
	    self=self.left
	 elif self.type=='reaf': #allways has success
	    for j in range(0,int(N[i])):
	       self=self.right
	       n+=1
	       if self.type == 'lleaf': #important for last elements
		  break #inner for-loop
	    if M!=m and self.type=='lleaf':
  	       return self.data2
	    else:
	       return self.data
	    break
	 elif self.type=='leaf':
	    if M==m:
	       return self.data
	    else:
	       self=self.left
	 else: #self.type=='lleaf':
	    if N[i]>0:
	       return self.data
	    else:
	       return self.data2

   def size(self):
      """ function that returns the size of the tree (number of data-points). 
      It is needed for temporary needs and debugging only since the size can be calculated
      by the input-argumens of '__init__' and 'fill'.
      """
      if self.type=='lleaf':
	 return 2
      if self.type=='reaf':
	 return self.right.size() + 1
      if self.type=='leaf':
	 return self.left.size() + 1
      if self.type=='node':
	 return self.left.size() + self.right.size()

   def compress(float32):
      """Code from davidejones (http://davidejones.com/blog/1413-python-precision-floating-point/)
      does this help reducing the memory-requirement of my function?
      """ 
      F16_EXPONENT_BITS = 0x1F
      F16_EXPONENT_SHIFT = 10
      F16_EXPONENT_BIAS = 15
      F16_MANTISSA_BITS = 0x3ff
      F16_MANTISSA_SHIFT =  (23 - F16_EXPONENT_SHIFT)
      F16_MAX_EXPONENT =  (F16_EXPONENT_BITS << F16_EXPONENT_SHIFT)
      a = struct.pack('>f',float32)
      b = binascii.hexlify(a)
      f32 = int(b,16)
      f16 = 0
      sign = (f32 >> 16) & 0x8000
      exponent = ((f32 >> 23) & 0xff) - 127
      mantissa = f32 & 0x007fffff
      if exponent == 128:
	 f16 = sign | F16_MAX_EXPONENT
      if mantissa:
	 f16 |= (mantissa & F16_MANTISSA_BITS)
      elif exponent > 15:
	 f16 = sign | F16_MAX_EXPONENT
      elif exponent > -15:
	 exponent += F16_EXPONENT_BIAS
	 mantissa >>= F16_MANTISSA_SHIFT
	 f16 = sign | exponent << F16_EXPONENT_SHIFT | mantissa
      else:
	 f16 = sign
      return f16

version=1.3
# End of Btree.py
