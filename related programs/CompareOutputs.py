#!/usr/bin/python
#include [[readLogs.py]]
import readLogs as rl
#include further dicts
import sys, re, mmap, numpy as np, os
Hartree2cm_1=219474.63 

def main(argv=None):
   numfiles=len(argv) #check, how many files are availible
   #ask, what should be done:
   #1:  Get Hessian of all and compare them
   #2:  Get Frequencies and eigenvectors and compare them
   #3: both
   print "Please type the actions to be done:\n"\
         " 1:  Get Hessian of all and compare them\n 2:  Get Frequencies"\
         " and eigenvectors and compare them\n 3: both\n"
   #get user-input:
   todo = raw_input("your choise: ")
   todo=int(todo)

   if todo==1 or todo==3:
      F=[]
      coord=[]
      E=[]
      print "get Hessian"
      for i in range(numfiles):
         print i,":"
         f,c,e=ReadF(argv[i])
         F.append(f)
         coord.append(c)
         E.append(e)
      print "have all",numfiles," Hessians"
      if numfiles==1:
         print "nothing to compare"
      else:
         print "                       norm(F[i]-F[j])     E[i]-E[j]   norm(coord[i]-coord[j]) i   j "
         for i in range(numfiles):
            for j in range(i+1,numfiles):
               #print "norm of difference:", np.linalg.norm(F[i]-F[j]), E[i]-E[j], np.linalg.norm(coord[i]-coord[j]), i,j
               print "norm of difference:", np.linalg.norm(np.matrix(F[i])-np.matrix(F[j])), E[i]-E[j],\
                     np.linalg.norm(np.matrix(coord[i])-np.matrix(coord[j])), i, j
         print "done"

   if todo==2 or todo==3:
      L=[]
      freq=[]
      for i in range(numfiles):
         f,Lmat=readLfreq(argv[i])
         freq.append(f)
         L.append(Lmat)
      if numfiles==1:
         print "nothing to compare"
      else:
         print "                       norm(f[i]-f[j])     norm(L[i]-L[j]) i   j "
         for i in range(numfiles):
            for j in range(i+1,numfiles):
               print "norm of difference:", np.linalg.norm(np.matrix(freq[i])-np.matrix(freq[j])),"\t",\
                     np.linalg.norm(np.matrix(L[i])-np.matrix(L[j])),"\t", i, j
               print 'directly compare frequencies:'
               for k in range(len(freq[i])):
                  print freq[i][k],"\t", freq[j][k]

def ReadF(infile):
   assert os.path.isfile(infile) and os.access(infile, os.R_OK),\
            initial[0]+' is not a valid file name or not readable.'
   #test, what kind of file was given: G09, GAMESS or NWChem

   log=open("calculation", "w")
   logging=[4,log ]
   with open(infile, "r+b") as f: #open file as mmap
      mapping = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
      for line in iter(mapping.readline, ""): #go through every line and test characteristic part
         if "GAMESS" in line: #if it is found: read important quantities from file
            print "GAMESS-file"
            dim, Coord, mass, A, E=rl.ReadGAMESS(logging, infile)
            return A, Coord, E #A is force constant matrix
         elif "Gaussian(R)" in line:
            print "Gaussian-file"
            dim, Coord, mass, A, E=rl.ReadG09(logging, infile)
            return A, Coord, E #A is force constant matrix
         elif "Northwest Computational Chemistry Package (NWChem)" in line:
            print "nwchem-file"
            dim, Coord, mass, A, E=rl.ReadNWChem(logging, infile)
            return A, Coord, E #A is force constant matrix
      else: 
         print "file type not recognised"
         return 0, 0, 0

def readLfreq(infile):
   assert os.path.isfile(infile) and os.access(infile, os.R_OK),\
            infile+' is not a valid file name or not readable.'
   #test, what kind of file was given: G09, GAMESS or NWChem
   with open(infile, "r+b") as f: #open file as mmap
      mapping = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
      dim=114
      for line in iter(mapping.readline, ""): #go through every line and test characteristic part
         if "GAMESS" in line: #if it is found: read important quantities from file
            print "GAMESS-file"
            f,L=rl.getGamessLf([infile], dim)
            f=f[0][6:]
            return f,L
         elif "Gaussian(R)" in line:
            print "Gaussian-file"
            L=rl.getGaussianL([infile], dim)
            freq=rl.getGaussianf([infile], dim)
            return freq[0]*Hartree2cm_1,L
         elif "Northwest Computational Chemistry Package (NWChem)" in line:
            print "nwchem-file"
            f, L=rl.getNwchemLf([infile], dim)
            return f[0], L
      else: 
         print "file type not recognised"
         return 0,0

if __name__ == "__main__":
   main(sys.argv[1:])

version=0.2
