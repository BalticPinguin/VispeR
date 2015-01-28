#!/usr/bin/python
#include further dicts
import sys, re, mmap, numpy as np
AMU2au=1822.88839

def main(argv=None):
  
   def replace(files, freq, L):
      """ This function creates a new file (determined by files, ending with 
      ".rep" and copies the log-file (files) into it, replacing the frequencies and 
      normal modes by those calculated by smallscript.
      The function is suited to test, whether these results coincide qualitatively with the gaussian's.
   
      ** PARAMETERS: **
      logging: specifies the level of logging-outputs
      log:   i
      files: file taken as basis ('.rep' added to be used for replacements)
      freq:  frequencies to be inserted
      L:     normal modes to be inserted  
   
      no return-statements (results are written to file)
   
      **NOTE:**
      The code is originally from stevaha (http://stackoverflow.com/questions/1597649/replace-strings-in-files-by-python)
      """
      with open(files) as f:
         out_fname = files + ".kn"
         out = open(out_fname, "w")
         s=0
         t=0
         u=-3
         for line in f:
            if re.search(r'Frequencies -- [\d .-]+', line) is not None:
               t=0 #reset t, when frequencies occur
               u+=3 #
               if len(freq[s:])> 2: # there are at least three more frequencies
                  out.write(re.sub(r'Frequencies -- [\d .-]+',
                        'Frequencies --'+'    '+str("%.4f" % freq[s])+'              '\
                        +str("%.4f" % freq[s+1])+'               '+str("%.4f" % freq[s+2]), line))
               elif len(freq[s:])== 2: # there are only two frequencies left
                  out.write(re.sub(r'Frequencies -- [\d .-]+',
                        'Frequencies --'+'    '+str("%.4f" % freq[s])+'               '\
                        +str("%.4f" % freq[s+1]), line))
               elif len(freq[s:])== 1: # there is just one additional freq
                  out.write(re.sub(r'Frequencies -- [\d .-]+',
                        'Frequencies --'+'   '+str("%.4f" % freq[s]), line))
               s+=3
            elif re.search(r'[ ]+\d+[ ]+\d+[ -]+\d.\d\d[ -]+\d.\d\d+[ \d.-]+', line) is not None:
               if len(L[t][u:].T)> 2: # there are at least three more modes
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

   def ReadLog(fileN, L):
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
      mass=np.zeros(( dim,dim )) # this is an integer since dim=3*N with N=atomicity
      for j in range(len(mtemp)):
         for k in range(len(mtemp[j])):
            mass[(k+foonum)*3][(k+foonum)*3]=np.sqrt(float(mtemp[j][k]))
            mass[(k+foonum)*3+1][(k+foonum)*3+1]=np.sqrt(float(mtemp[j][k]))
            mass[(k+foonum)*3+2][(k+foonum)*3+2]=np.sqrt(float(mtemp[j][k]))
         foonum+=len(mtemp[j])
      L=mass.dot(L)
      for i in range(len(L)):
         #renormalise matrix
         norm=L[i].T.dot(L[i])
         L[i]/=norm
      return L[6:].T

   assert len(argv)==3, 'exactly one argument required.'
   try:
      infile=open(argv[0], "r")
      L=mmap.mmap(infile.fileno(), 0, prot=mmap.PROT_READ)
      infile.close()
      infile=open(argv[1], "r")
      f=mmap.mmap(infile.fileno(), 0, prot=mmap.PROT_READ)
      infile.close()
   except IOError:
      print "file", inputf, "not found."
      return 2
   #first, read f, I
   #vector consisting of lines of input-file
   F=re.findall(r"[\d \- \.]+",f)[6:]
   for i in range(len(F)):
      F[i]=float(F[i])
   #matrix consisting of input-file
   temp=re.findall(r"[\d\t\- e.]+(?=\n)",L)
   #check dimensionality of L, F
   assert len(temp)==len(F)+6, "an error occured; number of frequencies does not coincide with number of normal modes {0} {1}".format(len(temp),len(F)+6)
   L=np.zeros(( len(temp), len(temp) ))
   for i in range(len(temp)):
      try:
         L[i]=re.findall(r"[\d\-e.]+",temp[i])
      except ValueError:
         print "error in mode", i
   #unmass-weight L:
   L=ReadLog(argv[2], L)
   #print L, F to output-file
   replace(argv[2], F, L)

   
if __name__ == "__main__":
   main(sys.argv[1:])
