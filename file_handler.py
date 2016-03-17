#!/usr/bin/python2
# filename: file_handler.py

# CHANGELOG
# =========
#in version 0.1.0:  
#  1) Intialised class
#

class logging():
   logfile=''
   level='0'
   numWarn=0
   width=5  # determines, how many rows of a matrix are writen aside each other.

   def __init__(self, level, logfile):
      def invokeLogging(logfile, mode="important"):
         """ initialises the logging-functionality
            **PARAMETERS**
            logfile   name of file to be used as log-file. It is expected to be an array of 
                   length 0 or one.
            mode:     5 different values are possible (see below); the lowest means: print 
                     much, higher values mean 
                     less printing

            **RETURNS**
            log:      opened file to write in
         """
         if logfile==[]:
            log=open("calculation.log", "a")
            self.logfile="calculation.log"
         else:
            s=logfile[-1].strip()
            log=open(s, "a")
            self.logfile=s
         
         #remove all white spaces and take only first character.
         mode= mode.strip()[0]
         #search, which option was set.
         if mode in ['a', "A", "0", 0]:
            logging=0
            log.write('use log-level all\n')
         elif mode in ['d', 'D', "1"]:
            logging=1
            log.write('use log-level detailed\n')
         elif mode in ['m', 'M','2']:
            logging=2
            log.write('use log-level medium\n')
         elif mode in ['i', 'I','3']:
            logging=3
         elif mode in ['s', 'S', '4']:
            logging=4
         else:
            logging=3
            log.write("logging-mode "+mode+" not recognized. Using 'important' instead\n")
         return logging, log
      
      self.logfile=logfile
      self.level, self.loghandler = invokeLogging(logfile, level)
      self.write("\n==================================================================\n"
               "=====================  output of Visper  =========================\n"
               "==================================================================\n\n")

   def __del__(self):
      """This simple function has the only task to close the log-file after 
         calculation as it is nice. 
         More over, it gives the user feedback about successfull finish;
         also nice if some warnings appeared before.
      """
      self.loghandler.close()
      #count warnings in self.loghandler:
      logf=open(self.logfile,"r")
      warnings=0
      for line in logf:
         if 'WARNING:' in line:
            warnings+=1
      logf.close()
      foo=open(self.logfile, 'a')
      if warnings==0:
         foo.write("\n==================================================================\n"
                  "=========== VISPER FINISHED OPERATION SUCCESSFULLY.  =============\n"
                  "==================================================================\n\n")
      else:
         foo.write("\n==================================================================\n"
                  "============== VISPER FINISHED WITH "+repr(warnings)+" WARNINGS.  ===============\n"
                  "==================================================================\n\n")
      foo.close()
   
   def write(self, text, level=70):
      if level>self.level:
         self.loghandler.write(text)

   def printVec(self,vec):
      """This funcion is no that tricky but better than rewriting it everywhere it is
         indeed.
      """
      num = self.width//2
      self.loghandler.write("\n")
      if len(vec)>num:
         for j in range(len(vec)//num):
            for k in range(num):
               self.loghandler.write("    %03d  %e \t"%(j+k*len(vec)//num+1, vec[j+k*len(vec)//num]))
            self.loghandler.write("\n")
      else:
         for k in range(len(vec)):
            self.loghandler.write("    %03d  %e \t"%(k+1, vec[k]))
      self.loghandler.write("\n")
   
   def printMat(self,mat):
      """Function to print matrices in a nice way to the log-file.
         To keep it general, the matrix needs to be given as argument.
      """
      k=range(0,len(mat[0]))
      s=0
      t=min(s+self.width,len(mat[0]))
      while s<len(mat[0]):
         self.loghandler.write("\n")
         for temp in range(s,t):
            self.loghandler.write("          %03d "%(k[temp]+1))
         self.loghandler.write("\n")
         for j in range(len(mat)):
            self.loghandler.write(" %03d"%(j+1))
            for temp in range(s,t):
               self.loghandler.write("  %+.5e"%(mat[j][k[temp]]))
            self.loghandler.write("\n")
         s=t
         t=min(s+self.width,len(mat[0]))

#version 0.1.0
# End of file_handler.py
