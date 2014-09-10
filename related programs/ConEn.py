#!/usr/bin/python
import re, mmap, sys, getopt

def usage():
   print "usage: ConEn.py -i <input-file name>"
   print "you can specify as many input-files as you want. But mark each by '-i'."

def main(argv=None):
   if argv is None:
      argv = sys.argv[1:]
   try:
      opts, args=getopt.getopt(argv, 'i:h', ["help", "input="])
   except getopt.GetoptError as err:
      print(err)
      usage()
      sys.exit(2)
   if opts in ['-h', '--help']:
      usage()
      sys.exit()
   inputs=[]
   print opts
   for opt,s in opts:
      if opt in ['-i', '--input']:
	 inputs.append(s)
      else:
	 usage()
	 sys.exit(2)
   print inputs
   En=[]
   name=[]
   for files in inputs:
      name.append(files)
      log=open(files, "r") #open file and map it for better working
      mapedlog=mmap.mmap(log.fileno(), 0, prot=mmap.PROT_READ) #map the file
      log.close
   
      bar=re.findall(r"Total Energy, [\w --()./=]+", mapedlog, re.I)
      mapedlog.close() # changes fileno() but not the result
      E=[]
      for i in range(len(bar)):# get all energies
	 E.append(re.findall(r"-+\d+.\d+", bar[i])[0]) # these are energies of the states
      En.append(E) # is this ok this way?
   
   
   for i in range(len(En)):
      for j in range(len(En[i])):
	 print repr(j) + '   ' + repr(float(En[i][j])) +'   '+ name[i]
      print ' '
      print ' '

if __name__ == "__main__":
   main(sys.argv[1:])


##equivalent bash-handling:
# export pattern="FC*.log"
# for file in $(ls ${pattern}); do grep "Total Energy" ${file} | awk '{print NR, $5}' >>plot; echo >> plot ; echo >> plot; done
