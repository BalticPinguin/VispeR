#include <vector>
#include <iostream>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif
unsigned int putN( double* intens, double* freq, int* mode, unsigned int n, unsigned int len);
//double** putN( double* intens, double* freq, int* mode, unsigned int n, unsigned int len);
#ifdef __cplusplus
}
#endif

struct transition{
   double intens;
   double freq;
   vector<unsigned int> mode;
};

struct Spect{
   vector<transition> trans;

   unsigned int n(unsigned int i){return trans[i].mode.size();}
   unsigned int size(){return trans.size();}
   void add(transition newtrans){
      this->trans.push_back(newtrans);
   }
   transition operator()(unsigned int i){
     return trans[i];
   }
   void operator=(Spect rhs){
      trans=rhs.trans;
   }
};

unsigned int min(vector<unsigned int> values ){
   unsigned int min=values[0];
   for(unsigned int i=1; i<values.size(); i++){
      if(min>values[i])
         min=values[i];
   }
   return min;
}

Spect putN2( Spect opa, unsigned int n){
   Spect npa;
   // parallelise this loop:
   #pragma omp parallel for
   for( unsigned int i=0; i<opa.size(); i++){
      if (opa(i).intens<5e-6)
         continue;
      // be carefull that only one process writes at the same time
      #pragma omp critical
      npa.add(opa(i));
      if (n<1)
         continue;
      Spect tempSpect;
      transition tempTrans;
      for(unsigned int j=0; j<opa.size(); j++){
         // add all possible combinations with previously computed
         if (opa.n(j)>1)
            // if transition j is already a combination transition
            continue;
         if (min(opa(j).mode)==0)
            //no combinations with 0-0 transitions
            continue;
         for(unsigned int k=0; k<opa.n(i); k++){
            if (opa.trans[i].mode[k]>=min(opa.trans[j].mode))
               // no combinations with transitions where current transition is involved.
               // Also, take care that no combination is counted twice, therefore 
               //   skip this if the involved modes are larger than those to come.
               continue;
         }
         tempTrans.intens=opa.trans[i].intens*opa.trans[j].intens;
         tempTrans.freq=opa.trans[i].freq+opa.trans[j].freq;
         tempTrans.mode.resize(opa.n(i)+opa.n(j));
         for(unsigned int k=0; k<opa.n(i); k++){
            tempTrans.mode[k]=opa.trans[i].mode[k];
         }
         for(unsigned int k=opa.n(i); k<opa.n(j)+opa.n(i); k++){
            tempTrans.mode[k]=opa.trans[j].mode[k];
         }
         tempSpect.add(tempTrans);
         // now, overwrite array tempSpect with the computed combinations.
         tempSpect=putN2(tempSpect, n-1);
         // be carefull that only one process writes at the same time
         #pragma omp critical
         {
         for(unsigned int k=opa.n(i); k<tempSpect.size(); k++){
            npa.add(tempSpect(k));
         }
         }
      }
   }
   return npa;
}

unsigned int putN( double* intens, double* freq, int* mode, unsigned int n, unsigned int len){
   // set-up new variables
   Spect opa;
   transition trans;
   trans.mode.resize(1);
   // add given quantities into structures needed here:
   for(unsigned int i=0; i<len; i++){
      trans.intens=intens[i];
      trans.freq=freq[i];
      trans.mode[0]=mode[i];
      opa.add(trans);
   }
   // calculate n-particle spectrum:
   Spect npa=putN2(opa, n);
   len=npa.size();
   //---> this is memory-magic. One should not do it this way!
   intens=(double*) realloc(intens,len*sizeof(double));
   freq=(double*) realloc(freq, len*sizeof(double));
   //NPA[3][0]=len; // I want to have information on length as well.
   for(unsigned int i=0; i<len; i++){
      intens[i]=npa.trans[i].intens;
      freq[i]=npa.trans[i].freq;
   }
   return len;
   // how can i free the memeory of spectrum? --> is there an implicit descructor?
}

int main(){
   return 0;
}
