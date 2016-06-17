#include <vector>
#include <stdlib.h>
#include <iostream>
#include <assert.h>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif
//void putN( double*, double*, double* , double**, double** , unsigned int, unsigned int*);
void putN(int* length, double** intens, double** freq, double** mode, int n);
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
   void print(){
      for(unsigned int i=0; i<this->size(); i++){
         std::cout<<this->trans[i].intens<<"   ";
         std::cout<<this->trans[i].freq<<"   ";
         for (unsigned int j=0; j<this->n(i); j++){
            std::cout<<this->trans[i].mode[j]<<" ";
         }
         std::cout<<std::endl;
      }
      std::cout<<std::endl;
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

Spect putN1( Spect opa, unsigned int n){
   Spect npa;
   std::cout<<n<<std::endl;
   opa.print();
   // parallelise this loop:
   #pragma omp parallel for
   for( unsigned int i=0; i<opa.size(); i++){
      if (opa(i).intens<5e-6)
         continue;
      // be carefull that only one process writes at the same time
      #pragma omp critical
      npa.add(opa(i));
      if (n<=1)
         continue;
      Spect tempSpect;
      tempSpect.add(opa(i)); // add one transition, all combinations with it later.
      transition tempTrans;
      for(unsigned int j=0; j<opa.size(); j++){
         // add all possible combinations with previously computed
         if (opa.n(j)>1){ //opa[j].mode.size()
            // if transition j is already a combination transition
            continue;
         }
         if (opa(j).mode[0]==0){
            //no combinations with 0-0 transitions
            continue;
         }
         bool equal=false;
         for(unsigned int k=0; k<opa.n(i); k++){
            if (opa.trans[i].mode[k]>=opa.trans[j].mode[0]){
               // no combinations with transitions where current transition is involved.
               // Also, take care that no combination is counted twice, therefore 
               //   skip this if the involved modes are larger than those to come.
               equal=true;
               break;
            }
         }
         if (equal)
            continue;
         // else: i /j are valid combination:
         tempTrans.intens=opa.trans[i].intens*opa.trans[j].intens;
         tempTrans.freq=opa.trans[i].freq+opa.trans[j].freq;
         tempTrans.mode.resize(opa.n(i)+opa.n(j));
         assert(tempTrans.mode.size()==opa.n(i)+1); // assured by first if() condition.
         tempTrans.mode[0]=opa.trans[j].mode[0]; // opa.n(j)==1
         for(unsigned int k=0; k<opa.n(i); k++){
            tempTrans.mode[k+1]=opa.trans[i].mode[k];
         }
         tempSpect.add(tempTrans);
         // now, overwrite array tempSpect with the computed combinations.
      }
      tempSpect=putN1(tempSpect, n-1);
      // be carefull that only one process writes at the same time
      #pragma omp critical
      {
      //tempSpect(0)=npa(i), therefore leave it away!
      for(unsigned int k=1; k<tempSpect.size(); k++)
         npa.add(tempSpect(k));
      }
   }
   return npa;
}

Spect putN2( Spect n_1pa, unsigned int n, Spect opa){
   Spect npa;
   //std::cout<<n<<std::endl;
   //n_1pa.print();
   // parallelise this loop:
   #pragma omp parallel for
   for( unsigned int i=0; i<n_1pa.size(); i++){
      if (n_1pa(i).intens<5e-6)
         continue;
      // be carefull that only one process writes at the same time
      #pragma omp critical
      npa.add(n_1pa(i));
      if (n<=1)
         continue;
      Spect tempSpect;
      transition tempTrans;
      for(unsigned int j=0; j<opa.size(); j++){
         // add all possible combinations with previously computed
         if (opa(j).mode[0]==0){
            //no combinations with 0-0 transitions
            continue;
         }
         bool equal=false;
         for(unsigned int k=0; k<n_1pa.n(i); k++){
            if (n_1pa.trans[i].mode[k]<=opa.trans[j].mode[0]){
               // no combinations with transitions where current transition is involved.
               // Also, take care that no combination is counted twice, therefore 
               //   skip this if the involved modes are smaller than those to come.
               // Also takes care that no combination with 0 is done.
               equal=true;
               break;
            }
         }
         if (equal)
            continue;
         // else: i /j are valid combination:
         tempTrans.intens=n_1pa.trans[i].intens*opa.trans[j].intens;
         tempTrans.freq=n_1pa.trans[i].freq+opa.trans[j].freq;
         tempTrans.mode.resize(n_1pa.n(i)+1);
         tempTrans.mode[0]=opa.trans[j].mode[0]; // opa.n(j)==1
         for(unsigned int k=0; k<n_1pa.n(i); k++){
            tempTrans.mode[k+1]=n_1pa.trans[i].mode[k];
         }
         tempSpect.add(tempTrans);
         // now, overwrite array tempSpect with the computed combinations.
      }
      if (tempSpect.size()>1){
         tempSpect=putN2(tempSpect, n-1, opa);
         // be carefull that only one process writes at the same time
         #pragma omp critical
         {
         for(unsigned int k=0; k<tempSpect.size(); k++)
            npa.add(tempSpect(k));
         }
      }
   }
   return npa;
}

void putN( int* length, double** intens, double** freq, double** mode, int n){
   // set-up new variables
   Spect opa;
   unsigned int leng = *length ;
   double* opa_int = *intens;
   double* opa_fre = *freq;
   double* opa_mod = *mode;
   transition trans;
   trans.mode.resize(1);
   // add given quantities into structures needed here:
   for(unsigned int i=0; i<leng; i++){
      trans.intens=opa_int[i];
      trans.freq=opa_fre[i];
      trans.mode[0]=opa_mod[i];
      opa.add(trans);
   }
   //std::cout<<"OPA:"<<std::endl;
   //opa.print();
   // calculate n-particle spectrum:
   Spect npa=putN2(opa, n, opa);
   //std::cout<<std::endl;
   //std::cout<<"NPA:"<<std::endl;
   //npa.print();
   
   leng=npa.size();
   double* npa_int, * npa_fre ;
   npa_int = (double*)malloc(sizeof(double) * leng);
   npa_fre = (double*)malloc(sizeof(double) * leng);
   for(unsigned int i=0; i<leng; i++){
      npa_int[i]=npa.trans[i].intens;
      npa_fre[i]=npa.trans[i].freq;
   }
   *length=leng;
   *intens=npa_int;
   *freq = npa_fre;
}
