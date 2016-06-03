#include <vector>
#include <stdlib.h>
#include <iostream>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif
//void putN( double*, double*, double* , double**, double** , unsigned int, unsigned int*);
void putN(int* length, double** intens, double** freq, double** mode, int n);
void test(int *size, double **out1, double **out2);
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

void putN3( double* intens, double* freq, double* mode,
           double** npa_i, double** npa_f, unsigned int n, unsigned int* leng){
   // set-up new variables
   Spect opa;
   unsigned int length = *leng ;
   std::cout<<length<<std::endl;
   transition trans;
   trans.mode.resize(1);
   // add given quantities into structures needed here:
   for(unsigned int i=0; i<length; i++){
      trans.intens=intens[i];
      trans.freq=freq[i];
      trans.mode[0]=mode[i];
      opa.add(trans);
   }
   // calculate n-particle spectrum:
   Spect npa=putN2(opa, n);
   length=npa.size();
   
   double* npa_int, * npa_fre ;
   npa_int = (double*)malloc(sizeof(double) * length);
   npa_fre = (double*)malloc(sizeof(double) * length);
   for(unsigned int i=0; i<length; i++){
      npa_int[i]=npa.trans[i].intens;
      npa_fre[i]=npa.trans[i].freq;
   }
   *leng=length;
   *npa_i = npa_int;
   *npa_f = npa_fre;
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
   // calculate n-particle spectrum:
   Spect npa=putN2(opa, n);
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

void test(int *size, double **out1, double **out2) {
    int i;
    * size=10;
    double *data1, *data2;
    data1 = (double *)malloc(sizeof(double) * *size);
    data2 = (double *)malloc(sizeof(double) * *size);
    for (i = 0; i < *size; i++){
        data1[i] = i;
        data2[i] = i * 2;
    }
    //*size = N;
    *out1 = data1;
    *out2 = data2;
}

//int main(){
   //return 0;
//}
