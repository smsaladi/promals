#include "header_cpp.h"
#include "smallnumber.h"

smallnumber operator +(const smallnumber & oper1, const smallnumber & oper2) {
    int i1, i2, k;
    double b;

    i1 = oper1.INT;
    i2 = oper2.INT;

    if(i1<i2) k = i2;
    else k = i1;
    b = exp(oper1.INT-k+oper1.DEC) + exp(oper2.INT-k+oper2.DEC);
    b = log(b);
    b = k + b;

    return smallnumber(b, 1);
    //fprintf(stdout, "%d %f\n", sum->INT, sum->DEC);
}

smallnumber operator *(const sm & a1, const sm & a2) {

    double b;

    b = a1.INT + a1.DEC + a2.INT + a2.DEC;

    return smallnumber(b, 1);

    //fprintf(stdout, "%d %f\n", mul->INT, mul->DEC);

}

smallnumber operator /(const sm & a1, const sm & a2) {
 
    double b;

    b = a1.INT + a1.DEC - a2.INT - a2.DEC;

    return smallnumber(b, 1);
}

int operator > (const sm & a1, const sm & a2) {

    double b;

    b = a1.INT + a1.DEC - a2.INT - a2.DEC;
	
    if(b>0) return 1;
    else return 0;
}


int operator <(const sm & a1, const sm & a2) {

    double b;

    b = a1.INT + a1.DEC - a2.INT - a2.DEC;
	
    if(b<0) return 1;
    else return 0;
}

int operator ==(const sm & a1, const sm & a2) {

    double b;

    b = a1.INT + a1.DEC - a2.INT - a2.DEC;
	
    if(b=0) return 1;
    else return 0;
}

int operator >=(const sm & a1, const sm & a2) {

    double b;

    b = a1.INT + a1.DEC - a2.INT - a2.DEC;
	
    if(b>=0) return 1;
    else return 0;
}

int operator <=(const sm & a1, const sm & a2) {

    double b;

    b = a1.INT + a1.DEC - a2.INT - a2.DEC;
	
    if(b<=0) return 1;
    else return 0;
}

// old stuff
void log2sm(double t, sm *a) {

    int i;
    double b;

    i =  (int) (t);
    b = t - i;
    a->INT = i;
    a->DEC = b;

    //fprintf(stdout, "%d %f\n", a->INT, a->DEC);

}  

void d2sm(double t, sm *a)  {
    int i;
    double b;

    if(t>0) t = log(t);
    else t = -500;

    i =  (int) (t);
    b = t - i;
    a->INT = i;
    a->DEC = b;

    //fprintf(stdout, "%d %f\n", a->INT, a->DEC);
}

double sm2d(sm *a) {

   double b;

   b = a->INT+a->DEC;
   return (exp(b) );
}

void sum2sm(sm *a1, sm *a2, sm *sum) {

    int i1, i2, k;
    double b;

    i1 = a1->INT;
    i2 = a2->INT;

    if(i1<i2) k = i2;
    else k = i1; 
    b = exp(a1->INT-k+a1->DEC) + exp(a2->INT-k+a2->DEC);
    b = log(b);
	
    b = k + b;

    log2sm(b, sum);

    //fprintf(stdout, "%d %f\n", sum->INT, sum->DEC);

}

void mul2sm(sm *a1, sm *a2, sm *mul) {

    double b;

    b = a1->INT + a1->DEC + a2->INT + a2->DEC;

    log2sm(b, mul);

    //fprintf(stdout, "%d %f\n", mul->INT, mul->DEC);

}

void div2sm(sm *a1, sm *a2, sm *div) {
 
    double b;

    b = a1->INT + a1->DEC - a2->INT - a2->DEC;

    log2sm(b, div);

}

void cpsm( sm *a1, sm *a2) {

     a1->INT = a2->INT;
     a1->DEC = a2->DEC;

}
