// smallnumber is decomposed to two parts: integer part 
// of the logarithm and decimal part of the logarithm

#ifndef __smallnumber_jpei_
#define __smallnumber_jpei_

typedef class smallnumber {
  

   public:
	int INT;
	double DEC;

	smallnumber() { INT = 0; DEC = 0; } // default is 1 
	// copy constructor
	smallnumber(const smallnumber & other_small_number) {
		INT = other_small_number.INT;
		DEC = other_small_number.DEC;
	}
	// construct a small number from a real number
	smallnumber(double x) { set(x); }
	// construct a small number from a logrithm real number
	smallnumber (double x, int LOG) { INT = (int)(x); DEC = x-INT; }
	~smallnumber() {}

	void set(double t) {
		if(t>0) {t = log(t); INT =  (int) (t); DEC = t - INT;}
		else {INT = -100000, DEC = 0;}
	}

	void setLog(double t) {
		INT = (int) (t);
		DEC = t - INT;
	}
		
	smallnumber & operator =(const smallnumber & oper2) {
		INT = oper2.INT; DEC = oper2.DEC;
		return (*this);
	}
   
	smallnumber & operator =(const double & oper2) {
		set(oper2);
		return (*this);
	}

	double real() { return (exp(INT+DEC) ); }

} sm;

smallnumber operator +(const smallnumber & oper1, const smallnumber & oper2);
smallnumber operator *(const smallnumber & a1, const smallnumber & a2);
smallnumber operator /(const sm & a1, const sm & a2);
int operator ==(const sm & a1, const sm & a2);
int operator >(const sm & a1, const sm & a2);
int operator <(const sm & a1, const sm & a2);
int operator >=(const sm & a1, const sm & a2);
int operator <=(const sm & a1, const sm & a2);

// old stuff
void log2sm(double t, sm *a);
void d2sm(double t, sm *a);
double sm2d(sm *a);
void sum2sm(sm *a1, sm *a2, sm *sum);
void mul2sm(sm *a1, sm *a2, sm *mul);
void div2sm(sm *a1, sm *a2, sm *div);
void cpsm( sm *a1, sm *a2);

#endif
