/////////////////////////////////////////////////////////////////
// ScoreType.h
//
// Routines for doing math operations in PROBCONS.
/////////////////////////////////////////////////////////////////

#ifndef SCORETYPE_H
#define SCORETYPE_H

#include <cmath>
#include <algorithm>
#include <cfloat>

typedef float ScoreType;

const float LOG_ZERO = -2e20;
const float LOG_ONE = 0.0;

/////////////////////////////////////////////////////////////////
// LOG()
//
// Compute the logarithm of x.
/////////////////////////////////////////////////////////////////

inline ScoreType LOG (ScoreType x){
  return log (x);
}

/////////////////////////////////////////////////////////////////
// EXP()
//
// Computes exp(x).
/////////////////////////////////////////////////////////////////

inline ScoreType EXP (ScoreType x){
  //return exp(x);
  if (x > -2){
    if (x > -0.5){
      if (x > 0)
	return exp(x);
      return (((0.03254409303190190000*x + 0.16280432765779600000)*x + 0.49929760485974900000)*x + 0.99995149601363700000)*x + 0.99999925508501600000;
    }
    if (x > -1)
      return (((0.01973899026052090000*x + 0.13822379685007000000)*x + 0.48056651562365000000)*x + 0.99326940370383500000)*x + 0.99906756856399500000;
    return (((0.00940528203591384000*x + 0.09414963667859410000)*x + 0.40825793595877300000)*x + 0.93933625499130400000)*x + 0.98369508190545300000;
  }
  if (x > -8){
    if (x > -4)
      return (((0.00217245711583303000*x + 0.03484829428350620000)*x + 0.22118199801337800000)*x + 0.67049462206469500000)*x + 0.83556950223398500000;
    return (((0.00012398771025456900*x + 0.00349155785951272000)*x + 0.03727721426017900000)*x + 0.17974997741536900000)*x + 0.33249299994217400000;
  }
  if (x > -16)
    return (((0.00000051741713416603*x + 0.00002721456879608080)*x + 0.00053418601865636800)*x + 0.00464101989351936000)*x + 0.01507447981459420000;
  return 0;
}

const float EXP_UNDERFLOW_THRESHOLD = -4.6;
const float LOG_UNDERFLOW_THRESHOLD = 7.5;

/////////////////////////////////////////////////////////////////
// LOOKUP()
//
// Computes log (exp (x) + 1), for 0 <= x <= 7.5.
/////////////////////////////////////////////////////////////////

inline float LOOKUP (float x){
  //return log (exp(x) + 1);
  assert (x >= 0.00f);
  assert (x <= LOG_UNDERFLOW_THRESHOLD);
  //return ((-0.00653779113685f * x + 0.09537236626558f) * x + 0.55317574459331f) * x + 0.68672959851568f;
  if (x <= 1.00f) return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
  if (x <= 2.50f) return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
  if (x <= 4.50f) return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;
  //assert (x <= LOG_UNDERFLOW_THRESHOLD);
  return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;

  //return (((0.00089738532761f * x - 0.01859488697982f) * x + 0.14415772028626f) * x + 0.49515490689159f) * x + 0.69311928966454f;
}

/////////////////////////////////////////////////////////////////
// LOOKUP_SLOW()
//
// Computes log (exp (x) + 1).
/////////////////////////////////////////////////////////////////

inline float LOOKUP_SLOW (float x){
  return log (exp (x) + 1);
}

/////////////////////////////////////////////////////////////////
// MAX()
//
// Compute max of three numbers
/////////////////////////////////////////////////////////////////

inline float MAX (float x, float y, float z){
  if (x >= y){
    if (x >= z)
      return x;
    return z;
  }
  if (y >= z)
    return y;
  return z;
}

/////////////////////////////////////////////////////////////////
// LOG_PLUS_EQUALS()
//
// Add two log probabilities and store in the first argument
/////////////////////////////////////////////////////////////////

inline void LOG_PLUS_EQUALS (float &x, float y){
  if (x < y)
    x = (x == LOG_ZERO || y - x >= LOG_UNDERFLOW_THRESHOLD) ? y : LOOKUP(y-x) + x;
  else
    x = (y == LOG_ZERO || x - y >= LOG_UNDERFLOW_THRESHOLD) ? x : LOOKUP(x-y) + y;
}

/////////////////////////////////////////////////////////////////
// LOG_PLUS_EQUALS_SLOW()
//
// Add two log probabilities and store in the first argument
/////////////////////////////////////////////////////////////////

inline void LOG_PLUS_EQUALS_SLOW (float &x, float y){
  if (x < y)
    x = (x == LOG_ZERO) ? y : LOOKUP_SLOW(y-x) + x;
  else
    x = (y == LOG_ZERO) ? x : LOOKUP_SLOW(x-y) + y;
}

/////////////////////////////////////////////////////////////////
// LOG_ADD()
//
// Add two log probabilities
/////////////////////////////////////////////////////////////////

inline float LOG_ADD (float x, float y){
  //if(x==y) { cout << "x is equal to y" << endl; }
  if (x < y) {
	//cout << "x less than y" << endl;
	return (x == LOG_ZERO || y - x >= LOG_UNDERFLOW_THRESHOLD) ? y : LOOKUP(y-x) + x;
  }
  //cout << "x equal or larger than y" << endl;
  //cout << setprecision(10) << (x-y) << endl;
  return (y == LOG_ZERO || x - y >= LOG_UNDERFLOW_THRESHOLD) ? x : LOOKUP(x-y) + y;
}


/////////////////////////////////////////////////////////////////
// LOG_ADD()
//
// Add three log probabilities
/////////////////////////////////////////////////////////////////

inline float LOG_ADD (float x1, float x2, float x3){
  return LOG_ADD (x1, LOG_ADD (x2, x3));
}

/////////////////////////////////////////////////////////////////
// LOG_ADD()
//
// Add four log probabilities
/////////////////////////////////////////////////////////////////

inline float LOG_ADD (float x1, float x2, float x3, float x4){
  return LOG_ADD (x1, LOG_ADD (x2, LOG_ADD (x3, x4)));
}

/////////////////////////////////////////////////////////////////
// LOG_ADD()
//
// Add five log probabilities
/////////////////////////////////////////////////////////////////

inline float LOG_ADD (float x1, float x2, float x3, float x4, float x5){
  return LOG_ADD (x1, LOG_ADD (x2, LOG_ADD (x3, LOG_ADD (x4, x5))));
}

/////////////////////////////////////////////////////////////////
// LOG_ADD()
//
// Add siz log probabilities
/////////////////////////////////////////////////////////////////

inline float LOG_ADD (float x1, float x2, float x3, float x4, float x5, float x6){
  return LOG_ADD (x1, LOG_ADD (x2, LOG_ADD (x3, LOG_ADD (x4, LOG_ADD (x5, x6)))));
}

/////////////////////////////////////////////////////////////////
// LOG_ADD()
//
// Add seven log probabilities
/////////////////////////////////////////////////////////////////

inline float LOG_ADD (float x1, float x2, float x3, float x4, float x5, float x6, float x7){
  return LOG_ADD (x1, LOG_ADD (x2, LOG_ADD (x3, LOG_ADD (x4, LOG_ADD (x5, LOG_ADD (x6, x7))))));
}

/////////////////////////////////////////////////////////////////
// ChooseBestOfThree()
//
// Store the largest of three values x1, x2, and x3 in *x.  Also
// if xi is the largest value, then store bi in *b.
/////////////////////////////////////////////////////////////////

inline void ChooseBestOfThree (float x1, float x2, float x3, char b1, char b2, char b3, float *x, char *b){
  if (x1 >= x2){
    if (x1 >= x3){
      *x = x1;
      *b = b1;
      return;
    }
    *x = x3;
    *b = b3;
    return;
  }
  if (x2 >= x3){
    *x = x2;
    *b = b2;
    return;
  }
  *x = x3;
  *b = b3;
}

#endif
