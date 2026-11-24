#include "capd/capdlib.h"
extern const char* VF;
const char* VF =
  "par:lr,li;"
  "var:p,q,k,m,a1r,a1i,a2r,a2i,a3r,a3i,b2r,b2i,b3r,b3i,c3r,c3i;"
  "fun:k*q, -1-k*p, -m, -0.25*p,"
  "-a2r, -a2i,"
  "-(a3r+b2r), -(a3i+b2i),"
  "-(-(5*k*m+0.25*q)*a1r - k*k*a2r + b3r), -(-(5*k*m+0.25*q)*a1i - k*k*a2i + b3i),"
  "-b3r, -b3i,"
  "-(-(0.25-3*m*m-k*p)*a1r + (lr*a1r - li*a1i) - k*k*b2r + c3r),"
  "-(-(0.25-3*m*m-k*p)*a1i + (lr*a1i + li*a1r) - k*k*b2i + c3i),"
  "-(-(0.25-3*m*m-k*p)*a2r + (lr*a2r - li*a2i) + (5*k*m+0.25*q)*b2r),"
  "-(-(0.25-3*m*m-k*p)*a2i + (lr*a2i + li*a2r) + (5*k*m+0.25*q)*b2i);";

