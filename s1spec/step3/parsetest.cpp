#include <iostream>
#include <fstream>
#include <string>
#include "capd/capdlib.h"
using namespace std; using namespace capd;
int main(int argc, char* argv[]) {
  ifstream in(argv[1]); string line; int skip = atoi(argv[2]);
  for (int i = 0; i < skip; ++i) getline(in, line);
  for (int k = 0; k < 3; ++k) {
    getline(in, line);
    string head = line.substr(0, line.find("fun:") + 4);
    string funs = line.substr(line.find("fun:") + 4);
    // split top-level by commas
    int depth = 0; string cur; int idx = 0;
    for (char ch : funs) {
      if (ch == '(') depth++; if (ch == ')') depth--;
      if ((ch == ',' || ch == ';') && depth == 0) {
        try { IMap f(("time:t;par:lam,lamt,T0,kk,mm,pp,dq,ca,cb;var:" + string("x") + ";fun:0;").c_str()); } catch (...) {}
        cur.clear(); idx++; continue;
      }
      cur += ch;
    }
    try { IMap f(line.c_str()); cout << "line " << k << " OK dim " << f.dimension() << "\n"; }
    catch (std::exception& e) { cout << "line " << k << " FAIL: " << e.what() << "\n"; }
  }
}
