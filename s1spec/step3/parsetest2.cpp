#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "capd/capdlib.h"
using namespace std; using namespace capd;
int main(int argc, char* argv[]) {
  ifstream in(argv[1]); string line; int skip = atoi(argv[2]);
  for (int i = 0; i < skip; ++i) getline(in, line);
  getline(in, line);
  string head = line.substr(0, line.find("fun:") + 4);
  string funs = line.substr(line.find("fun:") + 4);
  vector<string> comps; int depth = 0; string cur;
  for (char ch : funs) {
    if (ch == '(') depth++; if (ch == ')') depth--;
    if ((ch == ',' || ch == ';') && depth == 0) { comps.push_back(cur); cur.clear(); continue; }
    cur += ch;
  }
  for (size_t j = 0; j < comps.size(); ++j) {
    string s = head;
    for (size_t i = 0; i < comps.size(); ++i) { s += (i == j ? comps[i] : string("0")); s += (i + 1 < comps.size() ? "," : ";"); }
    try { IMap f(s.c_str()); } catch (std::exception& e) { cout << "component " << j << " FAIL: " << comps[j].substr(0, 120) << "\n"; }
  }
  cout << comps.size() << " components checked\n";
}
