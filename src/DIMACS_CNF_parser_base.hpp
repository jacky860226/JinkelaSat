#pragma once
#ifndef _DIMACS_CNF_PARSER_BASE_
#define _DIMACS_CNF_PARSER_BASE_
#include <cstdio>
#include <istream>
#include <sstream>
#include <string>
#include <vector>
template <class litTy> class DIMACS_CNF_parser_base {
  virtual litTy intToLit(int var) const = 0; // user define
  virtual void readVariableNum(unsigned variableNum) {
    // user define
  }
  virtual void readClauseNum(unsigned clauseNum) {
    // user define
  }
  virtual void readClause(const std::vector<litTy> &clause) {
    // user define
  }
  virtual void init() {
    // user define
  }

public:
  void parse(std::istream &input) {
    this->init();
    std::string str, all;
    while (std::getline(input, str)) {
      if (str[0] == 'c')
        continue;
      if (str[0] == 'p') {
        unsigned variableNum, clauseNum;
        sscanf(str.c_str(), "p cnf %u %u", &variableNum, &clauseNum);
        this->readVariableNum(variableNum);
        this->readClauseNum(clauseNum);
        continue;
      }
      all += " " + str + " ";
    }
    std::stringstream ss(all);
    std::vector<litTy> clause;
    int v;
    while (ss >> v) {
      if (v == 0) {
        clause.shrink_to_fit();
        this->readClause(clause);
        clause.clear();
        continue;
      }
      clause.emplace_back(this->intToLit(v));
    }
  }
};
#endif /* _DIMACS_CNF_PARSER_BASE_ */