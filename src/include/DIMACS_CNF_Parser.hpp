#pragma once
#ifndef _DIMACS_CNF_PARSER_
#define _DIMACS_CNF_PARSER_
#include <cstdio>
#include <istream>
#include <sstream>
#include <string>
#include <vector>

namespace DIMACS {
namespace CNF {

class Input {
  const unsigned VariableNum, ClauseNum;
  const std::vector<std::vector<int>> Clauses;

public:
  Input(unsigned VariableNum, unsigned ClauseNum,
        std::vector<std::vector<int>> &&Clauses)
      : VariableNum(VariableNum), ClauseNum(ClauseNum),
        Clauses(std::move(Clauses)) {}
  Input(Input &&other)
      : VariableNum(other.VariableNum), ClauseNum(other.ClauseNum),
        Clauses(std::move(other.Clauses)) {}
  unsigned getVariableNum() const { return VariableNum; }
  unsigned getClauseNum() const { return ClauseNum; }
  const std::vector<std::vector<int>> &getClauses() const { return Clauses; }
};

class Parser {
public:
  Input parse(std::istream &input) {
    unsigned VariableNum, ClauseNum;
    std::vector<std::vector<int>> Clauses;
    std::string Line, Buffer;
    while (std::getline(input, Line)) {
      if (Line[0] == 'c')
        continue;
      if (Line[0] == 'p') {
        sscanf(Line.c_str(), "p cnf %u %u", &VariableNum, &ClauseNum);
        continue;
      }
      Buffer += " " + Line + " ";
    }
    std::stringstream SS(Buffer);
    std::vector<int> Clause;
    int Lit;
    while (SS >> Lit) {
      if (Lit == 0) {
        Clause.shrink_to_fit();
        Clauses.emplace_back(std::move(Clause));
        Clause.clear();
        continue;
      }
      Clause.emplace_back(Lit);
    }
    return Input(VariableNum, ClauseNum, std::move(Clauses));
  }
};

} // namespace CNF
} // namespace DIMACS

#endif /* _DIMACS_CNF_PARSER_ */