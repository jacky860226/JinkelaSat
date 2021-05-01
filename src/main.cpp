#include "DIMACS_CNF_Parser.hpp"
#include "JinkelaSat.hpp"
#include <iostream>
#include <ostream>

void adapteInput(MiniSat::JinkelaSat &Solver, const DIMACS::CNF::Input &Input) {
  for (auto &ClauseRaw : Input.getClauses()) {
    std::vector<MiniSat::Literal> Clause;
    for (int RawLit : ClauseRaw) {
      Clause.emplace_back(std::abs(RawLit) - 1, RawLit < 0);
    }
    Solver.addClause(std::move(Clause));
  }
}

void to_ostream(MiniSat::JinkelaSat &Solver, const MiniSat::Status Result,
                std::ostream &Out) {
  if (Result == MiniSat::Status::True) {
    Out << "SAT\n";
    size_t VarCnt = 1;
    for (auto VarStatus : Solver.getModel()) {
      Out << (VarStatus == MiniSat::Status::False ? "-" : "") << VarCnt++
          << ' ';
    }
    Out << "0\n";
  } else {
    Out << "UNSAT\n";
  }
}

int main() {
  DIMACS::CNF::Parser Parser;
  auto Input = Parser.parse(std::cin);

  MiniSat::JinkelaSat Solver;
  adapteInput(Solver, Input);

  auto Result = Solver.solve();
  to_ostream(Solver, Result, std::cout);

  return 0;
}