#include"src/JinkelaSat.hpp"
#include"src/DIMACS_CNF_parser_base.hpp"

#include<iostream>
using std::cin;
using std::cout;

template<class litTy=SimplifyMiniSat::Lit>
class MyParser: public DIMACS_CNF_parser_base<litTy>{
    litTy intToLit(int var) const override{
        return litTy(std::abs(var)-1, var < 0);
    }
    void readClause(const std::vector<litTy> &clause) override{
        solver.addClause(clause);
    }
    SimplifyMiniSat::JinkelaSat &solver;
public:
    MyParser(SimplifyMiniSat::JinkelaSat &solver):solver(solver){}
};

int main(){
    SimplifyMiniSat::JinkelaSat solver;
    MyParser<> parser(solver);
    parser.parse(cin);
    auto res = solver.solve();
    if(res == SimplifyMiniSat::Status::True){
        cout << "SAT\n";
        size_t vid = 1;
        for(auto v: solver.getModel()){
            cout << (v==SimplifyMiniSat::Status::False ? "-" : "") << vid++ << ' ';
        }
        cout << "0\n";
    }else{
        cout << "UNSAT\n";
    }
    return 0;
}