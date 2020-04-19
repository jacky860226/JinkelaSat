#pragma once
#ifndef _DIMACS_CNF_PARSER_BASE_
#define _DIMACS_CNF_PARSER_BASE_
#include<istream>
#include<vector>
#include<string>
#include<sstream>
#include<cstdio>
template<class litTy>
class DIMACS_CNF_parser_base{
    virtual litTy intToLit(int var) const = 0; // user define
    virtual litTy literalGetNeg(litTy lit) const = 0; // user define
    virtual void readVariableNum(unsigned variableNum){
        // user define
    }
    virtual void readClauseNum(unsigned clauseNum){
        // user define
    }
    virtual void readClause(const std::vector<litTy> &clause){
        // user define
    }
    virtual void init(){
        // user define
    }
public:
    void parse(std::istream &input){
        this->init();
        std::string str, all;
        while(std::getline(input, str)){
            if(str[0]=='c') continue;
            if(str[0]=='p'){
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
        bool redundant = false;
        while(ss >> v){
            if(v==0){
                clause.shrink_to_fit();
                if(!redundant) this->readClause(clause);
                redundant = false;
                clause.clear();
                continue;
            }else if(redundant){
                continue;
            }
            litTy literal = this->intToLit(v);
            bool repeat = false;
            for(litTy l: clause){
                if(l == literal){
                    repeat = true;
                    break;
                }
                if(l == this->literalGetNeg(literal)){
                    redundant = true;
                    break;
                }
            }
            if(!repeat) clause.emplace_back(literal);
        }
    }
};
#endif /* _DIMACS_CNF_PARSER_BASE_ */