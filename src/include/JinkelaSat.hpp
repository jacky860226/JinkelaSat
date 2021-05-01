/************************************************************
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010  Niklas Sorensson
Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:
The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 ************************************************************/
#pragma once
#ifndef _JINKELA_SAT_
#define _JINKELA_SAT_
#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
namespace MiniSat {

// Luby, A. Sinclair and D. Zuckerman, Optimal speedup of Las Vegas algorithms,
// Information Processing
inline double luby(double y, int x) {
  int size, seq;
  for (size = 1, seq = 0; size < x + 1; seq++, size = 2 * size + 1)
    ;
  while (size - 1 != x) {
    size = (size - 1) >> 1;
    seq--;
    x = x % size;
  }
  return pow(y, seq);
}

struct Literal {
  int Var;
  Literal() : Var(-1) {}
  Literal(int Var, bool Sign = false) : Var(Var * 2 + (int)Sign) {}
  Literal operator^(bool b) const {
    Literal Result;
    Result.Var = Var ^ (int)b;
    return Result;
  }
  Literal operator~() const { return operator^(1); }
  bool sign() const { return Var & 1; }
  int var() const { return Var >> 1; }
  operator int() const { return Var; }
  bool isUndef() const { return Var == -1; }
};

class Clause : public std::vector<Literal> {
  bool _learnt;
  float act;

public:
  Clause(const std::vector<Literal> &ps, bool learnt)
      : std::vector<Literal>(ps), _learnt(learnt), act(0) {}
  float &activity() { return act; }
  bool learnt() const { return _learnt; }
};
using ClausePtr = std::shared_ptr<Clause>;

enum Status : int8_t { False = 0, True = 1, Undef = 2 };

class JinkelaSat {
  struct VarData {
    ClausePtr reason;
    size_t level;
    VarData(ClausePtr r, size_t l) : reason(r), level(l) {}
  };

  struct Watcher {
    ClausePtr cref;
    Literal blocker;
    Watcher() {}
    Watcher(ClausePtr cr, Literal p) : cref(cr), blocker(p) {}
    bool operator==(ClausePtr cr) const { return cref == cr; }
  };

  size_t qhead;
  std::vector<Literal> trail;
  std::vector<ClausePtr> clauses;
  std::vector<ClausePtr> learnts;
  std::vector<int> trail_lim;
  std::vector<Status> assigns;
  std::vector<VarData> vardata;
  std::vector<double> activity;
  std::vector<bool> polarity;
  std::vector<uint8_t> seen;
  std::set<std::pair<double, int>> order_heap;
  std::unordered_map<int, std::vector<Watcher>> watches;

  std::vector<Status> model;

  const double learntsize_factor = 1.0 / 3.0;
  const double restart_inc = 2.0;
  const size_t restart_first = 100;
  const size_t learntsize_adjust_start_confl = 100;
  const double learntsize_adjust_inc = 1.5;
  const double learntsize_inc = 1.1;
  const double var_decay = 0.95;
  const double clause_decay = 0.999;
  double cla_inc;
  double var_inc;

  double max_learnts;
  double learntsize_adjust_confl;
  size_t learntsize_adjust_cnt;

  int simpDB_assigns;
  int64_t simpDB_props;

  uint64_t next_var, num_clauses, clauses_literals, learnts_literals;
  bool ok;

  void insertVarOrder(int x) {
    if (!order_heap.count({activity[x], x}))
      order_heap.emplace(activity[x], x);
  }

  int newVar() {
    int v = next_var++;
    assigns.emplace_back(Status::Undef);
    vardata.emplace_back(nullptr, 0);
    activity.emplace_back(0);
    seen.emplace_back(0);
    polarity.emplace_back(true);
    insertVarOrder(v);
    return v;
  }

  size_t decisionLevel() const { return trail_lim.size(); }
  void newDecisionLevel() { trail_lim.emplace_back(trail.size()); }
  size_t nAssigns() const { return trail.size(); }
  Status value(int x) const { return assigns[x]; }
  Status value(Literal p) const {
    if (value(p.var()) == Status::Undef)
      return Status::Undef;
    return Status(int8_t(value(p.var())) ^ int8_t(p.sign()));
  }
  ClausePtr reason(int x) const { return vardata[x].reason; }
  size_t level(int x) const { return vardata[x].level; }
  bool locked(const ClausePtr c) const {
    return value((*c)[0]) == Status::True && reason((*c)[0].var()) &&
           reason((*c)[0].var()) == c;
  }

  void claBumpActivity(Clause &c) {
    if ((c.activity() += cla_inc) > 1e20) {
      for (auto l : learnts)
        l->activity() *= 1e-20;
      cla_inc *= 1e-20;
    }
  }
  void varBumpActivity(int v) {
    auto ori = activity[v];
    if ((activity[v] += var_inc) > 1e100) {
      for (auto &act : activity)
        act *= 1e-100;
      var_inc *= 1e-100;
      rebuildOrderHeap();
    } else {
      if (order_heap.count({ori, v})) {
        order_heap.erase({ori, v});
        order_heap.emplace(activity[v], v);
      }
    }
  }
  void claDecayActivity() { cla_inc *= (1 / clause_decay); }
  void varDecayActivity() { var_inc *= (1 / var_decay); }

  void rebuildOrderHeap() {
    order_heap.clear();
    for (size_t v = 0; v < nVars(); ++v)
      if (value(v) == Status::Undef)
        order_heap.emplace(activity[v], v);
  }

  void watcheErase(std::vector<Watcher> &watche, ClausePtr cr) {
    auto it = std::find(watche.begin(), watche.end(), cr);
    if (it == watche.end())
      return;
    watche.erase(it);
  }

  void attachClause(ClausePtr cr) {
    const Clause &c = *cr;
    watches[~c[0]].emplace_back(cr, c[1]);
    watches[~c[1]].emplace_back(cr, c[0]);
    if (c.learnt())
      learnts_literals += c.size();
    else
      num_clauses++, clauses_literals += c.size();
  }

  void detachClause(ClausePtr cr) {
    const Clause &c = *cr;
    watcheErase(watches[~c[0]], cr);
    watcheErase(watches[~c[1]], cr);
    if (c.learnt())
      learnts_literals -= c.size();
    else
      num_clauses++, clauses_literals -= c.size();
  }
  void removeClause(ClausePtr &cr) {
    Clause &c = *cr;
    detachClause(cr);
    if (locked(cr))
      vardata[c[0].var()].reason = nullptr;
    cr = nullptr;
  }

  bool satisfied(const Clause &c) const {
    for (auto l : c)
      if (value(l) == Status::True)
        return true;
    return false;
  }
  void removeSatisfied(std::vector<ClausePtr> &cs) {
    size_t j = 0;
    for (auto c : cs) {
      if (satisfied(*c))
        removeClause(c);
      else {
        for (size_t k = 2; k < c->size(); ++k)
          if (value(c->at(k)) == Status::False) {
            c->at(k--) = c->back();
            c->pop_back();
          }
        cs[j++] = c;
      }
    }
    cs.resize(j);
  }

  void uncheckedEnqueue(Literal p, ClausePtr from = nullptr) {
    assigns[p.var()] = Status(!p.sign());
    vardata[p.var()] = VarData(from, decisionLevel());
    trail.emplace_back(p);
  }

  void cancelUntil(size_t level) {
    if (decisionLevel() > level) {
      for (int c = int(trail.size()) - 1; c >= trail_lim[level]; c--) {
        int x = trail[c].var();
        assigns[x] = Status::Undef;
        polarity[x] = trail[c].sign();
        insertVarOrder(x);
      }
      qhead = trail_lim[level];
      trail.resize(trail_lim[level]);
      trail_lim.resize(level);
    }
  }

  ClausePtr propagate() {
    ClausePtr confl = nullptr;
    int num_props = 0;
    while (qhead < trail.size()) {
      Literal p = trail[qhead++];
      std::vector<Watcher> &ws = watches[p];
      std::vector<Watcher>::iterator i, j;
      for (i = j = ws.begin(); i != ws.end();) {
        Literal blocker = i->blocker;
        if (value(blocker) == Status::True) {
          *j++ = *i++;
          continue;
        }
        auto cr = i->cref;
        Clause &c = *cr;
        if (c[0] == ~p)
          std::swap(c[0], c[1]);
        ++i;

        Literal first = c[0];
        Watcher w = Watcher(cr, first);
        if (first != blocker && value(first) == Status::True) {
          *j++ = w;
          continue;
        }

        for (size_t k = 2; k < c.size(); k++) {
          if (value(c[k]) != Status::False) {
            std::swap(c[1], c[k]);
            watches[~c[1]].emplace_back(w);
            goto NextClause;
          }
        }
        *j++ = w;
        if (value(first) == Status::False) {
          confl = cr;
          qhead = trail.size();
          while (i != ws.end())
            *j++ = *i++;
        } else
          uncheckedEnqueue(first, cr);
      NextClause:;
      }
      ws.resize(j - ws.begin());
    }
    simpDB_props -= num_props;
    return confl;
  }

  bool litRedundant(Literal p, std::vector<Literal> &analyze_toclear) {
    enum {
      seen_undef = 0,
      seen_source = 1,
      seen_removable = 2,
      seen_failed = 3
    };
    auto c = reason(p.var());
    std::vector<std::pair<uint32_t, Literal>> Stack;
    for (size_t i = 1;; ++i) {
      if (i < c->size()) {
        Literal l = (*c)[i];
        if (level(l.var()) == 0 || seen[l.var()] == seen_source ||
            seen[l.var()] == seen_removable)
          continue;
        if (reason(l.var()) == nullptr || seen[l.var()] == seen_failed) {
          Stack.emplace_back(0, p);
          for (auto v : Stack)
            if (seen[v.second.var()] == seen_undef) {
              seen[v.second.var()] = seen_failed;
              analyze_toclear.emplace_back(v.second);
            }
          return false;
        }
        Stack.emplace_back(i, p);
        i = 0, p = l;
        c = reason(p.var());
      } else {
        if (seen[p.var()] == seen_undef) {
          seen[p.var()] = seen_removable;
          analyze_toclear.emplace_back(p);
        }
        if (Stack.empty())
          break;
        std::tie(i, p) = Stack.back();
        c = reason(p.var());
        Stack.pop_back();
      }
    }
    return true;
  }

  void analyze(ClausePtr confl, std::vector<Literal> &out_learnt,
               int &out_btlevel) {
    int pathC = 0;
    Literal p;
    out_learnt.emplace_back();
    int index = int(trail.size()) - 1;
    do {
      Clause &c = *confl;
      if (c.learnt())
        claBumpActivity(c);
      for (size_t j = p.isUndef() ? 0 : 1; j < c.size(); ++j) {
        Literal q = c[j];
        if (!seen[q.var()] && level(q.var()) > 0) {
          varBumpActivity(q.var());
          seen[q.var()] = 1;
          if (level(q.var()) >= decisionLevel())
            pathC++;
          else
            out_learnt.emplace_back(q);
        }
      }
      while (!seen[trail[index--].var()])
        ;
      p = trail[index + 1];
      confl = reason(p.var());
      seen[p.var()] = 0;
      pathC--;
    } while (pathC > 0);
    out_learnt[0] = ~p;
    auto analyze_toclear = out_learnt;
    size_t j = 1;
    for (size_t i = 1; i < out_learnt.size(); ++i)
      if (reason(out_learnt[i].var()) == nullptr ||
          !litRedundant(out_learnt[i], analyze_toclear))
        out_learnt[j++] = out_learnt[i];
    out_learnt.resize(j);
    if (out_learnt.size() > 1) {
      int max_i = 1;
      for (size_t i = 2; i < out_learnt.size(); ++i)
        if (level(out_learnt[i].var()) > level(out_learnt[max_i].var()))
          max_i = i;
      std::swap(out_learnt[1], out_learnt[max_i]);
      out_btlevel = level(out_learnt[1].var());
    } else
      out_btlevel = 0;
    for (auto v : analyze_toclear)
      seen[v.var()] = 0;
  }

  bool simplify() {
    if (propagate() != nullptr)
      return false;
    if (int(nAssigns()) == simpDB_assigns || (simpDB_props > 0))
      return true;
    removeSatisfied(learnts);
    rebuildOrderHeap();
    simpDB_assigns = nAssigns();
    simpDB_props = clauses_literals + learnts_literals;
    return true;
  }

  static bool reduceDB_cmp(ClausePtr x, ClausePtr y) {
    return x->size() > 2 && (y->size() == 2 || x->activity() < y->activity());
  }
  void reduceDB() {
    size_t j = 0;
    double extra_lim = cla_inc / learnts.size();
    std::sort(learnts.begin(), learnts.end(), reduceDB_cmp);
    for (size_t i = 0; i < learnts.size(); i++) {
      auto c = learnts[i];
      if (c->size() > 2 && !locked(c) &&
          (i < learnts.size() / 2 || c->activity() < extra_lim)) {
        removeClause(c);
      } else
        learnts[j++] = c;
    }
    learnts.resize(j);
  }

  Literal pickBranchLit() {
    const int VarUndef = -1;
    int next = VarUndef;
    while (next == VarUndef || value(next) != Status::Undef)
      if (order_heap.empty()) {
        next = VarUndef;
        break;
      } else {
        next = (*--order_heap.end()).second;
        order_heap.erase(--order_heap.end());
      }
    if (next == VarUndef)
      return Literal();
    return Literal(next, polarity[next]);
  }

  Status search(int nof_conflicts) {
    int conflictC = 0;
    for (;;) {
      auto confl = propagate();
      if (confl != nullptr) {
        ++conflictC;
        if (decisionLevel() == 0)
          return Status::False;
        std::vector<Literal> learnt_clause;
        int backtrack_level;
        analyze(confl, learnt_clause, backtrack_level);
        cancelUntil(backtrack_level);
        if (learnt_clause.size() == 1) {
          uncheckedEnqueue(learnt_clause[0]);
        } else {
          auto cr = std::make_shared<Clause>(learnt_clause, true);
          learnts.emplace_back(cr);
          attachClause(cr);
          claBumpActivity(*cr);
          uncheckedEnqueue(learnt_clause[0], cr);
        }
        varDecayActivity();
        claDecayActivity();
        if (--learntsize_adjust_cnt == 0) {
          learntsize_adjust_confl *= learntsize_adjust_inc;
          learntsize_adjust_cnt = (int)learntsize_adjust_confl;
          max_learnts *= learntsize_inc;
        }
      } else {
        if (nof_conflicts >= 0 && conflictC >= nof_conflicts) {
          cancelUntil(0);
          return Status::Undef;
        }
        if (decisionLevel() == 0 && !simplify())
          return Status::False;
        if (int(learnts.size()) - int(nAssigns()) >= max_learnts)
          reduceDB();
        Literal next = pickBranchLit();
        if (next.isUndef())
          return Status::True;
        newDecisionLevel();
        uncheckedEnqueue(next);
      }
    }
  }

public:
  JinkelaSat()
      : qhead(0), cla_inc(1.0), var_inc(1), simpDB_assigns(-1), simpDB_props(0),
        next_var(0), num_clauses(0), clauses_literals(0), learnts_literals(0),
        ok(true) {}
  size_t nClauses() const { return num_clauses; }
  size_t nVars() const { return next_var; }
  bool addClause(std::vector<Literal> ps) {
    if (!ok)
      return false;
    for (Literal p : ps)
      while (p.var() >= int(nVars()))
        newVar();
    std::sort(ps.begin(), ps.end());
    Literal p;
    size_t j = 0;
    for (Literal l : ps) {
      if (value(l) == Status::True || l == ~p) {
        return true;
      } else if (value(l) != Status::False && l != p) {
        ps[j++] = p = l;
      }
    }
    ps.resize(j), ps.shrink_to_fit();
    if (ps.empty())
      return ok = false;
    else if (ps.size() == 1) {
      uncheckedEnqueue(ps[0]);
      return ok = (propagate() == nullptr);
    } else {
      ClausePtr cr = std::make_shared<Clause>(ps, false);
      clauses.emplace_back(cr);
      attachClause(cr);
    }
    return true;
  }
  Status solve() {
    if (!ok)
      return Status::False;
    max_learnts = nClauses() * learntsize_factor;
    learntsize_adjust_confl = learntsize_adjust_start_confl;
    learntsize_adjust_cnt = learntsize_adjust_start_confl;
    Status status = Status::Undef;
    int curr_restarts = 0;
    do {
      double rest_base = luby(restart_inc, curr_restarts);
      status = search(rest_base * restart_first);
      ++curr_restarts;
    } while (status == Status::Undef);
    model = assigns;
    cancelUntil(0);
    return status;
  }
  const std::vector<Status> &getModel() const { return model; }
};

} // namespace MiniSat
#endif /* _JINKELA_SAT_ */