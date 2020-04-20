# JinkelaSat
A Header-Only CDCL SAT Solver based on MiniSat

MiniSat: https://github.com/niklasso/minisat

## 1. How to Build
**Step 1:** Download the source code. For example,
~~~
$ git clone https://github.com/jacky860226/JinkelaSat.git
~~~

**Step 2:** Go to the project root and build by
~~~
$ cd JinkelaSat
$ g++ -std=c++11 -O3 main.cpp -o JinkelaSatSolver
~~~

### 1.1. Dependencies

* [GCC](https://gcc.gnu.org/) (version >= 7.5.0) or other working c++ compliers

## 2. How to solve cnf

Read DIMACS cnf data from **standard input**.

```
$ ./JinkelaSatSolver < SampleCNF.cnf
```

### 2.1. SATISFIABLE

If Sat Solver proves a given problem is SATISFIABLE(SAT), it prints SAT and a model for the problem.

```
SAT
-1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 13 -14 15 16 0
```

The model is separated by a space and ended with zero(0).

-   A positive integer represents `True`.
-   A negative integer represents `False`.

For example,`-1` represents `X1=False`,`13` represents `X13=True`.

### 2.2. UNSATISFIABLE

It only prints UNSAT.

```
UNSAT
```

---

This project is an assignment for **NTHU 10820CS 538100 Applied Mathematical Logic**, plagiarism is prohibited.
