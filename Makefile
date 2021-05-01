CXX := g++
TARGET := JinkelaSat
CXXFLAGS := -std=c++11 -O3 -Wall -Wextra
INCLUDE := src/include
SRC_DIRS := src
SRCS := $(wildcard $(SRC_DIRS:=/*.cpp))
OBJS := $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I $(INCLUDE) -MMD -c $< -o $@

clean:
	rm -rf $(TARGET) $(OBJS) $(DEPS)

BENCHMARK_DIRS := benchmarks/SAT\
                  benchmarks/UNSAT
BENCHMARK_DIRS := $(wildcard $(BENCHMARK_DIRS:=/*.cnf))

test: $(TARGET)
	@for files in $(BENCHMARK_DIRS) ; do \
        echo ; \
        echo test on case \"$$files\" ; \
		time -v ./$(TARGET) < $$files ; \
    done

.PHONY: all clean test
-include $(DEPS)