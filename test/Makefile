CC=clang
CXX=clang++
CXX1XFLAGS=-I../src -std=c++11 -O3 -march=native -mno-avx512f -mpopcnt 

# CXX1XFLAGS=-I../src -std=c++11 -g -O2 -march=native -mpopcnt -mno-avx2 -mno-avx -fno-vectorize -DBACKWARD_HAS_DW -DCOMPILE_CPU
# CXX1XFLAGS=-I../src -std=c++11 -O3 -march=native -mpopcnt -Rpass=loop-vectorize -DCOMPILE_CPU
# CXX1XFLAGS=-I../src -std=c++11 -g -O3 -march=native -fno-vectorize -mno-avx -mno-avx2 -mpopcnt -Rpass=loop-vectorize

LDFLAGS=-Wl -O3
# LDFLAGS=-Wl -O3 -ldw

 # CC=gcc
 # CXX=g++
 # CXX1XFLAGS=-I../src -std=c++11 -O3 -march=native -mpopcnt 
 # LDFLAGS=-Wall -O3

CXXFLAGS=$(CXX1XFLAGS)
RM=rm -f
LDLIBS= -pthread

# -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize

# SRCS=harness.cpp tools.cpp test.cpp ../src/join_base.cpp
SRCS=harness.cpp test.cpp ../src/join_base.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: bench

bench: $(OBJS)
	$(CXX) $(LDFLAGS) -o bench $(OBJS) $(LDLIBS) 

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS) bench

distclean: clean
	$(RM) *~ .depend

include .depend
