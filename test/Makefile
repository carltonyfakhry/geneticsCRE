CC=clang
CXX=clang++
RM=rm -f
CXX1XFLAGS=-std=c++11 -g -O3 -march=native -mno-avx -mno-avx2 -mpopcnt -I../src -Rpass=loop-vectorize
CXXFLAGS=$(CXX1XFLAGS)
# LDFLAGS=-g,-Wl,-O3,--sort-common,--as-needed,-z,relro
LDFLAGS=-g -Wl -O3
LDLIBS=

# -Rpass=loop-vectorize -Rpass-missed=loop-vectorize 

SRCS=harness.cpp $(wildcard ../src/method*.cpp)
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
