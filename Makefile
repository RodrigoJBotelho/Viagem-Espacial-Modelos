CXX ?= c++
CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra

all: bin/radiacao bin/calor1d

bin/radiacao: radiacao/radiacao.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

bin/calor1d: reentrada/calor1d.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

run: all
	@mkdir -p saida
	cd saida && ../bin/radiacao && ../bin/calor1d

exemplos: run
	cp saida/*.csv docs/exemplos/

clean:
	rm -rf bin saida

.PHONY: all run exemplos clean
