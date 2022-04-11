gcc_options = -std=c++17 -Wall --pedantic-errors -DMKL_ILP64  -I"${MKLROOT}/include" -g
l_b = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

program : main.o make_hamiltonian.o spm.o smp.o szz.o
	g++ -o $@ $^ $(l_b)

main.o : main.cpp
	g++ -c $(gcc_options) $< $(l_b)


make_hamiltonian.o : ./make_hamiltonian/make_hamiltonian.cpp
	g++ -c $(gcc_options) $< $(l_b)

spm.o : ./make_hamiltonian/spm.cpp
	g++ -c $(gcc_options) $< $(l_b)

smp.o : ./make_hamiltonian/smp.cpp
	g++ -c $(gcc_options) $< $(l_b)

szz.o : ./make_hamiltonian/szz.cpp
	g++ -c $(gcc_options) $< $(l_b)


run : program
	./program

clean:
	rm -f ./program

.PHONY : run clean