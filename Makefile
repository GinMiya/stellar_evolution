CPP := g++
OPTION := -O3 -ffast-math -funroll-loops -march=native -fopenmp -mavx

SOURCE := dens_temp_mass.cpp
FILE := dens_temp_mass
FILE:
	$(CPP) $(SOURCE) -o $(FILE) $(OPTION)

clean:
	rm -f $(FILE)
	rm -r -f data_HHeZ
	rm -r -f data_CO
	rm -r -f data_tmax_HHeZ
	rm -r -f data_tmax_CO
	rm -r -f data_mmin
