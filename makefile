# Définition des variables
FC = gfortran
FFLAGS = -O2 -Wall -g -fcheck=all -fbacktrace -ffpe-trap=zero,overflow,underflow
OBJ = main.o fonction.o
EXEC = program

# Règle pour compiler l'exécutable
$(EXEC): $(OBJ)
	$(FC) $(OBJ) -o $(EXEC)

# Règle pour compiler main.o, en s'assurant que fonction.o est compilé avant
main.o: main.f90 fonction.o
	$(FC) $(FFLAGS) -c main.f90

# Règle pour compiler fonction.o et générer fonction.mod
fonction.o: fonction.f90
	$(FC) $(FFLAGS) -c fonction.f90

# Règle pour nettoyer les fichiers objets et l'exécutable
clean:
	rm -f $(OBJ) $(EXEC) *.mod
