C = gcc
OPTS = -c -fPIC
LOADER = gcc

OBJECTS = bayesian_probit_regression_cfunctions.c

 bayesian_probit_regression_cfunctions.so: $(OBJECTS)
	R CMD SHLIB -o bayesian_probit_regression_cfunctions.so $(OBJECTS)

.c.o: ; $(C) $(OPTS) -c $<

clean:
	-rm *.o *.so
