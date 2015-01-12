# MTFA

MTFA is Python module for solving Poisson equation with Modified Tomas-Fermi approximation (MTFA).

## How to use it
```
import mtfa

# define structure
InN = Material('InN')
s = Structure(InN, Nd=1.9e19)

s.V0 = -0.735 #surface potential

# solve Poisson equation
s.initGuess()
s.solve()

plt.plot(s.z, s.sol) #plot solution
```

## Material properties

Properties of semiconductor materials is defined in database.py.