# MTFA

MTFA is Python module for solving Poisson equation with Modified Tomas-Fermi approximation (MTFA).

## How to use it
```python
import mtfa

# define structure
InN = Material('InN')
s = Structure(InN, Nd=1.9e19)

s.V0 = -0.735 #set surface potential

# solve Poisson equation
s.initGuess()
s.solve()

plt.plot(s.z, s.sol) #plot solution
```

## Material properties

Properties of semiconductor materials is defined in database.py.

## References

1. [Inversion and accumulation layers at InN surfaces](http://www.sciencedirect.com/science/article/pii/S0022024805014600)
TD Veal, LFJ Piper, WJ Schaff, CF McConville - Journal of crystal growth, 2006
2. [A novel self-consistent theory of the electronic structure of inversion layers in InSb MIS structures](http://onlinelibrary.wiley.com/doi/10.1002/pssb.2221340245/abstract)
J.-P. Zöllner, H. Übensee, G. Paasch, T. Fiedler andG. Gobsch - physica status solidi (b), 1986