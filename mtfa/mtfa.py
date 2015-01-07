# -*- coding: utf-8 -*-
from math import pi, sqrt
from numpy import sinc, exp, inf
from scipy.integrate import quad, quadrature
import numpy as np
from scipy.optimize import newton_krylov, newton
from numpy import zeros_like,ones_like, zeros, exp
import matplotlib.pyplot as plt
np.seterr(over='ignore') #ignore overflow error 

from database import Material, materialproperty, kb, hb, q0, m0, eps0

def diff2(f, a, x, left, right=0):
    """discrete form of second order derivative d(a*df) 
    with non uniform mesh
    """
    h = x[1:] - x[:-1]
    f = np.insert(f, 0, left)
    f = np.append(f, right)

    
    d2f = (a[:-2]+a[1:-1])/(h[:-1]*(h[:-1]+h[1:]))*f[:-2] -\
          f[1:-1]/(h[:-1]+h[1:])*\
          ((a[:-2]+a[1:-1])/h[:-1]+(a[2:]+a[1:-1])/h[1:])+\
          (a[2:]+a[1:-1])/(h[1:]*(h[:-1]+h[1:]))*f[2:]
    
    return d2f


class Material(Material):
    pass
        
class Structure(object):
    """Class for structure properties.
    
    attributes:
    V0: surface potential
    T: temperate
    Nd: donor concentration in cm-3
    Na: acceptor concentration in cm-3
    length: length of structure
    n: number of points
    Eg: energy band gap
    m_eff: effective mass
    eps: dielectric constant
    Ei: ionization energy
    g: degeneracy factor 
    F: electric field in V/m
    nss: surface sheet charge in cm-2
    """
    def __init__(self, mat, V0=0, T=300, Nd=0, Na=0, length=5e-8, n=100):
        self.material = mat
        self.V0 = V0
        self.T = T
        self.Nd = Nd*1e6
        self.Na = Na*1e6
        self.length = length
        self.n = n
        self.h = length/(n+1)
        self.z = self.gen_mesh()
        self.Eg = mat.Eg(self.T)
        self.eps = self.gen_array(mat.eps)*eps0
        self.all_ionized = mat.all_ionized
        if not self.all_ionized:
            self.Ei = mat.Ei
            self.g = mat.g
        self.m_eff = mat.m_eff
        # nonparabolicity factor        
#        self.alpha = (1-self.m_eff/m0)**2/self.Eg
        self.alpha = 1/self.Eg
        self.epsrel = 0.01
        self.Ef = self.fermi()
                
    def gen_mesh(self):
        return np.linspace(0, self.length, self.n)
    
    def gen_array(self, param):
        return param*ones_like(self.z)
    
    def cdos(self, E, Ec):
        """ density of states in conduction band"""
        C = (2*self.m_eff/hb**2)**(3./2)/(2*pi**2)
        E = E - Ec
        return C*sqrt(E)*sqrt(1+self.alpha*E)*(1+2*self.alpha*E)
        
    def mcdos(self, E, Ec, z):
        """ Modified density of states in conduction band"""
        if self.Ef < 0:
            L = hb/sqrt(2*self.m_eff*kb*self.T)
        else:
            L = hb/sqrt(2*self.m_eff*self.Ef) # Fermi length
        # correction factor
        f = 1 - sinc((2*z/L)*sqrt((E-Ec)/kb/self.T)*sqrt(1+self.alpha*(E-Ec))/pi)
        return self.cdos(E, Ec)*f
    
    def fd(self, E, Ef):
        """ Fermi-Dirac function"""
        return 1/(1+exp((E-Ef)/kb/self.T))
    
    def cdos_fd(self, E, Ec, Ef):
        """multiply of two function for integral calculation"""
        return self.cdos(E, Ec)*self.fd(E, Ef)

    def mcdos_fd(self, E, Ec, z, Ef):
        """multiply of two function for integral calculation"""
        return self.mcdos(E, Ec, z)*self.fd(E, Ef)
        
    def ced(self, Ec, Ef):
        """conduction electron density"""
        return quad(self.cdos_fd, Ec, inf, args=(Ec, Ef), 
                    epsrel=self.epsrel)[0]

    def mced(self, Ec, z, Ef):
        """depth distribution of modified conduction electron density"""
        return quad(self.mcdos_fd, Ec, inf, args=(Ec, z, Ef), 
                    epsrel=self.epsrel)[0] 
                    
    def ionized(self, Ec, Ef):
        """densities of ionized shallow donors and acceptors"""
        if not self.all_ionized:
            return self.Nd/(1+self.g*np.exp((Ef-Ec+self.Ei)/kb/self.T))
        else:
            return self.Nd
            
    def charge(self, Ef):
        """equation for charge neutrality calculation"""
        return self.ionized(0, Ef) - self.ced(0, Ef)
    
    def fermi(self):
        """find fermi level using nonlinear solver"""
        return newton(self.charge, 0, maxiter=200)
    
    def poisson(self, V):
        n = zeros_like(V)
        Nd  = zeros_like(V)
        d2V = diff2(V, self.eps, self.z, self.V0)*self.h

        for i, z in enumerate(self.z[1:-1]):
            n[i] = self.mced(V[i], z, self.Ef)
            Nd[i] = self.ionized(V[i], self.Ef)
            
            
        return d2V - q0*(Nd - n)*self.h
    
    def initGuess(self):
        """ initial guess"""
        self.guess = zeros(self.n-2, float)
    
    def solve(self):
        """ solve a Poisson equation"""
        self.sol = newton_krylov(self.poisson, self.guess, 
                                 method='lgmres', verbose=1,
                                 f_tol=2e-5, maxiter=20)
        # set initial guess for next iteration in C-V
        self.guess = self.sol
        #add boundaries
        self.sol = np.insert(self.sol, 0, self.V0)
        self.sol = np.append(self.sol, 0)
        self.F = (self.sol[0]-self.sol[1])/(self.z[0]-self.z[1])
        self.nss = self.F*self.eps[0]/q0/1e4
        