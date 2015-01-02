# -*- coding: utf-8 -*-
from math import pi, sqrt
from numpy import sinc, exp, inf
from scipy.integrate import quad, quadrature
import numpy as np
from scipy.optimize import newton_krylov, newton
from numpy import zeros_like,ones_like, zeros, exp
import matplotlib.pyplot as plt
np.seterr(over='ignore') #ignore overflow error 

import database

#Defining constants
kb = 8.61735E-5     # Boltzmann constant in eV/K
hb = 6.58212E-16    # Plank's constant(h) in eV*s
q0 = 1.60218E-19    # electron charge (C)
m0 = 9.10938E-31    # electron mass (kg)
m0 = m0/q0          # scale mass to eV
eps0 = 8.85419E-12  # vacuum permittivity in F/m


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


class Material(object):
    """Class for material properties.
    
    attributes: 
    Eg0: energy band gap at 0K
    Eg_alpha: parameter  for Eg temperature dependence
    Eg_betta: parameter  for Eg temperature dependence
    m_e: conduction band effective mass (relative to electron mass)
    m_hh: heavy hole band effective mass 
    m_lh: light hole band effective mass
    m_eff: effective mass
    eps: dielectric constant
    Ei: ionization energy
    g: degeneracy factor
    dop_type: type of dopant (donor or acceptor)
    """
    def __init__(self, material, dopant, database=database):
        self.material = material
        self.dopant = dopant
        self.material_property = database.materialproperty
        matprops = self.material_property[material]
        self.Eg0 =  matprops['Eg0']
        self.Eg_alpha =  matprops['Eg_alpha']
        self.Eg_betta =  matprops['Eg_betta']
        self.m_e = matprops['m_e']*m0
        self.m_hh = matprops['m_hh']*m0
        self.m_lh = matprops['m_lh']*m0
        self.eps = matprops['eps']
        self.Ei = matprops['impurity'][dopant]['Ei']
        self.g = matprops['impurity'][dopant]['g']
        self.dop_type = matprops['impurity'][dopant]['type']
        if self.dop_type == 'donor':
            self.m_eff = self.m_e
    
    def Eg(self, T):
        """Temperature dependence of energy bandgap"""
        if self.Eg0 == None:
            return self.Eg
        else:
            return self.Eg0 - self.Eg_alpha*T**2/(T + self.Eg_betta) 
        
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
        self.Ei = mat.Ei
        self.g = mat.g
        self.m_eff = mat.m_eff
        self.Ef = self.fermi()
    
    def gen_mesh(self):
        return np.linspace(0, self.length, self.n)
    
    def gen_array(self, param):
        return param*ones_like(self.z)
    
    def cdos(self, E, Ec, z):
        """ Modified density of states in conduction band"""
        C = (2*self.m_eff/(hb**2))**(3./2)/(2*pi**2)
        L = hb/sqrt(2*self.m_eff*kb*self.T) # Fermi length
        # nonparabolicity factor        
        alpha = (1-self.m_eff/m0)**2/self.Eg 
        E = E - Ec
        # correction factor
        f = 1 - sinc((2*z/L)*sqrt(E/kb/self.T)*sqrt(1+alpha*E)/pi)
        return C*sqrt(E)*sqrt(1+alpha*E)*(1+2*alpha*E)*f
    
    def fd(self, E, Ef):
        """ Fermi-Dirac function"""
        return 1/(1+exp((E-Ef)/kb/self.T))
    
    def cdos_fd(self, E, Ec, z, Ef):
        """multiply of two function for integral calculation"""
        return self.cdos(E, Ec, z)*self.fd(E, Ef)
    
    def ced(self, Ec, z, Ef):
        """depth distribution of conduction electron density"""
        return quad(self.cdos_fd, Ec, inf, args=(Ec, z, Ef), 
                    epsrel=0.01)[0]
    
    def ionized(self, Ec, Ef):
        """densities of ionized shallow donors and acceptors"""
        if self.material.dop_type == 'donor':
            return self.Nd/(1+self.g*np.exp((Ef-Ec+self.Ei)/kb/self.T))
    
    def charge(self, Ef):
        """equation for charge neutrality calculation"""
        return self.ionized(0, Ef) - self.ced(0, 1e-3, Ef)
    
    def fermi(self):
        """find fermi level using nonlinear solver"""
        return newton(self.charge, 0, maxiter=200)
    
    def poisson(self, V):
        n = zeros_like(V)
        Nd  = zeros_like(V)
        d2V = diff2(V, self.eps, self.z, self.V0)*self.h

        for i, z in enumerate(self.z[1:-1]):
            n[i] = self.ced(V[i], z, self.Ef)
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
        self.guess = self.sol
        #add boundaries
        self.sol = np.insert(self.sol, 0, self.V0)
        self.sol = np.append(self.sol, 0)
        # set initial guess for next iteration in C-V
        