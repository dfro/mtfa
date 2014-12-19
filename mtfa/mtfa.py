# -*- coding: utf-8 -*-
from math import pi, sqrt
from numpy import sinc, exp, inf
from scipy.integrate import quad, quadrature
import numpy as np
from scipy.optimize import newton_krylov
from numpy import zeros_like, mgrid, zeros, exp
import matplotlib.pyplot as plt

import database

#Defining constants
kb = 8.61735E-5     # Boltzmann constant in eV/K
hb = 6.58212E-16    # Plank's constant(h) in eV*s
q0 = 1.60218E-19    # electron charge (C)
m0 = 9.10938E-31    # electron mass (kg)
m0 = m0/q0          # scale mass to eV
eps0 = 8.85419E-12  # vacuum permittivity in F/m

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
    Ea: ionization energy
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
        self.Ea = matprops['impurity'][dopant]['Ea']
        self.g = matprops['impurity'][dopant]['g']
        self.dop_type = matprops['impurity'][dopant]['type']
        if self.dop_type == 'donor':
            self.m_eff = self.m_e
    
    def Eg(self, T):
        """Temperature dependence of energy bandgap"""
        return self.Eg0 - self.Eg_alpha*T**2/(T + self.Eg_betta) 
        
class Structure(object):
    """Class for structure properties.
    
    attributes:
    T: temperate
    Nd: donor concentration
    Na: acceptor concentration
    Eg: energy band gap
    m_eff: effective mass
    eps: dielectric constant
    Ea: ionization energy
    g: degeneracy factor    
    """
    def __init__(self, mat, T=300, Nd=0, Na=0):
        self.material = mat
        self.T = T
        self.Nd = Nd
        self.Na = Na
        self.Eg = mat.Eg(self.T)
        self.eps = mat.eps
        self.Ea = mat.Ea
        self.g = mat.g
        self.m_eff = mat.m_eff
     
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
    
    def ced(self, Ec, z, Ef):
        """depth distribution of  conduction electron density"""
        