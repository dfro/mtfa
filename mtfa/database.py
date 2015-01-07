# -*- coding: utf-8 -*-
# Ref. http://www.ioffe.ru/SVA/NSM/Semicond 
# for material properties

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
    Ei: ionization energy
    g: degeneracy factor
    dop_type: type of dopant (donor or acceptor)
    """
    def __init__(self, material, dopant=None):
        self.material = material
        self.dopant = dopant
        self.material_property = materialproperty
        matprops = self.material_property[material]
        self.Eg_300 = matprops['Eg']
        self.Eg0 =  matprops['Eg0']
        self.Eg_alpha =  matprops['Eg_alpha']
        self.Eg_betta =  matprops['Eg_betta']
        self.m_e = matprops['m_e']*m0
        self.m_hh = matprops['m_hh']*m0
        self.m_lh = matprops['m_lh']*m0
        self.eps = matprops['eps']
        if dopant == None:
            self.all_ionized = True
        else:
            self.all_ionized = False
            self.Ei = matprops['impurity'][dopant]['Ei']
            self.g = matprops['impurity'][dopant]['g']
            self.dop_type = matprops['impurity'][dopant]['type']
        
        self.m_eff = self.m_e
    
    def Eg(self, T):
        """Temperature dependence of energy bandgap"""
        if self.Eg0 == None:
            return self.Eg_300
        else:
            return self.Eg0 - self.Eg_alpha*T**2/(T + self.Eg_betta)

materialproperty = {
'InAs':{
'Eg':0.354, # (eV) at 300K
'Eg0':0.415, # (eV) band gap at 0K
'Eg_alpha':2.76e-4, # (eV/K) parameter  for Eg temperature dependence
'Eg_betta':83, # (K) parameter  for Eg temperature dependence
'm_e':0.023, #conduction band effective mass (relative to electron mass)
'm_hh':0.41, #heavy hole band effective mass 
'm_lh':0.026, #light hole band effective mass
'eps':15.15, #dielectric constant
'impurity':{ #impurity properties
    'S':{
        'Ei':0.001, # (eV) ionization energies
        'g':2, # degeneracy factor
        'type':'donor'
        },
    },
},

'InN':{
'Eg':0.503, #0.642
'Eg0':None, #
'Eg_alpha':2.5e-4, 
'Eg_betta':624, 
'm_e':0.045, #0.11 
'm_hh':1.63, 
'm_lh':0.27,
'eps':15.3, 
'impurity':{ 
    'X':{
        'Ei':0.01,
        'g':2, 
        'type':'donor'
        },
    },
},

'Si':{
'Eg':1.12, 
'Eg0':1.17,
'Eg_alpha':4.73e-4, 
'Eg_betta':636, 
'm_e':1.18, #  effective mass of the density of states for 6 valleys
'm_hh':0.49, 
'm_lh':0.16,
'eps':11.7, 
'impurity':{ 
    'P':{
        'Ei':0.1, #0.045
        'g':2, 
        'type':'donor'
        },
    'B':{
        'Ei':0.045, 
        'g':4, 
        'type':'acceptor'
        },
    },
},

'SnO2':{
'Eg':3.6, 
'Eg0':None,
'Eg_alpha':None, 
'Eg_betta':None, 
'm_e':0.27,
'm_hh':0, 
'm_lh':0,
'eps':12.2, 
'impurity':{ 
    'Sb':{
        'Ei':10e-3, 
        'g':2, 
        'type':'donor'
        },
    },
},

}