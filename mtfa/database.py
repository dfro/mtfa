# -*- coding: utf-8 -*-
# Ref. http://www.ioffe.ru/SVA/NSM/Semicond
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
'm_e':1.18, #  effective mass of the density of states with for 6 valleys
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