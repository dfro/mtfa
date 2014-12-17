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
        'Ea':0.001, # (eV) ionization energies
        'g':2, # degeneracy factor
        'type':'donor'
        },
    },

},
}