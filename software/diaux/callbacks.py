import numpy as np
def extinction_event(t, 
                      params,
                      args):
       """
       Callback for determining if an extinction has occurred and the 
       integration should be terminated. 
       """
       extinction_event.terminal = True
       _dynamics = params[:-args['num_nutrients']]
       num_params = 4 + args['num_nutrients']
       _, masses, total_mass = _unpack_masses(params, args['num_species'],
                                              args['num_nutrients'])
       masses = [np.sum(m) for m in masses]
       freqs = masses / total_mass
       if (freqs <= args['extinction_threshold']).any():
           for i, f in enumerate(freqs):
               if f <= args['extinction_threshold']:
                   args['species'][i].extinct = True
           return 0.0
       else:    
           return -1.0

def biomass_bottleneck_event(t, 
                             params,
                             args):
        
    """
    Callback for triggering a biomass bottleneck.  
    """ 
    biomass_bottleneck_event.termnal = True
    extinction_event.direction = 1
    _dynamics = params[:-args['num_nutrients']]
    num_params = 4 + args['num_nutrients']
    _, _, total_mass = _unpack_masses(params, args['num_species'], 
                                      args['num_nutrients'])
    if total_mass >= args['bottleneck_mass']:   
        return 0.0
    else:
         return 1.0

def _unpack_masses(params, num_species, num_nutrients):
    """
    Unpacks parameters in given order and returns the
    proteinaceous mass for each species and the total mass. 
    """
    _dynamics = params[:-num_nutrients]
    num_params = 4 + num_nutrients
    masses =[] 
    species = []
    total_mass = 0
    for i in range(num_species):
        _species = _dynamics[i*num_params:num_params * (i + 1)]
        _mass = _species[:-2]
        species.append(_species)
        masses.append(_mass) 
        total_mass += np.sum(_mass)
 
    return [species, masses, total_mass]