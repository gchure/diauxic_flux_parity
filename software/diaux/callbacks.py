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
       masses = np.zeros(args['num_species'])
       total_mass = 0
       for i in range(args['num_species']):
           _mass = _dynamics[i*num_params:num_params * (i + 1)][:-2]    
           masses[i] = np.sum(_mass)
           total_mass += masses[i] 
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
    masses = np.zeros(args['num_species'])
    total_mass = 0
    for i in range(args['num_species']):
           _mass = _dynamics[i*num_params:num_params * (i + 1)][:-2]    
           masses[i] = np.sum(_mass)
           total_mass += masses[i] 
    if total_mass >= args['bottleneck_mass']:   
        return 0.0
    else:
         return 1.0
