import numpy as np
import scipy.integrate
import pandas as pd
from .callbacks import extinction_event, _unpack_masses
OD_CONV = 1.15E17 # Amino acids per OD600 unit (approximately)

#TODO: Set `reset_properties` function to clear and reninitialize all properties
class FluxParityAllocator:
    """Base class for a self replicator obeying flux-parity allocation."""
    def __init__(self, 
                 suballocation, 
                 constants={}, 
                 label=0):
        """
        Instantiates a self replicating organism that undergoes flux-parity regulation 
        of its resource allocation. 

        Parameters
        ----------
        suballocation : dict 
            A dictionary defining the metabolic suballocation strategy of the 
            self replicator. Must have the following keys:
                strategy: str, `'dynamic'` or `'static'`
                    The type of suballocation. If `'dynamic'`, the suballocation
                    of each metabolic sector follows a Monod function with a 
                    Monod constant `K` and sensitivty `n`. If `'static'`,
                    suballocation is deemed to be at a fixed value, supplied
                    with key `alpha.`    
                nu_max: numpy.ndarray of floats [0, inf)
                    The metabolic rate for consumption of each nutrient in units 
                    of inverse hours.
                frac_useful: numpy.ndarray
                    The fraction of each metabolic class that is deemed to be 
                    useful and contributes to the net metabolic flux.
                    Default is assumed to be 1 (all is useful.)
                hierarchy : numpy.ndarray of ints
                    The hierarchy order of the nutrients beginning at 0. This is only used if 
                    `strategy: 'dynamic'` is supplied. If not supplied, hierarchy 
                    is assumed to be in the order of the provided nutrients.
                alpha : `numpy.ndarray` of float, [0, 1]
                    The fixed suballocation of the metabolic strategy. Values 
                    must sum to `1.0`. This is only used if `strategy: 'static'` 
                    is supplied.
                K : numpy.ndarray of float [0, inf)
                    The Monod constant for the suballocation, in units of concentration.
                    This is only used if `strategy: 'dynamic'` is supplied.
                n : numpy.ndarray of int [1, inf) 
                    The sensitivity parameter of the Monod function. This is only 
                    used if `strategy: 'dynamic'` is supplied.

               
        constants : dict, optional
            A dictionary of the constant model parameters. Default values
            (described below) can be overridden by providing a key,value pair.
                gamma_max : float
                    The maximum translation rate in inverse hours. Default 
                    value is 9.65 inv. hr. 
                phi_O : float [0, 1]
                    The allocation towards all non-ribosomal and non-metabolic 
                    proteins. Default value is 0.55. 
               tau : float [0, inf)                               
                    The sensitivity parameter of the flux-parity ribosomal 
                    allocation function. Default is 1.0
               kappa_max : float [0, inf)
                    The maximum transcription rate of uncharged tRNA in 
                    mass abundance units per hour. Default is 1.15E-3 inverse hr.
               Km_u : float [0, inf)                
                    The Michaelis-Menten constant for binding of uncharged tRNA 
                    to the metabolic components in units of relative mass abundance. 
                    Default value is 3E-5.
               Km_c : float [0, inf)                
                    The Michaelis-Menten constant for binding of charged tRNA 
                    to the ribosomal components in units of relative mass abundance. 
                    Default value is 3E-5.
               Km : numpy.ndarray [0, inf)
                    The Monod constant for metabolic activity as a function of
                    the external nutrient concentration in units of concentration. 
                    Default is the same for each nutrient with 5E-6 M.
               Y : numpy.ndarray [0, inf)
                    The yield coefficient for growth on each nutrient. Default 
                    is 2.95E19 for each nutrient supplied.              
               death_rate : float
                    The death rate of the species in units of mass per unit time 
                    (typically hours). Default is 0.

        """
        self.num_metab = len(suballocation['nu_max'])
        self.label = label
        self.extinct = False
        self._properties = False
        # Set the properties of the self replicator.
        _constants = {'gamma_max': 9.65,
                     'phi_O': 0.55,
                     'tau': 1,
                     'kappa_max':1.15E-3,
                     'Km_u': 3E-5,
                     'Km_c': 3E-5,
                     'Km': [5E-6 for _ in range(self.num_metab)],
                     'Y': [2.95E19 for _ in range(self.num_metab)],
                     'nu_max': suballocation['nu_max'],
                     'death_rate': 0}

        self._overridden_pars = {}
        for k, v in constants.items():
            _constants[k] = v
            self._overridden_pars[k] = v
        if 'frac_useful' not in suballocation:
            self.frac_useful = np.array([1 for _ in self.num_metab])
        else:
            self.frac_useful = np.array(suballocation['frac_useful'])
        # Set attributes from the constant dictionary.
        for k, v in _constants.items():
            setattr(self, k, v)

        if 'hierarchy' not in suballocation.keys():
            self.hierarchy = np.arange(self.num_metab).astype(int)
        else:
            self.hierarchy = np.array([int(k) for k in suballocation['hierarchy']]).astype(int)

        # Set the suballocation details
        self.strategy = suballocation['strategy']
        if self.strategy == 'static':
            if np.sum(suballocation['alpha']) != 1:
                raise ValueError(f"Suballocation parameters must sum to 1!")
            self.alpha = suballocation['alpha'] 
        elif self.strategy == 'dynamic':
            self.K = suballocation['K']
            self.n = suballocation['n']
        else:
            raise ValueError("Supplied metabolic strategy must be either `static` or `dynamic`.") 

    def __repr__(self):
        rep = f"""
=============== Self Replicator #{self.label} ================
extinct?                    : {self.extinct}
metabolic_strategy          : {self.strategy}
number of metabolic classes : {self.num_metab}
metabolic rates             : {self.nu_max}
hierarchy                   : {self.hierarchy}
"""
        if self.strategy == 'dynamic':
            rep += f"\tK's [M] : {self.K}"
            rep += f"\n\tn's : {self.n}"
        elif self.strategy == 'static':
            rep += f"\talpha's: {self.alpha}"
        if len(self._overridden_pars) > 0:
            rep += f"\nOverridden Default Parameters:" 
            for k, v in self._overridden_pars.items():
                rep += f"\n\t{k}: {v}"
        if self._properties:
            rep+=f"""\n
---------------------Allocation------------------------------
tRNA ratio (charged/uncharged) : {self.ratio}
gamma (translation rate)       : {self.gamma} [inv. hr.]
nu (metabolic rate)            : {self.nu} [inv. hr.]
kappa (transcription rate)     : {self.kappa} [inv. hr.]
phi_Rb (ribosomal allocation)  : {self.phi_Rb}   
phi_Mb (metabolic allocation)  : {self.phi_Mb} 
            """ 
        return rep 

    def compute_properties(self, tRNA_c, tRNA_u, nutrients):
        """
        Computes the self replicator properties, including the allocation parameters 
        and the corresponding rates. 

        # Parameters
        # ----------
        # tRNA_c : float [0, inf)
        #     The charged tRNA concentration in relative mass abundance units.
        # tRNA_u : float[0, inf) 
        #     The uncharged tRNA concentration in relative mass abundance units.
        # nutrients : numpy.ndarray float 
        #     The concentration of the nutrients in the environment for calculation
        #     of rates. 
        """
        self.tRNA_u = tRNA_u
        self.tRNA_c = tRNA_c
        if tRNA_u == 0: #Include less than zero for numerical underflow
            if tRNA_c == 0:
                self.phi_Rb = 0
                self.ratio = 0
            else:
                self.phi_Rb = (1 - self.phi_O)
                self.ratio = np.inf
        else:
            self.ratio = tRNA_c / tRNA_u
            self.phi_Rb = (1 - self.phi_O) * (self.ratio / (self.ratio + self.tau))

        # Adjust the flux parity regulation on transcription
        self.kappa = self.kappa_max * self.phi_Rb / (1 - self.phi_O)

        # Set the translation rate
        self.gamma = self.gamma_max * tRNA_c / (tRNA_c + self.Km_c)
        
        # Set the metabolic rates
        self.nu  = np.zeros_like(self.nu_max)
        metab_factor = tRNA_u / (tRNA_u + self.Km_u)
        for i in range(self.num_metab):
            env_factor = nutrients[i] / (nutrients[i] + self.Km[i])
            self.nu[i] = self.nu_max[i] * env_factor * metab_factor

        # Compute the metabolic allocation
        self.phi_Mb = np.zeros(self.num_metab)
        if self.strategy == 'static':
            # Set the fixed suballocations as prescribed by the alpha parameters.
            for i in range(self.num_metab):
                self.phi_Mb[i] = (1 - self.phi_O - self.phi_Rb) * self.alpha[i] 

        elif self.strategy == 'dynamic':
            self.alpha = np.zeros(self.num_metab)
            # set the suballocation based on the nutrient concentrations
            occupied_phi_Mb = 0 
            occupied_suballocation = 0
            idx_sort = self.hierarchy.argsort()
            for i, v in enumerate(idx_sort):  
                if self.hierarchy[v] == (self.num_metab - 1):
                    # Evaluate the remaining suballocation to the final nutrient in the hierarchy.
                    self.alpha[i] = (1 - occupied_suballocation)
                    self.phi_Mb[i] = (1 - self.phi_O - self.phi_Rb - occupied_phi_Mb) 
                else:
                    # Set the suballocation given the supplied monod function properties
                    numer = (nutrients[v] / self.K[i])**self.n[v]
                    factor = numer / (numer + 1)
                    self.alpha[i] = (1 - occupied_suballocation) * factor
                    self.phi_Mb[v] =  (1 - self.phi_O - self.phi_Rb - occupied_phi_Mb) * factor
                    occupied_suballocation += self.alpha[i]
                    occupied_phi_Mb += self.phi_Mb[i]
        self._properties = True

    def compute_derivatives(self, 
                            masses, 
                            nutrients):
        """
        Computes the mass and concentration derivatives given the internal flux-parity 
        allocation state.

        Parameters
        ----------
        masses : numpy.ndarray float
            The masses of the protein sectors and the concentrations of the charged 
            and uncharged tRNAs. This should be in the order of the metabolic masses,
            ribosomal mass, "other" protein masses, uncharged tRNAs, and charged tRNAs

        nutrients: numpy.ndarray float 
            The concentrations of the external nutrients.

        Returns
        -------
        masses_dot : numpy.ndarray float 
            The time derivative of the masses given the self replicator dynamics.
        nutrients_dt : numpy.ndarray float
            The time derivative the nutrient concentrations consumed by a the 
            self replicator given the dynamics.    
        """

        # Unpack the variables 
        self.M_Mb = masses[:self.num_metab]
        self.M_Rb = masses[self.num_metab]
        self.M_O = masses[self.num_metab + 1]
        self.M = np.sum(masses[:self.num_metab + 2])
        self.tRNA_u = masses[self.num_metab + 2]
        self.tRNA_c = masses[self.num_metab + 3]
        self.nutrients = np.array(nutrients) * (np.array(nutrients) >= 0)

        # Enforce positivity of all variables. 
        self.M_Mb *= self.M_Mb > 0
        self.M_Rb *= self.M_Rb > 0
        self.M_O *= self.M_O > 0
        self.tRNA_u *= self.tRNA_u #> self.min_tRNA_u
        self.tRNA_c *= self.tRNA_c #> self.min_tRNA_c

        # Evaluate the properties
        self.compute_properties(self.tRNA_c, self.tRNA_u, self.nutrients)

        # Evalutate the derivatives
        dM_Rb = self.phi_Rb * self.gamma * self.M_Rb - self.M_Rb * self.death_rate
        dM_O = self.phi_O * self.gamma * self.M_Rb  -  self.M_O * self.death_rate
        dM_Mb = [phi_Mb * self.gamma * self.M_Rb - self.M_Mb[i] * self.death_rate for i, phi_Mb in enumerate(self.phi_Mb)]
        dtRNA_u = self.kappa + self.gamma * self.M_Rb * (1 - self.tRNA_u) / self.M - np.sum(self.nu * self.M_Mb * self.frac_useful) / self.M - self.tRNA_u * self.death_rate
        dtRNA_c = np.sum(self.nu * self.M_Mb * self.frac_useful) / self.M - self.gamma * self.M_Rb * (1 + self.tRNA_c) / self.M - self.tRNA_c * self.death_rate
        dc_nt = -self.nu * self.M_Mb / self.Y
        masses_dt = []
        for d in dM_Mb:
            masses_dt.append(d)
        masses_dt.append(dM_Rb)
        masses_dt.append(dM_O)
        masses_dt.append(dtRNA_u)
        masses_dt.append(dtRNA_c)
        return [np.array(masses_dt), dc_nt]


class Ecosystem:
    """Base class for an ecosystem of self-replicators"""
    def __init__(self, 
                 species, 
                 nutrients):
        """
        Set up an ecosystem of J self replicators growing on K nutrients

        Parameters
        ----------
        species : list of `FluxParityAllocator`
            A list of the consitituent species of the ecosystem. Each species
            should be specified as a FluxParityAllocator object.

        nutrients: list of dict 
            A list of nutrients with dictionary specifying the details about the 
            nutrient loads with the following keys.
                init_concs : numpy.ndarray
                    The initial concentrations of the nutrients in the environment.
                feed_concs : numpy.ndarray
                    The concentration of the nutrients in the feed stock. If not 
                    provided, it will be assumed to be the same as the initial 
                    concentration.
                inflow_rates : numpy.ndarray 
                    The inflow rates of each nutrient in units of concentration 
                    per unit time.  If not supplied, it will be assumed to be 0.
                degradation_rates : numpy.ndarray
                    The degradation rates of each nutrient in units of concentraiton 
                    per unit time. If not supplied, it will be assumed to be 0.
 
        """
        self.species = species
        self.num_species =  len(species)
        self.num_nutrients = len(nutrients['init_concs'])
        self.extinction = False
        self._seed = False

        # Set the defaults
        if 'feed_concs' not in nutrients.keys():
            nutrients['feed_concs'] = nutrients['init_concs']
        for k in ['inflow_rates', 'degradation_rates']:
            if k not in nutrients.keys():
                nutrients[k] = np.zeros(self.num_nutrients) 
        for k, v in nutrients.items():
            setattr(self, k, v)

    # TODO: Allow for approximate steady state initiation. 
    def seed(self, 
             freqs=None, 
             steadystate=True,
             od_init=0.04,
             max_iter=1000):
        """
        Initialize the ecosystem by setting the starting masses and tRNA 
        concentrations with a prescribed starting frequency. 

        Parameters
        ----------
        freqs : numpy.ndarray or None
            The desired starting frequencies for each species. If `None`, all 
            species are assumed to start with equal frequency.
        steadystate: bool
            If `True`, the system is initialized in the steay-state regime 
            with approximately fixed nutrient concentrations as provided by 
            `init_concs`.  Steady-state is defined when dtRNA_u and dtRNA_c 
            are within 1E-9 units of 0 and phi_Rb is equivalent to the 
            ribosomal content (M_Rb/M) to within 1E-3.
        od_init : float
            The approximate initial optical density of the culture. This is 
            converted to mass of amino acids using a conversion factor of 
            1.15E17 AA per OD600. Default value is 0.04.
        max_iter : int
            The maximum number of iterations initialization undergoes to 
            approximate the steady state configuration. This only matters if 
            `steadystate = True`.
        """
        if freqs is None:
            freqs = np.ones(self.num_species) / self.num_species
        M0 = od_init * OD_CONV * freqs
        params = []
        _num_species = self.num_species
        for i, _species in enumerate(self.species):
            # Set the concentrations of the charged and uncharged tRNA and 
            # compute the initial properties
            _species.compute_properties(1E-5, 1E-5, self.init_concs)
            # If the approximate initial steady-state is desired, approximate
            if steadystate:
                # Reset the number of species as a convenience while equilibrating 
                # to steady state. This gets reset back to the original value upon
                # completion of the loop. 
                self.num_species = 1
                self.extinction_threshold = None
                t_range = [0, 1]
                print(f'Finding the the approximate steady state for species {_species.label}...', end='') 
                ss = False
                for j in range(max_iter):
                        p0 = [phi_Mb for _, phi_Mb in enumerate(_species.phi_Mb)]
                        p0.append(_species.phi_Rb)
                        p0.append(_species.phi_O)
                        p0.append(_species.tRNA_u)
                        p0.append(_species.tRNA_c)
                        for c in self.init_concs:
                            p0.append(c)

                        _, _ = self._integrate(t_range, p0)
                        last_soln = self.last_soln.y[:, -1]
                        # Compute the derivatives for each species
                        d_masses, _ = _species.compute_derivatives(last_soln[:self.num_nutrients + 4], last_soln[-self.num_nutrients:])
                        dtRNA_u, dtRNA_c = d_masses[-2:]
                        dtRNA_u *= np.abs(dtRNA_u) <= 1E-9
                        dtRNA_c *= np.abs(dtRNA_c) <= 1E-9
                        total_mass = np.sum(self.last_soln.y[:-1][:-(self.num_nutrients + 2)])
                        rb_content_diff = np.abs(_species.phi_Rb - (self.last_soln.y[:,-1][self.num_nutrients+1]/total_mass))
                        if (dtRNA_c == 0) & (dtRNA_u == 0) & (rb_content_diff <= 1E-3):
                            print('done!')
                            ss = True
                            break
                        
                if ss == False:
                    print('FAILED. No steady-state found. Increase `max_iter`. Proceeding using last state.')                        
                delattr(self, "extinction_threshold")
                delattr(self, "last_soln")
            # Set the initial parameters in the order of metabolic proteins,
            # ribosomal proteins, other proteins, uncharged tRNA, and charged tRNA
            for _, phi_Mb in enumerate(_species.phi_Mb):
                params.append(phi_Mb * M0[i])
            params.append(_species.phi_Rb * M0[i])
            params.append(_species.phi_O * M0[i])
            if steadystate:
                params.append(_species.tRNA_u)
                params.append(_species.tRNA_c)
            else:
                params.append(1E-5)
                params.append(1E-5)
        for c in self.init_concs:
            params.append(c)  
        self.num_species = _num_species
        self._seed = params 


    def _dynamical_system(self, t, 
                          params, 
                          args):
        """
        Computes the derivatives of the masses of the community members.
        """ 

        # Re-access dimensional information
        num_params = 4 + self.num_nutrients 
        nutrients = params[-self.num_nutrients:]  

        # Re-access ecosystem parameters
        inflow = self.inflow_rates
        deg = self.degradation_rates
        feed = self.feed_concs

        # Unpack the parameters for easy iteration
        if self.num_species > 1:
            masses = [params[i*num_params:i*num_params + num_params] for i in range(self.num_species)]
        else:
            masses = [params[:-self.num_nutrients]]

        # Initialize storage vectors for the derivatives
        derivs = []
        nutrient_derivs = np.zeros_like(nutrients)
        for i, m in enumerate(masses):
            species = args['species'][i]
            _masses, _nutrients = species.compute_derivatives(m, nutrients)
            for _m in _masses:
                derivs.append(_m) 
            nutrient_derivs += _nutrients 
    
        for i, n in enumerate(nutrient_derivs):
            # Account for chemostatic change
            _n = inflow[i] * feed[i] + n - deg[i] * nutrients[i] 
            derivs.append(_n)
        return derivs 
    
    def _parse_soln(self, 
                    soln, 
                    tol=1E-10, 
                    tshift=0):
        """
        Unpacks the Bunch object into two dataframes -- one for the species 
        abundances and one for the nutrients
        """
        num_params  = self.num_nutrients + 4 

        # Parse the solution object for the nutrient dynamics and pack into a 
        # pandas DataFrame
        nutrient_df = pd.DataFrame([])
        nutrient_dynamics = soln.y[-self.num_nutrients:]
        nutrient_dynamics *= np.abs(nutrient_dynamics) >= tol
        species_dynamics = soln.y[:-self.num_nutrients]
        for i in range(self.num_nutrients):
            _nutrient = nutrient_dynamics[i, :]
            _df = pd.DataFrame(_nutrient.T, columns=['env_conc'])
            _df['time_hr'] = soln.t + tshift
            _df['feed_conc'] = self.feed_concs[i]
            _df['inflow_rate'] = self.inflow_rates[i]
            _df['degradation_rate'] = self.degradation_rates[i]
            _df['nutrient_label'] = i+1
            nutrient_df = pd.concat([nutrient_df, _df])

        # Compute the total mass of the ecosystem to calculate the frequencies.
        total_mass = np.zeros_like(soln.t)
        for i in range(self.num_species):
            _species = species_dynamics[i*num_params:i*num_params + 4]
            total_mass += np.sum(_species, axis=0)
        self.total_mass = total_mass
        # Parse the solution result for the masses of the individual species.
        species_df = pd.DataFrame([])
        colnames = [f'M_Mb_{i+1}' for i in range(self.num_nutrients)]
        for _c in ['M_Rb', 'M_O', 'tRNA_u', 'tRNA_c']:
            colnames.append(_c)
        for i in range(self.num_species):
            FPA = self.species[i]
            # Chunk the result into the individual species
            _species = species_dynamics[i*num_params:num_params * (i + 1)]
        
            # Convert small negative values to 0 with the given tolerance
            _species *= np.abs(_species) >= tol

            # Pack into a dataframe and label
            _df = pd.DataFrame(np.abs(_species.T), columns=colnames)
            _df['M'] = np.sum(_species[:-2], axis=0) # Computed total mass
            _df['frequency'] = _df['M'] / total_mass
            _df['ribosome_content'] = _df['M_Rb'] / _df['M']

            # Set up storage components for tracking the allocation and flux-parity 
            # details for each time point
            alloc = {'phi_Rb':[], 
                     'kappa':[], 
                     'gamma':[]}
            for j in range(self.num_nutrients):
                alloc[f'phi_Mb_{j+1}'] = []
                alloc[f'alpha_{j+1}'] = []
                alloc[f'nu_{j+1}'] = []

            # Iterate through each time point to calculate the flux-parity 
            # details
            for j in range(len(soln.t)):
                FPA.compute_properties(_species[-1, j], _species[-2, j], nutrient_dynamics[:, j])
                alloc['phi_Rb'].append(FPA.phi_Rb)
                alloc['kappa'].append(FPA.kappa)
                alloc['gamma'].append(FPA.gamma)
                for k in range(self.num_nutrients):
                    alloc[f'phi_Mb_{k+1}'].append(FPA.phi_Mb[k])
                    alloc[f'alpha_{k+1}'].append(FPA.alpha[k])
                    alloc[f'nu_{k+1}'].append(FPA.nu[k])

            # Update the dataframe with the computed properties
            for k, v in alloc.items():
                _df[k] = v

            # Add auxilliary information and store
            _df['time_hr'] = soln.t + tshift
            _df['death_rate'] = self.species[i].death_rate
            _df['species_label'] = self.species[i].label
            species_df = pd.concat([species_df, _df])

        return [species_df, nutrient_df]
 
    def _integrate(self, 
                   time_range, 
                   p0,
                   tol=1E-10, 
                   solver_kwargs={},
                   tshift=0):
        """
        Integrate the temporal dynamics of the ecosystem. 
        """
        args = {'species': self.species}
        events = []
        if self.extinction_threshold is not None: 
            events.append(extinction_event)
            args['num_nutrients'] = self.num_nutrients
            args['num_species'] = self.num_species
            args['extinction_threshold'] = self.extinction_threshold

        if len(events) > 0:
            soln = scipy.integrate.solve_ivp(self._dynamical_system, time_range, p0,
                                             args=(args,), events=events, **solver_kwargs)
            if soln.status == 1:
                self.extinction = True
            else:
                self.extinction = False
        else:
            soln = scipy.integrate.solve_ivp(self._dynamical_system, time_range, p0,
                                             args=(args,), **solver_kwargs)

        self.last_soln = soln
        # Parse the output and return the dataframes
        species_df,  nutrient_df = self._parse_soln(soln, tol, tshift=tshift)
        return species_df, nutrient_df

    def grow(self,
             time, 
             extinction_thresh=None, 
             bottleneck = {}, 
             dt = 0.01,
             tol=1E-10, 
             solver_kwargs={}):
        """
        Grows the ecosystem with a seeded community and returns the 
        dynamics as Pandas DataFrames.

        Parameters
        ----------
        time : float
            The total time for the integration to take place. This should be 
            given in units of hours. 
        extinction_thresh : float or None
            The threshold frequency at which a species extinction takes place. 
            Default is 1E-3 (0.1% by mass) which terminates the integration. If
            None, the extinction callback is not applied.
         bottleneck : dict of dicts
            A dictionary specifying the types of bottleneck that should take place 
            during the integration. 
            
            
            # The following callbacks `type`'s are supported. 
            # Each carries with it another dictionary of `args`:
            #     'time_dilution' : dict
            #         A time interval by which the community should be diluted. 
            #         Must provide the time 'interval' where the dilution takes 
            #         place. Must also provide *either* `factor`, which prescribes
            #         the magnitude of the dilution (e.g. 0.1 is a 10 fold dilution),
            #         or `target`, which is the target approximate optical density 
            #         the culture is diluted to.
            #     `biomass_dilution` : dict 
            #         A target biomass at which the community should be diluted. 
            #         This must provide the target biomass `maximum` in approximate 
            #         OD units and the `minimum` OD unit which the sets the dilution
            #         factor.
        tol : float, optional
            The precision to tolerate for the integration. Values below this tolerance 
            will be cast to 0. Default value is 1E-10.
        solver_kwargs : dict 
            Keyword arguments to be passed to `scipy.integrate.solve_ivp`.

        Returns
        -------
        species_df : pandas DataFrame
            A dataframe with abundance information of the various species in 
            the culture.
        nutrient_df : pandas DataFrame
            A Dataframe with the concentrations of the nutrients in the culture.
        """

        if self._seed == False:
            raise RuntimeError("Ecosystem must first be seeded before growth. Call `seed()` method.")

        if extinction_thresh is not None:
            self.extinction_threshold = extinction_thresh
        else:
            self.extinction_threshold = None
        if len(bottleneck) != 0:
           species_df = pd.DataFrame([])
           nutrient_df = pd.DataFrame([])
           if bottleneck['type'] == 'time':    
               interval = bottleneck['interval']
               num_dil = int(np.floor(time / interval))

               # Partition the desired time range into intervals
               time_range = [[0, interval]]
               _total_time = interval
               for i in range(1, num_dil):
                   start = time_range[i-1][1]
                   if _total_time >= time:
                       end = time
                   else:
                       end = _total_time + interval
                   time_range.append([start, end])
                   _total_time += interval
               self._time_range = time_range 

               # Iterate through each dilution cycle and integrate. 
               print("Integrating dilution series:")
               for i, t in enumerate(self._time_range):
                   print(f"Integrating growth round {i+1} of {num_dil}...")

                   # Determine how to pack the initial state
                   if i == 0:
                       p0 = self._seed    
                   else:
                       _p0 = self.last_soln.y[:, -1]
                       p0 = [] 
                       _species, _, _total_mass = _unpack_masses(_p0, 
                                                                       self.num_species,
                                                                       self.num_nutrients)
                       if 'factor' in bottleneck.keys():
                           factor = bottleneck['factor']
                       elif 'target' in bottleneck.keys():
                           factor = (bottleneck['target'] * OD_CONV) / _total_mass
                       else:
                           raise RuntimeError("Must provide either a dilution factor or a target biomass minimum.")
 
                       for s in _species:
                           s[:-2] *= factor
                           for _s in s:
                               p0.append(_s)
                       for j in range(self.num_nutrients):
                           p0.append(self.init_concs[j])

                   # Integrate the time interval and store the resulting dataframes
                   if dt is not None:
                        solver_kwargs['t_eval'] = np.arange(t[0], t[1], dt)
                   _species_df,  _nutrient_df = self._integrate(t, p0,
                                                               solver_kwargs=solver_kwargs,
                                                               tshift=0,
                                                               tol=tol)

                   species_df = pd.concat([species_df, _species_df], sort=False)
                   nutrient_df = pd.concat([nutrient_df, _nutrient_df], sort=False)
               print('done!')
               return [species_df, nutrient_df]
        else:                     
            print("Integrating growth cycle...", end='')
            time_range = [0, time]
            if dt is not None:
                solver_kwargs['t_eval'] = np.arange(time_range[0], time_range[1], dt)
            species_df, nutrient_df = self._integrate(time_range, self._seed,
                                                  solver_kwargs=solver_kwargs,
                                                  tol=tol)
            print('done!', end='')
            if self.extinction:
                print('Exinction event has occured.')
        return species_df, nutrient_df

