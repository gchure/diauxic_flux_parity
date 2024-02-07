import numpy as np
import pandas as pd


def flux_parity_system(params, args):

    return False

def unpack_params(params, args): 
    n_nuts = len(args['n_nutrients'])
    n_species = len(args(['n_species']))
    n_params = n_nuts + 4 # 4 is hard coded to account for M, M_o, TAA, TAA_star

    # Unpack the nutrient concentrations
    nutrients = [c for c in params[-n_nuts:]]
    masses = params[:-n_nuts]

    # Separate the params into species
    species = []
    for i in range(n_species):
        species.append(masses[i*n_params:n_params * (i + 1)])

    return nutrients, species
    

def self_replicator_FPM(params,
                        time,
                        args):
    """
    Defines the system of ordinary differenetial equations (ODEs) which describe 
    the self-replicator model with ppGpp regulation.

    Parameters
    ----------
    params: list, [M, Mr, Mp, (c_nt), T_AA, T_AA_star]
        A list of the parameters whose dynamics are described by the ODEs.
        M : positive float 
            Total biomass of the system
        M_Rb : positive float, must be < M 
            Ribosomal protein biomass of the system
        M_Mb : positive float, must be < M
            Metabolic protein biomass of the system 
        c_nt : positive float, optional
            The nutrient concentration in the environment. This should only 
            be provided if 'nutrients' is not False in the supplied arguments.
        T_AA_star : positive float
            Concentration of charged tRNAs in the culture. This is normalized to 
            total protein biomass.
        T_AA : positive float
            Concentration of uncharged tRNAs in the culture. This is normalized to 
            total protein biomass.
    time : float
        Evaluated time step of the system.
    args: dict 
        Dictionary of argument terms as follows
        gamma_max: positive float 
            The maximum translational capacity in units of inverse time.
        nu_max : positive float
            The maximum nutritional capacity in units of inverse time. 
        Kd_TAA : positive float
            The effective dissociation constant for uncharged tRNA to the metabolic 
            machinery. In units of abundance.
        Kd_TAA_star: positive float
            The effective dissociation constant for charged tRNA to the ribosome complex.
            In units of abundance
        kappa_max : positive float
            The maximum rate of uncharged tRNA synthesis in abundance units per unit time.
        phi_O : float, [0, 1], optional
            The fraction of the proteome occupied by 'other' protein mass.
        nutrients: bool or dict
            If False, nutrients will not be explicitly modeled and will be taken to 
            be saturating. If a dictionary is supplied, nutrients will be modeled 
            with following parameters
            Kd_cnc : float [0, inf)
                The effective dissociation constant of nutrients in the 
                to the metabolic machinery. 
            Y : float [0, inf)
                The yield coefficient of turning nutrients into precursors.

        dynamic_phiRb: bool or dict
            If True, phiRb will dynamically adjusted in reponse to charged/uncharged
            tRNA balance. If a dictionary is provided, seeded phiRb will be used.  
                phiRb: float [0, 1]
                    The seeded phiRb to be used.
        tRNA_regulation: bool
            if True, tRNA abundance will be regulated the same way as dynamic_phiRb.
            If False, kappa_max will be used directly. 
        antibiotic: bool
            If False, antiboitic presence will not be modeld and the fraction 
            of active ribosomes will be taken to be unity. If a dictionary is 
            provided with the following terms, the influence of antibiotics will
            be explicitly modeled.
                drug_conc : float [0, inf)
                    The concentration of the applied antibiotic
                Kd_drug : float [0, inf)
                    The effective dissociation constant of the drug to the ribosome.
        f_a : float, [0, 1]
            The faraction of ribosomes actively translating. 
        dil_approx: bool
            If True, then the approximation is made that the dilution of charged-tRNAs
            with growing biomass is negligible.

    Returns
    -------
    out: list, [dM_dt, dM_Rb_dt, dM_Mb_dt, (dc_nt_dt), dT_AA_dt, dT_AA_star_dt]
        A list of the evaluated ODEs at the specified time step.

        dM_dt : The dynamics of the total protein biomass.
        dM_Rb_dt : The dynamics of the ribosomal protein biomass.
        dM_Mb_dt : the dynamics of the metabolic protein biomass.
        dc_nt_dt : The dynamics of the nutrients in the environment, if modeled.
        dT_AA_dt : The dynamics of the uncharged tRNA concentration.
        dT_AA_star_dt : The dynamics of the uncharged tRNA concentration.
    """

    # Unpack the parameters
    if 'nutrients' in args.keys():
        M, M_Rb, M_Mb, c_nt, T_AA, T_AA_star = params
    else:
        M, M_Rb, M_Mb, T_AA, T_AA_star = params

    # Compute the capacities
    gamma = args['gamma_max'] * (T_AA_star / (T_AA_star + args['Kd_TAA_star']))
    if 'nutrients' in args.keys():
        pref = c_nt / (c_nt + args['nutrients']['Kd_cnt'])
    else:
        pref = 1
    nu = pref * args['nu_max'] * (T_AA / (T_AA + args['Kd_TAA']))

    # Compute the active fraction
    ratio = T_AA_star / T_AA

    fa = 1
    if 'antibiotic' in args.keys():
        fa -= args['antibiotic']['c_drug'] / \
            (args['antibiotic']['c_drug'] + args['antibiotic']['Kd_drug'])

    # Biomass accumulation
    if 'f_a' not in args.keys():
        f_a = 1
    else:
        f_a = args['f_a']
    dM_dt = f_a * gamma * M_Rb

    # Resource allocation
    if 'ansatz' in args.keys():
        if args['ansatz'] == 'binding':
            allocation = T_AA_star / (T_AA_star + T_AA)
    else:
        allocation = ratio / (ratio + args['tau'])

    if 'phiRb' not in args.keys():
        phiRb = (1 - args['phi_O']) * allocation
        kappa = args['kappa_max'] * allocation
    else:
        phiRb = args['phiRb']
        kappa = phiRb * args['kappa_max'] / (1 - args['phi_O'])

    dM_Rb_dt = phiRb * dM_dt
    dM_Mb_dt = (1 - phiRb - args['phi_O']) * dM_dt

    # Core tRNA dynamics
    dT_AA_star_dt = (nu * M_Mb - dM_dt) / M
    dT_AA_dt = (dM_dt - nu * M_Mb) / M

    # Dilution terms
    dT_AA_star_dt -= T_AA_star * dM_dt / M
    dT_AA_dt += kappa - (T_AA * dM_dt) / M

    if 'nutrients' in args.keys():
        dcnt_dt = -nu * M_Mb / args['nutrients']['Y']
        out = [dM_dt, dM_Rb_dt, dM_Mb_dt, dcnt_dt, dT_AA_dt, dT_AA_star_dt]
    else:
        out = [dM_dt, dM_Rb_dt, dM_Mb_dt, dT_AA_dt, dT_AA_star_dt]
    return out