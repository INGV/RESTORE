if __name__ == '__main__':

    """
    This module loads the seismic catalog and the input parameters, then runs RESTORE.py

    |--- CATALOG ---|

    By default, the seismic catalog is loaded from the working directory. It must be a CSV file (in ZMAP format). 
    Alternatively, it can be downloaded from web services based on FDSN specification. 

    |--- INPUT PARAMETERS ---|

    RESTORE.py requires the following parameters:

    1. Name of the catalog file
    2. Delimiter of the catalog file
    3. size: moving-window size (in number of events per window) --> 1000 by default (Mignan and Woessner, 2012)
    4. step: moving-window step (in number of events per step) --> 250 by default (Mignan and Woessner, 2012)
    5. st_dev_multiplier: multiplies the Mc standard deviation (Mc Sigma), controls the confidence level for STAI gaps identification; 
       STAI gaps are windows where mc >= mc_ok + n*sigma, where n = st_dev_multiplier; increasing the value of st_dev_multiplier 
       results in a more conservative approach in detecting temporary deviations of Mc
    5. Sigma: smoothing distance of the Gaussian kernel, controls the spread of the smoothing;
       smaller values will result in sharper and more localized smoothing
    6. sbin: bin in the latitude and longitude direction (degrees), controls the grid resolution
    7. fault_length_multiplier: multiplies the fault length rupture N times, controls the areal extent of the 
       subcatalog where inferences about Mc trend with time are made; 
       smaller values will result in higher resolution of STAI gaps detection, as local seismicity is less diluted
    8. Time of big shock causing STAI
    9. b [optional]: b-value of the Gutenberg-Richter law (alternatively, it can be estimated with the function provided in RESTORE)
    10. alpha: significance level for the Lilliefors test
    11. mc [optional]: reference value for the magnitude of completeness, if set by the user (by default, 
        it is estimated with the function provided in RESTORE)
    12. depth_distribution: [scipy.optimize.curve_fit] distribution to fit to hypocenter depths 
        --> normal, poisson, lognormal, beta, bimodal
    13. p0: [scipy.optimize.curve_fit] initial guess for the parameters of the hypo depth distribution:
        'mu', 'sigma' (normal), 
        'mu' (poisson), 
        'mu', 'sigma' (lognormal), 
        'a', 'b' (beta), 
        'mu1', 'sigma1', 'A1', 'mu2', 'sigma2', 'A2' [bimodal]
    """
    import config
    import RESTORE

    # |---------------------------|
    # |--- GET INPUT PARAMETERS --|
    # |---------------------------|
    
    config_dict = config.load_config('input_file.txt')

    name_catalog = config_dict['name_catalog']
    delimiter = config_dict['delimiter']
    size = config_dict['size']
    step = config_dict['step']
    st_dev_multiplier = config_dict['st_dev_multiplier']
    Sigma = config_dict['Sigma']
    sbin = config_dict['sbin']
    fault_length_multiplier = config_dict['fault_length_multiplier']
    t_end_quiet = config_dict['t_end_quiet']
    b = config_dict['b']
    alpha = config_dict['alpha']
    mc = config_dict['mc']
    depth_distribution = config_dict['depth_distribution']
    p0 = config_dict['p0'] 

    print(f"Confidence level for STAI gaps = {st_dev_multiplier} * Mc Sigma")
    print(f"Sigma (smoothing kernel) = {Sigma}")
    print(f"sbin (spatial grid resolution) = {sbin}°")
    print(f"alpha Lilliefors test = {alpha}")
    print(f"Hypo depth distribution to fit = {depth_distribution}")
    print(f"Hypo depth distribution intial params = {p0}")
    
    # ------------------------------------------------------------------------------------------------------------- #
    # --------------------------------------------- RUN RESTORE --------------------------------------------------- #
    # ------------------------------------------------------------------------------------------------------------- #

  
    catalog_sel = RESTORE.read_catalog_zmap(name_catalog, delimiter)
    print(f"----- Loading catalog: {name_catalog} -----")

    MagnColumn = 5
    
    serial_times = RESTORE.serial_time(catalog_sel)  # time stamps
    t_start = RESTORE.serial_time(t_end_quiet)  # timestamp for the starting time of the seismic sequence
    

    replenished_catalog = RESTORE.replenished_catalog(catalog_sel, t_start, mc, b, size, step,
                                                       st_dev_multiplier, alpha, serial_times, Sigma, sbin, 
                                                       fault_length_multiplier, depth_distribution, p0)
                                                    
    print("----- DONE -----")






