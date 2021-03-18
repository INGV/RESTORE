if __name__ == '__main__':

    """
    This module loads the seismic catalog and the input parameters, then runs RESTORE.py

    |--- CATALOG ---|

    By default, the seismic catalog is loaded from the working directory. It must be a CSV file (in ZMAP format). 
    Alternatively, it can be downloaded from web services based on FDSN specification. 

    |--- INPUT PARAMETERS ---|

    RESTORE.py requires the following parameters:

    1. Moving-window size (in number of events per window) --> 1000 by default (Mignan and Woessner, 2012)
    2. Moving-window step (in number of events per step) --> 250 by default (Mignan and Woessner, 2012)
    3. Space domain for the map plot
    4. Starting time of the seismic sequence (i.e. ending time of the seismically quiescent period)
    5. Reference value for the magnitude of completeness, if set by the user (by default, it is estimated with the 
       function provided in RESTORE)
    6. b-value for the Gutenberg-Richter law, which by default must be provided by the user (alternatively, it could be 
       estimated with the function provided in RESTORE)
    7. alpha, the significance level for the Lilliefors test

    """
    
    # |---------------------------|
    # |--- SET INPUT PARAMETERS --|
    # |---------------------------|
    
    # ZMAP catalog name and delimiter
    name_catalog = ""
    delimiter = ""

    # Moving window parameters
    size = 1000
    step = 250

    # Spatial map domain limits: lower-left longitude, lower-left latitude, upper-right longitude, upper-right latitude
    llcrnrlon =
    llcrnrlat =
    urcrnrlon =
    urcrnrlat =
    
    # Starting time of the seismic sequence (e.g. "2016,8,24,0,0,0.0")
    tseq = ""
    
    # b-value: comment this row if you want RESTORE to estimate the b-value for you! - DEFAULT: user's input
    b =

    # mc_ok: uncomment this row if you want to provide your reference value for the magnitude of completeness!
    # - DEFAULT: estimated with RESTORE.lilliefors(magcat_presequence, alpha)
    # mc_ok =

    # Significance level for the Lilliefors test
    alpha = 0.05


    # ------------------------------------------------------------------------------------------------------------- #
    # ------------------------------------------------- RESTORE --------------------------------------------------- #
    # ------------------------------------------------------------------------------------------------------------- #

    import RESTORE

    # |---------------------------|
    # |--- LOAD (ZMAP) CATALOG ---|
    # |---------------------------|
  
    catalog_sel = RESTORE.read_catalog_zmap(name_catalog, delimiter)

    # |-----------------------------------------------------------------|
    # |--- REFERENCE VALUE FOR THE MAGNITUDE OF COMPLETENESS (mc_ok) ---|
    # |-----------------------------------------------------------------|

    MagnColumn = 5
    magcat = catalog_sel[:, MagnColumn]
    
    serial_times = RESTORE.serial_time(catalog_sel)  # time stamps
    starting_time_sequence = RESTORE.serial_time(tseq)  # timestamp for the starting time of the seismic sequence

    idx_pre_sequence = [idx for idx, val in enumerate(serial_times) if val < starting_time_sequence]

    magcat_presequence = magcat[idx_pre_sequence]

    if 'mc_ok' not in globals():
        mc_ok = RESTORE.lilliefors(magcat_presequence, alpha)

    # |-------------------------------------------------|
    # |--- COMPLETE CATALOG (m >= mc_ok) AND B-VALUE ---|
    # |-------------------------------------------------|

    magcat_compl_idx = [idx for idx, val in enumerate(magcat) if val >= mc_ok]
    catalog_sel_compl = [catalog_sel[i, :] for i in magcat_compl_idx]

    if 'b' not in globals():
        b, _, _, _, _ = RESTORE.b_value(magcat, mc_ok)
    
    # |-----------------------|
    # |--- RUN RESTORE.py ----|
    # |-----------------------|

    replenished_catalog = RESTORE.replenished_catalog(catalog_sel, magcat, magcat_presequence, mc_ok, b, size, step,
                                                       llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, alpha, serial_times)
                                                       
    # |---------------------------|
    # |--- LILLIEFORS CHECK.py ---|
    # |---------------------------|
    
    mc_original_cat = RESTORE.lilliefors(catalog_sel[:, MagnColumn], alpha)
    mc_replenished_cat = RESTORE.lilliefors(replenished_catalog[:, MagnColumn], alpha)






