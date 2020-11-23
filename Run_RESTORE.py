if __name__ == '__main__':

    """
    This module loads the seismic catalog, sets the input parameters and runs RESTORE.py

    |--- CATALOG ---|

    The seismic catalog can be downloaded from web services based on FDSN specification. Alternatively, the user
    can directly load the catalog as a CSV file (in ZMAP format)

    |--- INPUT PARAMETERS ---|

    RESTORE.py requires the following parameters:

    1. Minimum magnitude in the seismic catalog

    2. Moving-window size (in number of events per window) --> 1000 by default (Mignan and Woessner, 2012)

    3. Moving-window step (in number of events per step) --> 250 by default (Mignan and Woessner, 2012)

    4. Space domain for the map plot
    
    5. Starting time of the seismic sequence (i.e. ending time of the seismically quiescent period)

    6. Reference value of the magnitude of completeness (estimated for the pre-sequence period):
       this could be either estimated with the Lilliefors test provided in RESTORE, or be set by the user

    7. b-value for the Gutenberg-Richter law: it could be either estimated with the function provided in RESTORE,
       or be set by the user

    """

    import RESTORE

    # |---------------|
    # |--- CATALOG ---|
    # |---------------|

    # -------------------- Download from a FDSN web service ------------------------------------------------- #

    # mmin =   # Minimum magnitude in the catalog
    # xmin =   # Minimum longitude
    # xmax =   # Maximum longitude
    # ymin =   # Minimum latitude
    # ymax =   # Maximum latitude
    # depthM =   # Maximum depth
    # time_start =   # String representing starting time in a recognizably valid format
    # time_end =   # String representing ending time in a recognizably valid format
    #
    # catalog_sel = RESTORE.acquisition_data(depthM, mmin, xmin, xmax, ymin, ymax, time_start, time_end)


    # -------------------- Load (ZMAP) catalog ------------------------------------------------- #

    name_catalog = " "
    delimiter = ' '
    catalog_sel = RESTORE.read_catalog_zmap(name_catalog, delimiter)

    # ------------------------------------------------------------------------------------------------------------- #

    # |------------------------|
    # |--- INPUT PARAMETERS ---|
    # |------------------------|

    mmin = 0.0  # Minimum magnitude in the catalog

    # Moving window parameters
    size = 1000
    step = 250

    # Spatial map domain limits:
    # lower-left longitude, lower-left latitude, upper-right longitude, upper-right latitude

    llcrnrlon = 00.00
    llcrnrlat = 00.00
    urcrnrlon = 00.00
    urcrnrlat = 00.00

    # ------------------------------------------------------------------------------------------------------------- #

    # |-----------------------------------------------------------------|
    # |--- REFERENCE VALUE FOR THE MAGNITUDE OF COMPLETENESS (mc_ok) ---|
    # |-----------------------------------------------------------------|

    MagnColumn = 5
    magcat = catalog_sel[:, MagnColumn]


    # ---------------- Use function lilliefors in RESTORE ------#

    serial_times = RESTORE.serial_time(catalog_sel)  # time stamps

    starting_time_sequence = RESTORE.serial_time('')  # starting time of the seismic sequence (time stamp)

    idx_pre_sequence = [idx for idx, val in enumerate(serial_times) if val < starting_time_sequence]

    magcat_presequence = magcat[idx_pre_sequence]

    alpha = 0.05  # significance level

    mc_ok, _ = RESTORE.lilliefors(magcat_presequence, mmin, alpha)  # reference value for mc

    # -------------------- User input ------------------------------------------------- #

    # mc_ok =   # reference value for mc

    # -------------------- m >= mc_ok ------------------------------------------------- #

    magcat_compl_idx = [idx for idx, val in enumerate(magcat) if val >= mc_ok]
    catalog_sel_compl = [catalog_sel[i, :] for i in magcat_compl_idx]

    # ------------------------------------------------------------------------------------------------------------- #

    # |---------------|
    # |--- B-VALUE ---|
    # |---------------|

    # -------------------- Use function 'b_value' in RESTORE ----- #

    # b, deltab, sigma, a, mag_cut = RESTORE.b_value(magcat, mc_ok)

    # -------------------- User input ----------------------------- #

    b = 1

    # ------------------------------------------------------------------------------------------------------------- #

    # |-----------------------|
    # |--- RUN RESTORE.py ---|
    # |-----------------------|

    replenished_catalog = RESTORE.replenished_catalog(catalog_sel, magcat, magcat_presequence, mmin, mc_ok, b, size, step,
                                                       llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, alpha, serial_times)

    mc_original_cat, _ = RESTORE.lilliefors(catalog_sel[:, MagnColumn], mc_ok, alpha)
    mc_replenished_cat, _ = RESTORE.lilliefors(replenished_catalog[:, MagnColumn], mc_ok, alpha)
