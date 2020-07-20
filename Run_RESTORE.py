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
    
    5. Reference value for the magnitude of completeness: this could be either estimated with the Lilliefors test 
       provided in RESTORE, or be set by the user
       
    6. b-value for the Gutenberg-Richter law: it could be either estimated with the function provided in RESTORE, 
       or be set by the user 
    
    """

    import RESTORE

    # |---------------|
    # |--- CATALOG ---|
    # |---------------|

    # -------------------- Download from a FDSN web service ------------------------------------------------- #

    mmin = 1.0  # Minimum magnitude in the catalog
    xmin = 13.10  # Minimum longitude
    xmax = 13.40  # Maximum longitude
    ymin = 42.20  # Minimum latitude
    ymax = 43.00  # Maximum latitude
    depthM = 30   # Maximum depth
    time_start = '2016-01-01T00:00:00'  # String representing starting time in a recognizably valid format
    time_end = '2016-09-30T00:00:00'  # String representing ending time in a recognizably valid format

    catalog_sel = RESTORE.acquisition_data(depthM, mmin, xmin, xmax, ymin, ymax, time_start, time_end)

    # -------------------- Load (ZMAP) catalog ------------------------------------------------- #

    #name_catalog = " "

    #catalog_sel = RESTORE.read_catalog_zmap(name_catalog)

    # ------------------------------------------------------------------------------------------------------------- #

    # |------------------------|
    # |--- INPUT PARAMETERS ---|
    # |------------------------|

    mmin = 1.0  # Minimum magnitude in the catalog

    # Moving window parameters
    size = 1000
    step = 250

    # Spatial map domain limits:
    # lower-left longitude, lower-left latitude, upper-right longitude, upper-right latitude

    # llcrnrlon = -123.5
    # llcrnrlat = 30.
    # urcrnrlon = -112.5
    # urcrnrlat = 37.

    llcrnrlon = 11.
    llcrnrlat = 41.
    urcrnrlon = 14.5
    urcrnrlat = 44.

    # ------------------------------------------------------------------------------------------------------------- #

    # |-----------------------------------------------------------------|
    # |--- REFERENCE VALUE FOR THE MAGNITUDE OF COMPLETENESS (mc_ok) ---|
    # |-----------------------------------------------------------------|

    MagnColumn = 5

    magcat = catalog_sel[:, MagnColumn]

    # -------------------- Use function lilliefors in RESTORE ------------------------------------------------- #

    mc_ok, _ = RESTORE.lilliefors(magcat, mmin)     # reference value for mc

    # -------------------- User input ------------------------------------------------- #

    #mc_ok = mmin

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

    replenished_catalog_Italy = RESTORE.replenished_catalog(catalog_sel, magcat, mmin, mc_ok, b, size, step,
                                                       llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat)
