if __name__ == '__main__':
    
    # |---------------------------|
    # |--- SET INPUT PARAMETERS --|
    # |---------------------------|
    
    # ZMAP catalog name and delimiter
    name_catalog = "ETAS_incomplete.txt"
    delimiter = ","

    # Moving window parameters
    size = 1000
    step = 250

    # Spatial map domain limits: lower-left longitude, lower-left latitude, upper-right longitude, upper-right latitude
    llcrnrlon = -123.5
    llcrnrlat = 30.
    urcrnrlon = -112.5
    urcrnrlat = 37.

    # b-value: comment this row if you want RESTORE to estimate the b-value for you! - DEFAULT: user's input
    b = 1

    # mc_ok: uncomment this row if you want to provide your reference value for the magnitude of completeness!
    # - DEFAULT: estimated with RESTORE.lilliefors(magcat_presequence, alpha)
    mc_ok = 2.0

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
    iidx_main = [j for j, v in enumerate(catalog_sel[:, MagnColumn]) if v == max(catalog_sel[:, MagnColumn])]
    starting_time_sequence = serial_times[iidx_main[0]]
    
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

