import RESTORE

# |---------------|
# |--- CATALOG ---|
# |---------------|

name_catalog = "ETAS_incomplete.txt"
delimiter = ','
catalog_sel = RESTORE.read_catalog_zmap(name_catalog, delimiter)

# |------------------------|
# |--- INPUT PARAMETERS ---|
# |------------------------|

mmin = 2.0  # Minimum magnitude in the catalog

# Moving window parameters
size = 1000
step = 250

# Spatial map domain limits:
# lower-left longitude, lower-left latitude, upper-right longitude, upper-right latitude

llcrnrlon = -123.5
llcrnrlat = 30.
urcrnrlon = -112.5
urcrnrlat = 37.

# |-----------------------------------------------------------------|
# |--- REFERENCE VALUE FOR THE MAGNITUDE OF COMPLETENESS (mc_ok) ---|
# |-----------------------------------------------------------------|

MagnColumn = 5
magcat = catalog_sel[:, MagnColumn]

serial_times = RESTORE.serial_time(catalog_sel)  # time stamps

iidx_main = [j for j, v in enumerate(catalog_sel[:, MagnColumn]) if v == max(catalog_sel[:, MagnColumn])]
starting_time_sequence = serial_times[iidx_main[0]]  # starting time of the seismic sequence (time stamp)

idx_pre_sequence = [idx for idx, val in enumerate(serial_times) if val < starting_time_sequence]
magcat_presequence = magcat[idx_pre_sequence]

alpha = 0.05  # significance level

mc_ok = 2.0  # reference value for mc

# -------------------- m >= mc_ok ------------------------------------------------- #

magcat_compl_idx = [idx for idx, val in enumerate(magcat) if val >= mc_ok]
catalog_sel_compl = [catalog_sel[i, :] for i in magcat_compl_idx]

# |---------------|
# |--- B-VALUE ---|
# |---------------|

b = 1

# |-----------------------|
# |--- RUN RESTORE.py ---|
# |-----------------------|

replenished_catalog = RESTORE.replenished_catalog(catalog_sel, magcat, magcat_presequence, mmin, mc_ok, b, size, step,
                                                  llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, alpha, serial_times)

mc_original_cat, _ = RESTORE.lilliefors(catalog_sel[:, MagnColumn], mc_ok, alpha)
mc_replenished_cat, _ = RESTORE.lilliefors(replenished_catalog[:, MagnColumn], mc_ok, alpha)