"""

RESTORE: REal catalogs STOchastic REplenishment

TO CITE:
"Stallone, A., & Falcone, G. (2021). Missing earthquake data reconstruction in the space-time-magnitude domain. Earth and Space Science, 8(8), e2020EA001481."
https://doi.org/10.1029/2020EA001481

For any comment, question or suggestion write to:
angela.stallone@ingv.it

This project has been founded by the Seismic Hazard Center
(Centro di Pericolosità Sismica, CPS, at the Istituto Nazionale di Geofisica e Vulcanologia, INGV)

|--------------------------|
 Latest revision: Aug, 2026
|--------------------------|

RESTORE (REal catalogs STOchastic REplenishment) implements a stochastic gap-filling method
to quickly simulate time, hypocenter coordinates and magnitude of earthquakes that have not
been detected due to short term aftershock incompleteness (STAI).
The algorithm assesses the temporal variability of the magnitude of completeness Mc
with a sliding window approach. Since the window has a fixed number of events k and its
shift dk is constant, estimates of Mc are elapsed by dk events.
A statistic-based approach is implemented to pinpoint those time intervals where a threshold value for the
magnitude of completeness is significantly exceeded ("STAI gaps").
The number and magnitude of missing events are estimated by calculating the difference between the
observed counts and the counts predicted by the Gutenberg-Richter relationship.
Occurrence times and hypocenter depths are reconstructed implementing Monte
Carlo techniques to sample from empirical distribution functions.

"""

import csv
import numpy as np
import json
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pandas.plotting import register_matplotlib_converters
from mpl_toolkits.basemap import Basemap
from scipy.stats import norm, lognorm, beta
from scipy.optimize import curve_fit
import os
import math
import mc_lilliefors
from matplotlib.patches import Rectangle

__version__ = "3.0.0"

# Column numbers in the catalog
LonColumn = 0
LatColumn = 1
YearColumn = 2
MonthColumn = 3
DayColumn = 4
MagnColumn = 5
DepthColumn = 6
HourColumn = 7
MinuteColumn = 8
SecondColumn = 9


class NoSTAIGapsFound(RuntimeError):
    pass


def _force_rng(rng=None):
    if rng is None:
        return np.random.default_rng()
    return rng


def _write_run_summary(summary, filename="RESTORE_run_summary.json"):
    with open(filename, "w", encoding="utf-8") as output:
        json.dump(summary, output, indent=2)


def read_catalog_zmap(name_catalog, delimiter):
    """
    Read a 10-column ZMAP catalog with delimiter specified by the user.

    NOTE: A magnitude grid of 0.1 is required.
    """

    lon = []
    lat = []
    year = []
    month = []
    day = []
    mag = []
    depth = []
    hour = []
    minute = []
    sec = []
    with open(name_catalog) as csvfile:
        read_csv = csv.reader(csvfile, skipinitialspace=False, delimiter=delimiter)
        for row in read_csv:
            lon.append(float(row[0]))
            lat.append(float(row[1]))
            year.append(int(float(row[2])))
            month.append(int(float(row[3])))
            day.append(int(float(row[4])))
            mag.append(float(row[5]))
            depth.append(float(row[6]))
            hour.append(int(float(row[7])))
            minute.append(int(float(row[8])))
            sec.append(round(float(row[9]), 2))
    nrows2 = len(lon)
    ncols = 10
    catalog = np.ndarray((nrows2, ncols), dtype=object)
    catalog[:, LonColumn] = lon[:]
    catalog[:, LatColumn] = lat[:]
    catalog[:, YearColumn] = year[:]
    catalog[:, MonthColumn] = month[:]
    catalog[:, DayColumn] = day[:]
    catalog[:, MagnColumn] = mag[:]
    catalog[:, DepthColumn] = depth[:]
    catalog[:, HourColumn] = hour[:]
    catalog[:, MinuteColumn] = minute[:]
    catalog[:, SecondColumn] = sec[:]
    return catalog


# A magnitude binning of 0.1 is required. 
# Changing the magnitude binning of the original catalog to 0.1 is up to the user.
MAG_BIN = 0.1
MAG_TOL = 1e-8
def validate_magnitude_bins(catalog):
    """
    Validate magnitude binning in the catalog.
    """

    magnitudes = np.asarray(catalog[:, MagnColumn], dtype=float)
    if magnitudes.size == 0:
        raise ValueError("The catalog contains no events.")
    if not np.all(np.isfinite(magnitudes)):
        invalid = magnitudes[~np.isfinite(magnitudes)][:10]
        raise ValueError(
            "Catalog magnitudes must be finite. "
            f"Invalid examples: {invalid.tolist()}"
        )

    scaled = magnitudes / MAG_BIN
    invalid_grid = ~np.isclose(scaled, np.rint(scaled), atol=MAG_TOL, rtol=0)
    if np.any(invalid_grid):
        examples = magnitudes[invalid_grid][:10]
        raise ValueError(
            "RESTORE requires magnitudes already discretized at ΔM = 0.1. It does not re-bin magnitudes."
            f"Invalid examples: {examples.tolist()}"
        )


def validate_magnitude_value(value, name):
    """
    Validate magnitude binning.
    """

    if value is None:
        return
    value = float(value)
    if not np.isclose(value / MAG_BIN, np.rint(value / MAG_BIN), atol=MAG_TOL, rtol=0):
        raise ValueError(
            f"{name} must align with the default magnitude binning (ΔM = 0.1); received {value}."
        )
    

def plot_stai_region(catalog, lon_min, lon_max, lat_min, lat_max, mainshock_lon, mainshock_lat):
    """
    Plot the area covered by the input catalog along with the STAI rectangular area,
    i.e., the smaller region used for Mc(t), STAI-gap detection, and reconstruction.
    """

    lon = np.asarray(catalog[:, LonColumn], dtype=float)
    lat = np.asarray(catalog[:, LatColumn], dtype=float)

    llcrnrlon = lon.min()
    urcrnrlon = lon.max()
    llcrnrlat = lat.min()
    urcrnrlat = lat.max()

    fig, ax = plt.subplots(figsize=(8, 7))
    m = Basemap(projection="merc", resolution="i", llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, ax=ax)
    try:
        m.drawcoastlines(linewidth=0.5)
    except Exception:
        pass
    m.drawcountries()
    m.drawmapboundary()
    lon_range = urcrnrlon - llcrnrlon
    lat_range = urcrnrlat - llcrnrlat
    lon_interval = max(round(0.2 * lon_range, 2), 0.01)
    lat_interval = max(round(0.2 * lat_range, 2), 0.01)

    m.drawparallels(np.arange(-90, 90, lat_interval), labels=[1, 0, 0, 0], fontsize=10, dashes=[1, 5] )
    m.drawmeridians(np.arange(-180, 180, lon_interval), labels=[0, 0, 0, 1], fontsize=10, dashes=[1, 5])
    m.fillcontinents(color="#F0F0F0", lake_color="#F0F0F0")

    # Input catalog
    x, y = m(lon, lat)
    ax.plot(x, y, "o", color="0.4", markersize=3, alpha=0.5, label="Input catalog")

    # STAI rectangle 
    x1, y1 = m(lon_min, lat_min)
    x2, y2 = m(lon_max, lat_max)

    rectangle = Rectangle((x1, y1), x2 - x1, y2 - y1, fill=False, edgecolor="tab:blue", linewidth=2, label="STAI region")
    ax.add_patch(rectangle)

    # Largest earthquake
    x_main, y_main = m(mainshock_lon, mainshock_lat)
    ax.plot(x_main, y_main, marker="*", color="red", markersize=14, markeredgecolor="k", markeredgewidth=0.7, linestyle="None", label="Largest earthquake")
    ax.legend(loc="best")
    plt.savefig("fig/STAI_Region.pdf", format="pdf")
    plt.close()


def filter_catalog(catalog, t_start, serial_times, fault_length_multiplier):
    """
    Extract the spatial region used for STAI analysis.

    The largest event defines the center. Its estimated rupture length, multiplied by ``fault_length_multiplier``, defines the rectangular bounds.
    The initial catalog is filtered to extract the subset of events occurred in the STAI region.
    The returned pre-sequence magnitudes are inside the same region and occur before ``t_start``.
    """

    idx_big_shock = np.argmax(catalog[:, MagnColumn])
    mag_big_shock = catalog[idx_big_shock, MagnColumn]
    print(f"Largest event magnitude: {mag_big_shock}")
    lat_big_shock = catalog[idx_big_shock, LatColumn]
    lon_big_shock = catalog[idx_big_shock, LonColumn]

    # From Mw to seismic moment
    moment_big_shock = 10 ** (3 / 2 * (mag_big_shock + 10.7)) * 10 ** (-7)
    # Rupture length estiamtion (from Mai and Beroza, 2000)
    l_big_shock = 10 ** (-5.20 + 0.35 * math.log10(moment_big_shock))  

    # Max dist in lon and lat direction is fault_length_multiplier * rupture length 
    ymax, xmax = _trainv(lat_big_shock, lon_big_shock, fault_length_multiplier * l_big_shock, fault_length_multiplier * l_big_shock)
    ymin, xmin = _trainv(lat_big_shock, lon_big_shock, -fault_length_multiplier * l_big_shock, -fault_length_multiplier * l_big_shock)

    lat_min = round(ymin, 2)
    lat_max = round(ymax, 2)
    lon_min = round(xmin, 2)
    lon_max = round(xmax, 2)
    print(
        "STAI region bounds: "
        f"lat {lat_min} to {lat_max}, lon {lon_min} to {lon_max}"
    )

    mask = (ymin <= catalog[:, LatColumn]) & (catalog[:, LatColumn] <= ymax) & \
           (xmin <= catalog[:, LonColumn]) & (catalog[:, LonColumn] <= xmax)

    filtered_catalog = catalog[mask, :]

    # Extract magnitudes pre-sequence
    mask_preseq = (np.array(serial_times) < t_start) & (ymin <= catalog[:, LatColumn]) & (catalog[:, LatColumn] <= ymax) & \
           (xmin <= catalog[:, LonColumn]) & (catalog[:, LonColumn] <= xmax)
    
    magcat_presequence = catalog[mask_preseq, MagnColumn]

    # Plot extracted STAI region
    plot_stai_region(catalog, xmin, xmax, ymin, ymax, lon_big_shock, lat_big_shock)

    return filtered_catalog, mask, magcat_presequence


def write_filtered_catalog(filtered_catalog):
    """ 
    Write the STAI-filtered catalog to file.
    """

    with open('STAI-Filtered_Catalog.txt', 'w') as txtfile:
        for i, _ in enumerate(filtered_catalog):
            lon = filtered_catalog[i, LonColumn]
            lat = filtered_catalog[i, LatColumn]
            year = filtered_catalog[i, YearColumn]
            month = filtered_catalog[i, MonthColumn]
            day = filtered_catalog[i, DayColumn]
            mag = filtered_catalog[i, MagnColumn]
            depth = filtered_catalog[i, DepthColumn]
            hour = filtered_catalog[i, HourColumn]
            minute = filtered_catalog[i, MinuteColumn]
            second = filtered_catalog[i, SecondColumn]
            txtfile.write(f"{lon},{lat},{year},{month},{day},{mag},{depth},{hour},{minute},{second}\n")


def serial_time(input_):
    """
    Convert the input date (string format) or date vectors (year, month, day, hour, minute, second) into timestamps.

    Origin time: Year 0, Month 0, Day 0 (same as Matlab)
    """

    secs_per_day = 24.0 * 60.0 * 60.0

    # input: catalog
    if isinstance(input_, np.ndarray):
        catalog = input_
        out1 = []
        for i in range(len(catalog[:, YearColumn])):
            microseconds = round((catalog[i, SecondColumn] - int(catalog[i, SecondColumn])) * 1000000.0)
            dates = datetime.datetime(int(catalog[i, YearColumn]), int(catalog[i, MonthColumn]),
                                      int(catalog[i, DayColumn]), int(catalog[i, HourColumn]),
                                      int(catalog[i, MinuteColumn]), int(catalog[i, SecondColumn]),
                                      microseconds)
            mdn = dates + datetime.timedelta(days=366)
            frac_seconds = (dates - datetime.datetime(dates.year, dates.month, dates.day, 0, 0, 0)).seconds / (
                secs_per_day)
            frac_microseconds = dates.microsecond / (secs_per_day * 1000000.0)
            out1.append(mdn.toordinal() + frac_seconds + frac_microseconds)
        return out1

    # input: date (as string)
    elif isinstance(input_, str):
        input_date = input_
        time = datetime.datetime.strptime(input_date, '%Y,%m,%d,%H,%M,%S.%f')

        mdn = time + datetime.timedelta(days=366)
        frac_seconds = (time - datetime.datetime(time.year, time.month, time.day, 0, 0, 0)).seconds / (secs_per_day)
        frac_microseconds = time.microsecond / (secs_per_day * 1000000)
        out2 = mdn.toordinal() + frac_seconds + frac_microseconds
        return out2


def _dynamic_date_formatter(dates):
    max_date = max(dates)
    min_date = min(dates)
    date_range_days = (max_date - min_date).days

    if date_range_days > 365:  # Over a year
        fmt = '%b-%Y'
    elif date_range_days > 30:  # Over a month
        fmt = '%b-%d'
    else:  # Shorter intervals
        fmt = '%b-%d'
    return mdates.DateFormatter(fmt)


def _equally_spaced_ticks(dates, num_ticks=6):
    min_date = min(dates)
    max_date = max(dates)
    return mdates.drange(min_date, max_date, (max_date - min_date) / (num_ticks - 1))


def timestamp_to_datetime(timestamp):
    """
    Convert custom timestamps to datetime objects
    """

    days = timestamp % 1
    hours = days % 1 * 24
    minutes = hours % 1 * 60
    seconds = minutes % 1 * 60

    date_out = datetime.datetime.fromordinal(int(timestamp)) + datetime.timedelta(days=int(days)) + \
               datetime.timedelta(hours=int(hours)) + datetime.timedelta(minutes=int(minutes)) + \
               datetime.timedelta(seconds=seconds) - datetime.timedelta(days=366)

    return date_out


# Minimun number of events for Mc estimation with the Lilliefors test
MIN_EVENTS = 50
def lilliefors(magcat, alpha, delta_m=MAG_BIN, b=1.0, n_dithers=50, min_events=MIN_EVENTS, rng=None):
    """
    Estimate Mc with the Lilliefors test using truncated-exponential dithering.

    The truncated-exponential dithering used here follows:
    "Stallone, A., & Spassiani, I. (2026). Correcting exponentiality test for binned earthquake magnitudes. Seismica, 5(1)."
    https://doi.org/10.26443/seismica.v5i1.2257
    """

    magcat = np.asarray(magcat, dtype=float)
    mc, _ = mc_lilliefors.estimate_mc_lilliefors(
        magcat,
        delta_m=delta_m,
        b_value=b,
        alpha=alpha,
        n_dithers=n_dithers,
        min_events=min_events,
        rng=rng,
    )
    if mc is None:
        raise ValueError("The Lilliefors estimator could not determine Mc.")
    validate_magnitude_value(mc, "Estimated Mc")

    return mc


def bootstrap_mc(mag_data, alpha, b, rng=None):
    """
    Estimate Mc uncertainty with bootstrap.

    Return Mc standard deviation.
    """

    print("Estimating Mc uncertainty with bootstrap")
    rng = _force_rng(rng)
    iterations = 50
    mc_bootstrap = []
    for _ in range(iterations):
        boot = rng.choice(mag_data, size=len(mag_data), replace=True)
        boot_mc = lilliefors(boot, alpha, b=b, rng=rng)
        mc_bootstrap.append(boot_mc)
    mc_sigma_bootstrap = np.std(mc_bootstrap)  # mc standard deviation

    return mc_sigma_bootstrap


def mc_vs_time(magcat, magcat_presequence, mc_ok, st_dev_multiplier, alpha, serial_times, size, step, b, rng=None):
    """
    Estimate Mc through time and identify critical temporal windows ("STAI gaps"). 

    A STAI gap is identified when its Mc is at least ``mc_ok + st_dev_multiplier * sigma_mc``, where
    ``mc_ok`` is the reference pre-sequence Mc and ``sigma_mc`` is the st dev estimated by bootstrap. 
    Short isolated intervals are removed.
    """

    register_matplotlib_converters()

    rng = _force_rng(rng)
    mc_ok_sigma = bootstrap_mc(magcat_presequence, alpha, b=b, rng=rng)
    print(f"Mc standard deviation: {mc_ok_sigma}")
    upper_lim = mc_ok + st_dev_multiplier * mc_ok_sigma

    # Moving window to estimate Mc vs time
    # Mc is estimated every "step" events, at t_window (the end time of the moving window)
    print("Estimating Mc through time in the STAI-filtered catalog")

    window_starts = np.arange(0, len(magcat) - size + 1, step)
    
    mc_time = []
    t_window = []  # time of the moving window (represented by the time of the last event within each window)
    for j in window_starts:
        window = magcat[j: j + size]
        mc_new = lilliefors(window, alpha, b=b, rng=rng)
        mc_time.append(mc_new)
        t_window.append(serial_times[j + size -1])

    # Find the temporal bounds of the STAI gaps (where Mc >= mc_ok + n * sigma)
    tmp_hole_lower_lim = []
    tmp_hole_upper_lim = []
    for i in range(len(t_window) - 1):

        # If these conditions are met, Mc is exceeding the upper limit (Mc >= mc_ok + n * sigma)
        # --> STAI gap starts
        if i == 0:
            if mc_time[i] >= upper_lim:
                # Set the STAI gap starting time to the time of the first event if STAI is occurring 
                # already in the first window
                tmp_hole_lower_lim.append(serial_times[0])

        else:
            if (mc_time[i] >= upper_lim and mc_time[i - 1] < upper_lim):
                # Set the STAI gap starting time to the time of the largest shock in the step,
                # which has caused the raise of the magnitude of completeness
                idx_critical = [idx for idx, val in enumerate(serial_times) if t_window[i - 1] <= val <= t_window[i]]
                val = magcat[idx_critical]
                iidx_main = [j for j, v in enumerate(val) if v == max(val)]
                idx_main = idx_critical[iidx_main[0]]
                tmp_hole_lower_lim.append(serial_times[idx_main])

        # If these conditions are met STAI gap ends
        if i == len(t_window) - 2:
            if mc_time[i] >= upper_lim and mc_time[i + 1] >= upper_lim:
                tmp_hole_upper_lim.append(t_window[i + 1])
        else:
            if mc_time[i] >= upper_lim and mc_time[i + 1] < upper_lim:
                tmp_hole_upper_lim.append(t_window[i])

    # Exclude too small gaps, which could simply arise from Mc statistical fluctuations 
    hole_lower_lim , hole_upper_lim = [], []
    for j in range(len(tmp_hole_lower_lim)):
        index = [idx for idx, val in enumerate(serial_times) if
                 tmp_hole_lower_lim[j] <= val < tmp_hole_upper_lim[j]]
        if len(index) >= 2 * step:
            hole_lower_lim.append(tmp_hole_lower_lim[j])
            hole_upper_lim.append(tmp_hole_upper_lim[j])

    print(f"STAI gaps found: {len(hole_lower_lim)}")
    if len(hole_lower_lim) == 0:
        raise NoSTAIGapsFound(
            "No STAI gaps found. Consider the following options:\n"
            "1. Decrease the fault_length_multiplier\n"
            "2. Decrease the st_dev_multiplier\n"
            "3. Check if the reference Mc is too high: either decrease the duration of the seismically quiescent period or set Mc manually in the input file"
        )

    # Extract mc and relative times within each STAI gap
    mc_times_hole, mc_hole = [], []
    for j in range(len(hole_lower_lim)):
        index = [idx for idx, val in enumerate(t_window) if hole_lower_lim[j] <= val < hole_upper_lim[j]]
        tmp_mc_times_hole = [t_window[i] for i in index]
        tmp_mc_hole = [mc_time[i] for i in index]
        mc_times_hole.append(tmp_mc_times_hole)
        mc_hole.append(tmp_mc_hole)

    # plot Mc vs time
    t_window_plot = [val for idx, val in enumerate(t_window) if val >= hole_lower_lim[0]]
    idx_t_window_plot = [idx for idx, val in enumerate(t_window) if val >= hole_lower_lim[0]]
    dates1 = [timestamp_to_datetime(t_window_plot[i]) for i in range(len(t_window_plot))]
    datenums1 = mdates.date2num(dates1)
    mc_time_plot = [mc_time[i] for i in idx_t_window_plot]

    fig, ax = plt.subplots(figsize=(10, 6))
    ticks = _equally_spaced_ticks(dates1)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(_dynamic_date_formatter(dates1))
    ax.plot(datenums1, mc_time_plot, label='Mc(t)', ls='-')
    ax.fill_between(datenums1, upper_lim, mc_time_plot,
                    where=upper_lim <= mc_time_plot, color='C1', alpha=0.7,
                    label=rf"$M_c \geq M^*_c$ + {st_dev_multiplier} $\sigma$")
    ax.set_xlabel("Time", labelpad=10, fontsize=14) 
    ax.set_ylabel("Mc", labelpad=10, fontsize=14)  
    ax.legend(loc='upper right')
    plt.savefig('fig/Mc_Time.pdf', format='pdf')
    plt.close()
    
    return hole_lower_lim, hole_upper_lim, mc_times_hole, mc_hole


def map_plot(catalog, lon1, lat1, lon2, lat2):

    llcrnrlon, llcrnrlat = min(catalog[:, LonColumn]), min(catalog[:, LatColumn])
    urcrnrlon, urcrnrlat = max(catalog[:, LonColumn]), max(catalog[:, LatColumn])

    m = Basemap(projection='merc', resolution='i', llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
    try:
        m.drawcoastlines(linewidth=0.5)
    except Exception:
        pass
    m.drawcountries()
    m.drawmapboundary()

    lon_range = urcrnrlon - llcrnrlon
    lat_range = urcrnrlat - llcrnrlat
    lat_interval = max(round(0.2 * lat_range, 2), 0.01)
    lon_interval = max(round(0.2 * lon_range, 2), 0.01)
  
    parallels = np.arange(-90, 90, lat_interval)
    meridians = np.arange(-180, 180, lon_interval)
    m.drawparallels(parallels, labels=[1, 1, 1, 1], fontsize=10, dashes=[1, 5])
    m.drawmeridians(meridians, labels=[1, 1, 1, 1], fontsize=10, dashes=[1, 5])
    m.fillcontinents(color='#F0F0F0', lake_color='#F0F0F0')
    x1, y1 = m(lon1, lat1)
    x2, y2 = m(lon2, lat2)
    plt.plot(x1, y1, 'o', color='0.3', alpha=.7, markersize=5, markeredgecolor='w')
    plt.plot(x2, y2, 'o', color='steelblue', alpha=.7, markersize=5, markeredgecolor='w')
    plt.savefig('fig/Spatial_map.pdf', format='pdf')
    plt.close()

    return m

def spatial_smoothing(catalog, xmin, xmax, ymin, ymax, Sigma, sbin):
    """
    Convert observed epicenters into grid cell sampling probabilities.

    Each event contributes an isotropic Gaussian kernel with bandwidth ``Sigma`` (km). 
    The normalized grid rates are later used to sample locations of the reconstructed events.
    """

    if not np.isfinite(Sigma) or Sigma <= 0:
        raise ValueError("Sigma must be positive and finite for spatial smoothing.")
    if not np.isfinite(sbin) or sbin <= 0:
        raise ValueError("sbin must be positive and finite for spatial smoothing.")

    pi = np.pi
    sbinx = sbin
    sbiny = sbin

    rows_cat, columns_cat = np.shape(catalog)
    X = np.arange(xmin, xmax + sbinx, sbinx)
    Y = np.arange(ymin, ymax + sbiny, sbiny)
    xi, yi = np.meshgrid(X, Y)
    Smooth = np.zeros((len(X) * len(Y), 5))
    Smooth[:, 0] = np.reshape(xi, (len(X) * len(Y)))
    Smooth[:, 1] = np.reshape(yi, (len(X) * len(Y)))
    onetoten = range(0, rows_cat)
    for i in onetoten:
        A = (np.sin(Smooth[:, 1] * pi / 180) * np.sin(float(catalog[i, 1]) * pi / 180) + np.cos(
            Smooth[:, 1] * pi / 180) * np.cos(float(catalog[i, 1]) * pi / 180) * np.cos(
            ((float(catalog[i, 0])) - Smooth[:, 0]) * pi / 180))
        A = np.clip(A, -1.0, 1.0)
        R = np.arccos(A) * 6371
        W = (1 / (2 * pi * pow(Sigma, 2))) * (np.exp((-pow(R, 2)) / (2 * pow(Sigma, 2))))
        Smooth[:, 3] = Smooth[:, 3] + W

    TotalRate = sum(Smooth[:, 3])
    if not np.isfinite(TotalRate) or TotalRate <= 0:
        raise ValueError(
            "Spatial smoothing produced a zero or non-finite total rate. "
            "Check Sigma, sbin, and the filtered catalog spatial extent."
        )
    Smooth[:, 4] = Smooth[:, 3] / TotalRate

    lon_cell_grid = Smooth[:, 0]
    lan_cell_grid = Smooth[:, 1]
    smoothed_rate = Smooth[:, 4]

    return lon_cell_grid, lan_cell_grid, smoothed_rate


def _trainv(lat1, long1, distx, disty):

    r0 = 6367.
    alfa = 57.29578
    fi1 = lat1 + alfa * disty / r0
    fr1 = fi1 / alfa
    lat2 = fi1 - alfa * distx * distx / (2. * r0 * r0) * math.tan(lat1 / alfa)
    lon2 = long1 + alfa * distx / (r0 * math.cos(fr1)) - alfa * distx ** 3 / (3. * r0 ** 3) * (
        math.tan(fr1)) ** 2 / math.cos(fr1)

    return lat2, lon2


def magnitude_distribution(magcat):
    """
    Return bin and cumulative counts.
    """
        
    magcat = np.asarray(magcat, dtype=float)
    if magcat.size == 0:
        raise ValueError("Cannot calculate a magnitude distribution from an empty catalog")

    # Map every magnitude to an integer index on the 0.1 grid
    mag_bins = np.rint(magcat / MAG_BIN).astype(int)

    min_bin = mag_bins.min()
    max_bin = mag_bins.max()
    all_bins = np.arange(min_bin, max_bin + 1)

    # Count how many events fall in each integer magnitude index
    # (subtracting min_bin shifts the first magnitude index to zero, as required by np.bincount)
    bin_counts = np.bincount(mag_bins - min_bin, minlength=len(all_bins))

    nonzero = bin_counts > 0
    bin_counts = bin_counts[nonzero]
    # Convert the integer magnitude indices back to magnitude values
    magnitudes = all_bins[nonzero] * MAG_BIN

    log_bin_counts = np.log10(bin_counts)

    # For each magnitude level M, calculate the number of events with magnitude greater than or equal to M
    cum_counts = np.cumsum(bin_counts[::-1])[::-1]
    log_cum_counts = np.log10(cum_counts)

    return magnitudes, log_bin_counts, log_cum_counts


# Define possible hypo depth prob distributions
def fit_normal(x, mu, sigma):
    return norm.pdf(x, mu, sigma)

def fit_lognormal(x, mu, sigma):
    return lognorm.pdf(x, sigma, scale=np.exp(mu))

def fit_beta(x, a, b):
    return beta.pdf(x, a, b)

def fit_bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
    return A1 * np.exp(-(x - mu1)**2 / (2 * sigma1**2)) / (sigma1 * np.sqrt(2 * np.pi)) + \
           A2 * np.exp(-(x - mu2)**2 / (2 * sigma2**2)) / (sigma2 * np.sqrt(2 * np.pi))

def fit_distribution(distribution_function, x, y, p0):
    params, _ = curve_fit(distribution_function, x, y, p0=p0)
    return params

def validate_depth_fit_density(y_fit, depth_distribution):
    """ 
    Validate the fitted hypocentral depth density before random sampling. 

    The fitted density is later normalized and used as the probability distribution for sampling reconstructed hypocentral depths. 
    """

    y_fit = np.asarray(y_fit, dtype=float)
    if y_fit.size == 0:
        raise ValueError(f"{depth_distribution} depth fit produced an empty density vector.")
    if not np.all(np.isfinite(y_fit)):
        raise ValueError(f"{depth_distribution} depth fit produced non-finite density values.")
    if np.any(y_fit < 0):
        raise ValueError(f"{depth_distribution} depth fit produced negative density values.")
    if np.sum(y_fit) <= 0:
        raise ValueError(f"{depth_distribution} depth fit produced a zero total density.")
    return y_fit


def depth_fit_evaluation_grid(x_edges, depth_distribution):
    """
    Return x values where the fitted depth distribution is evaluated.

    For beta fits, avoid exact 0 and 1 because valid beta shapes can have
    infinite density at the endpoints.
    """

    if depth_distribution == "beta":
        return np.linspace(x_edges[1], x_edges[-2], 500)
    return np.linspace(x_edges.min(), x_edges.max(), 500)

    
def select_distribution(user_input):
    print("Fitting hypocentral depth distribution")
    if user_input == 'bimodal':
        print("Parameters to fit: 'mu1', 'sigma1', 'A1', 'mu2', 'sigma2', 'A2'")
        return fit_bimodal
    elif user_input == 'normal':
        print("Parameters to fit: 'mu', 'sigma'")
        return fit_normal
    elif user_input == 'lognormal':
        print("Parameters to fit: 'mu', 'sigma'")
        return fit_lognormal
    elif user_input == 'beta':
        print("Parameters to fit: 'a', 'b'")
        return fit_beta
    else:
        raise ValueError("Invalid distribution choice")


def ghost_element(catalog, t_start, mc_user, b_user, size, step, st_dev_multiplier,
                  alpha, serial_times, Sigma, sbin, fault_length_multiplier,
                  depth_distribution, p0, rng=None):
    """
    Detect STAI gaps and simulate missing events in the space-time-magnitude domain.

    Steps:
    - catalog spatial filtering to the STAI region
    - reference and time-dependent Mc estimation
    - Gutenberg-Richter count and magnitude reconstruction
    - stochastic sampling of times, epicenters, and hypocentral depths
    """

    rng = _force_rng(rng)

    filtered_catalog, mask, magcat_presequence = filter_catalog(catalog, t_start, serial_times, fault_length_multiplier)

    if len(magcat_presequence) < MIN_EVENTS:
        raise ValueError(
            f"The filtered pre-sequence catalog contains fewer than {MIN_EVENTS} events, which is the min number of events required for the Lilliefors test"
        )

    # Save filtered dataset to file 
    write_filtered_catalog(filtered_catalog)

    # If the user has not provided a b-value, use the standard default b = 1.0
    if b_user is None:
        b = 1.0
        print("b-value: using default b = 1.0")
    else:
        b = float(b_user)
        print(f"b-value: user-defined b = {b}")

    # If the user has not provided a reference Mc, estimate it from the pre-sequence magnitudes
    if mc_user is None:
        ref_mc = lilliefors(magcat_presequence, alpha, b=b, rng=rng)
        print(f"Reference Mc: estimated Mc = {ref_mc}")
    else:
        ref_mc = float(mc_user)
        print(f"Reference Mc: user-defined Mc = {ref_mc}")
        validate_magnitude_value(ref_mc, "User-defined reference Mc")

    magcat_filtered = filtered_catalog[:, MagnColumn]
    serial_times_filtered = np.array(serial_times)[mask]

    # "hole" --> STAI gap
    hole_lower_lim, hole_upper_lim, mc_times_hole, mc_hole = mc_vs_time(magcat_filtered, magcat_presequence, ref_mc,
                                                                        st_dev_multiplier, alpha, serial_times_filtered, size, step, b=b, rng=rng)

    
    # ----- Fit hypo depth distribution -----

    z_data = np.asarray(filtered_catalog[:, DepthColumn], dtype=float)

    if depth_distribution == "beta":
        # Beta distribution is defined on [0, 1]: normalize observed depths before fitting
        beta_depth_min = np.min(z_data)
        beta_depth_max = np.max(z_data)
        beta_depth_range = beta_depth_max - beta_depth_min
        if beta_depth_range <= 0:
            raise ValueError("Beta fitting requires hypocentral depths with a non-zero range.")
        z_fit_data = (z_data - beta_depth_min) / beta_depth_range
    else:
        z_fit_data = z_data

    bins = np.linspace(np.min(z_fit_data), np.max(z_fit_data), 101)

    if depth_distribution == "bimodal":
        y, x_edges = np.histogram(z_fit_data, bins=bins)
    else:
        y, x_edges = np.histogram(z_fit_data, bins=bins, density=True)

    # Histogram bin centers used for fitting
    x = (x_edges[1:] + x_edges[:-1]) / 2

    distribution_function = select_distribution(depth_distribution)
    params = fit_distribution(distribution_function, x, y, p0=p0)
    print(f"Fitted depth-distribution parameters: {params}")

    # Evaluate the fitted distribution over the same domain used for fitting
    x_fit = depth_fit_evaluation_grid(x_edges, depth_distribution)
    y_fit = validate_depth_fit_density(distribution_function(x_fit, *params), depth_distribution)

    print("Reconstructing missing events")

    # Convert the reference completeness magnitude to integer grid index on the 0.1 grid to avoid floating point truncation errors
    ref_mc_bin = int(np.rint(ref_mc / MAG_BIN)) 

    # Map every event magnitude to integer index on the 0.1 grid to avoid floating point truncation errors
    # Examples:
    # 1.0 -> 10
    # 1.1 -> 11
    # 1.4 -> 14
    mag_bins = np.rint(np.asarray(magcat_filtered, dtype=float) / MAG_BIN).astype(int)
    
    m_ghost, t_ghost, n_ghost_hole = [], [], []
    for i in range(len(hole_lower_lim)):  # iterate over all the STAI gaps

        # Number of ghost events estimated for the first (which is likely too short) step.
        # This value is temporarily stored and will later be merged with the second step.
        first_step_ghosts = 0

        n_ghost_step_list = []
        for j in range(len(mc_hole[i])):  # iterate over all the steps in the STAI gap

            mc_step = mc_hole[i][j]  # Mc value for a given step

            # index = event indices belonging to the current temporal step
            if j == 0:
                index = [idx for idx, val in enumerate(serial_times_filtered) if hole_lower_lim[i] < val <= mc_times_hole[i][j]]
            else:
                index = [idx for idx, val in enumerate(serial_times_filtered) if mc_times_hole[i][j - 1] < val <= mc_times_hole[i][j]]

            # ---- Calculate the total number of bins in the step ----

            # Convert the mc for the current step to integer index on the 0.1 grid to avoid floating point truncation errors
            mc_step_bin = int(np.rint(mc_step / MAG_BIN)) 
            # Calculate how many 0.1 magnitude bins fall between ref_mc and mc_step
            # nbin = mc_step_bin - ref_mc_bin 

            index_array = np.asarray(index, dtype=int)

            # Observed number of events complete above the current step Mc
            n_m2 = np.count_nonzero(mag_bins[index_array] >= mc_step_bin)

            if mc_step_bin <= ref_mc_bin:
                n_ghost_step = 0
            else:
                n_ghost_bin_list = []
                # Process all discrete magnitude levels below mc_step down to ref_mc
                # Example:
                # mc_step_bin = 15 and ref_mc_bin = 10
                # current_bin takes the values:
                # 14, 13, 12, 11, 10
                for current_bin in range(mc_step_bin - 1, ref_mc_bin - 1, -1):  # iterate over all the bins in the step

                    # Expected cumulative number at the next lower threshold
                    n_m1 = int(n_m2 * 10 ** (b * MAG_BIN))
                    # Expected number in the current magnitude bin 
                    n_exp_bin = n_m1 - n_m2

                    # Observed number in the current magnitude bin 
                    n_obs_bin = np.count_nonzero(mag_bins[index_array] == current_bin)

                    # Estimated missing events in the current magnitude bin
                    n_ghost_bin = n_exp_bin - n_obs_bin

                    if n_ghost_bin > 0:
                        ghost_magnitude = round(current_bin * MAG_BIN, 1)
                        m_ghost.extend([ghost_magnitude] * n_ghost_bin)
                        n_ghost_bin_list.append(n_ghost_bin)

                    # The predicted cumulative count at the lower threshold becomes the upper threshold count for the next iteration
                    n_m2 = n_m1

                n_ghost_step = sum(n_ghost_bin_list)  # number of ghost events in the current step

            # Times random sampling

            # Special case: the STAI gap consists of a single Mc step
            # --> sample all ghost events occurrence time uniformly over the whole gap
            if len(mc_hole[i]) == 1:
                n_to_sample = n_ghost_step
                lower_time = hole_lower_lim[i]
                upper_time = mc_times_hole[i][0]

            # The first step immediately after the mainshock may contain very few observed events
            # --> store its ghost event count and merge it with the second step 
            elif j == 0:
                first_step_ghosts = n_ghost_step
                n_ghost_step_list.append(n_ghost_step)
                continue

            # Merge first and second steps (covering [hole_lower_lim, mc_times_hole[i][1]])
            elif j == 1:
                n_to_sample = first_step_ghosts + n_ghost_step
                lower_time = hole_lower_lim[i]
                upper_time = mc_times_hole[i][1]

            # Case for all subsequent steps
            else:
                n_to_sample = n_ghost_step
                lower_time = mc_times_hole[i][j - 1]
                upper_time = mc_times_hole[i][j]

            # Occurrence times are sampled from a uniform distribution whose support are the limits of the step
            for _ in range(n_to_sample):
                u = rng.random()
                random_t = lower_time + u * (upper_time - lower_time)
                t_ghost.append(random_t)

            # Store the number of ghost events estimated for the current step
            # (the first step was already stored before the 'continue')
            n_ghost_step_list.append(n_ghost_step)

        print(
            f"STAI gap {i + 1}/{len(hole_lower_lim)}: "
            f"reconstructed {sum(n_ghost_step_list)} events"
        )

        n_ghost_hole.append(sum(n_ghost_step_list))  # append the number of ghost events in the current STAI gap

    # ---------------------    
    # ---- LON and LAT ----
    # --------------------- 
        
    xmin, ymin = min(filtered_catalog[:, LonColumn]), min(filtered_catalog[:, LatColumn])
    xmax, ymax = max(filtered_catalog[:, LonColumn]), max(filtered_catalog[:, LatColumn])
    lon_ghost, lat_ghost, z_ghost = [], [], []
    for i in range(len(hole_lower_lim)):  # iterate over all the STAI gaps
        lon_ghost_hole = []
        lat_ghost_hole = []

        index = [idx for idx, val in enumerate(serial_times_filtered) if hole_lower_lim[i] <= val <= hole_upper_lim[i]]
        subcatalog = filtered_catalog[index]
        lon_cell_grid, lat_cell_grid, smoothed_rate = spatial_smoothing(subcatalog, xmin, xmax, ymin, ymax, Sigma, sbin)

        # Simulate lat and lon of ghost eqs
        probabilities = np.asarray(smoothed_rate, dtype=float)
        if not np.all(np.isfinite(probabilities)) or probabilities.sum() <= 0:
            raise ValueError(
                "Spatial sampling probabilities are zero or non-finite."
            )
        probabilities = probabilities / probabilities.sum()

        # Number of eqs to simulate
        n_to_sample = n_ghost_hole[i]

        selected_cells = rng.choice(len(probabilities), size=n_to_sample, p=probabilities)
        for selected_cell in selected_cells:
            # Sample uniformly within the selected spatial grid cell
            rand_lon = (lon_cell_grid[selected_cell] + rng.random() * sbin)
            rand_lat = (lat_cell_grid[selected_cell] + rng.random() * sbin)

            lon_ghost_hole.append(round(rand_lon, 4))
            lat_ghost_hole.append(round(rand_lat, 4))

        lon_ghost.extend(lon_ghost_hole)
        lat_ghost.extend(lat_ghost_hole)

        # ---------------   
        # ---- DEPTH ----
        # ---------------

        samples = rng.choice(x_fit, size=n_to_sample, p=y_fit / np.sum(y_fit))
        if depth_distribution == "beta":
            # Convert beta samples from [0, 1] back to km
            samples = beta_depth_min + samples * beta_depth_range

        z_ghost.extend([round(value, 1) for value in samples])
        
    # For the beta plot, convert the fitted x-axis and density back to km
    if depth_distribution == "beta":
        x_fit_plot = (beta_depth_min + x_fit * beta_depth_range)
        y_fit_plot = y_fit / beta_depth_range
    else:
        x_fit_plot = x_fit
        y_fit_plot = y_fit
        
    # Plot hypo depth distribution

    # For beta, convert bins back to km for plotting
    if depth_distribution == "beta":
        bins_plot = beta_depth_min + bins * beta_depth_range
    else:
        bins_plot = bins

    data_min = np.min(z_ghost)
    data_max = np.max(z_ghost)
    # Zoom around the reconstructed hypocentral depth range
    margin = (data_max - data_min) * 0.95  
    depth_min, depth_max = data_min - margin, data_max + margin

    fig, ax = plt.subplots(figsize=(8, 6))
    if depth_distribution == 'bimodal':
        ax.hist(z_data, bins=bins_plot, alpha=0.3, color='0.3', label='Hypo depths (filtered) original')
        ax.plot(x_fit_plot, y_fit_plot, color='red', label=f"Fitted distribution: {depth_distribution}")
        ax.hist(z_ghost, bins=bins_plot, color='steelblue', label='Hypo depths replenished')
    else:
        ax.hist(z_data, bins=bins_plot, density=True, alpha=0.3, color='0.3', label='Hypo depths (filtered) original')
        ax.plot(x_fit_plot, y_fit_plot, color='red', label=f"Fitted distribution: {depth_distribution}")
        ax.hist(z_ghost, bins=bins_plot, density=True, color='steelblue', label='Hypo depths replenished')
    ax.set_xlim(depth_min, depth_max)
    ax.set_xlabel('Hypocenter Depth [km]', labelpad=10, fontsize=14)
    ax.legend()
    plt.savefig('fig/Hypo_Depth_Distribution.pdf', format='pdf')
    plt.close()
    
    ####    

    # plots
    # Plot original (m >= Mc) + simulated eqs for the meaningful region only 
    # i.e. where replenishment has been performed

    idx_magcatok_filtered = [idx for idx, val in enumerate(magcat_filtered) if val >= ref_mc]
    magcat_compl_filtered = [magcat_filtered[i] for i in idx_magcatok_filtered]
    serial_times_ok_filtered = [serial_times_filtered[i] for i in idx_magcatok_filtered]

    dates1 = [timestamp_to_datetime(serial_times_ok_filtered[i]) for i in range(len(serial_times_ok_filtered))]
    dates2 = [timestamp_to_datetime(t_ghost[i]) for i in range(len(t_ghost))]
    datenums1 = mdates.date2num(dates1)
    datenums2 = mdates.date2num(dates2)

    fig, ax = plt.subplots(2, figsize=(10, 6))
    fmt = mdates.DateFormatter('%b %Y')
    ax[0].plot(datenums1, magcat_compl_filtered, 'o', markersize=3, markerfacecolor='0.3', markeredgecolor='None',
               label='Original (filtered) data')
    ax[0].plot(datenums2, m_ghost, 'o', markersize=3, markerfacecolor='steelblue', markeredgecolor='None',
               label='Replenished data')
    ax[0].xaxis_date()
    ax[0].xaxis.set_major_formatter(fmt)
    ax[0].set_xlabel("Time")
    ax[0].set_ylabel("Magnitude")
    ax[0].legend(loc='upper left')
    ax[1].plot(datenums1, magcat_compl_filtered, 'o', markersize=3, markerfacecolor='0.3', markeredgecolor='None',
               label='Original (filtered) data')
    ax[1].xaxis_date()
    ax[1].xaxis.set_major_formatter(fmt)
    ax[1].set_xlabel("Time")
    ax[1].set_ylabel("Magnitude")
    ax[1].legend(loc='upper left')
    plt.subplots_adjust(hspace=0.5)
    plt.savefig('fig/Magnitude_Time.pdf', format='pdf')
    plt.close()

    # Magnitude-Time plots

    new_times = np.array(serial_times_ok_filtered + t_ghost)
    new_mag = np.array(magcat_compl_filtered + m_ghost)
    sorted_new_mag = new_mag[np.argsort(new_times)]

    x1 = np.arange(1, len(magcat_compl_filtered) + 1)
    x2 = np.arange(1, len(sorted_new_mag) + 1)

    fig, ax = plt.subplots()
    ax.set_xlabel("Sequential number")
    ax.set_ylabel("Magnitude")
    ax.set_title('Original (filtered) catalog')
    plt.plot(x1, magcat_compl_filtered, 'o', markersize=3, markerfacecolor='0.3', markeredgecolor='None')
    plt.savefig('fig/Magnitude_SeqNumbers_Original.pdf', format='pdf')
    plt.close()
    

    fig, ax = plt.subplots()
    ax.set_xlabel("Sequential number")
    ax.set_ylabel("Magnitude")
    ax.set_title('Replenished (filtered) catalog')
    plt.plot(x2, sorted_new_mag, 'o', markersize=3, markerfacecolor='0.3', markeredgecolor='None')
    plt.savefig('fig/Magnitude_SeqNumbers_Replenished.pdf', format='pdf')
    plt.close()
    
    # Map
    # Plot original (m >= Mc) + simulated eqs

    lat = filtered_catalog[:, LatColumn]
    lon = filtered_catalog[:, LonColumn]
    lat_compl = [lat[i] for i in idx_magcatok_filtered]
    lon_compl = [lon[i] for i in idx_magcatok_filtered]

    map_plot(filtered_catalog, lon_compl, lat_compl, lon_ghost, lat_ghost)
    
    stai_durations_days = [
        float(upper - lower)
        for lower, upper in zip(hole_lower_lim, hole_upper_lim)
    ]
    run_summary = {
        "reference_mc": float(ref_mc),
        "b_value": float(b),
        "b_value_source": "default" if b_user is None else "user",
        "n_stai_gaps": len(hole_lower_lim),
        "stai_durations_days": stai_durations_days,
        "total_stai_duration_days": float(np.sum(stai_durations_days)),
        "n_events_reconstructed": len(m_ghost),
    }

    return ref_mc, m_ghost, t_ghost, lon_ghost, lat_ghost, z_ghost, run_summary


def replenished_catalog(catalog, t_start, mc_user, b_user, size, step, st_dev_multiplier,
                        alpha, serial_times, Sigma, sbin, fault_length_multiplier,
                        depth_distribution, p0, rng=None):
    """
    Create and write one replenished catalog realization.
    
    Original events below the reference Mc are excluded; reconstructed events
    are merged with the remaining original events, sorted by occurrence time,
    and marked with ``flag=1``. Original events have ``flag=0``.
    """

    rng = _force_rng(rng)

    if size < MIN_EVENTS:
        raise ValueError(
            f"The moving-window size must be at least min_events = {MIN_EVENTS}")

    # Validate magnitude binning (must be 0.1)
    validate_magnitude_bins(catalog)

    path = "fig"
    if not os.path.exists(path):
        os.mkdir(path)

    ref_mc, m_ghost, t_ghost, lon_ghost, lat_ghost, z_ghost, run_summary = ghost_element(catalog, t_start, mc_user, b_user,
                                                           size, step, st_dev_multiplier, alpha, serial_times, Sigma, sbin, fault_length_multiplier,
                                                           depth_distribution, p0, rng=rng)

    magcat = catalog[:, MagnColumn]
    magcat_compl_idx = [idx for idx, val in enumerate(magcat) if val >= ref_mc]

    magcat_compl = [magcat[i] for i in magcat_compl_idx]
    serial_times_compl = [serial_times[i] for i in magcat_compl_idx]

    new_times = np.array(serial_times_compl + t_ghost)
    sorted_new_times = new_times[np.argsort(new_times)]

    new_mag = np.array(magcat_compl + m_ghost)
    sorted_new_mag = new_mag[np.argsort(new_times)]

    lat = catalog[:, LatColumn]
    lon = catalog[:, LonColumn]
    depth = catalog[:, DepthColumn]
    lat_compl = [lat[i] for i in magcat_compl_idx]
    lon_compl = [lon[i] for i in magcat_compl_idx]
    depth_compl = [depth[i] for i in magcat_compl_idx]

    new_lat = np.array(lat_compl + lat_ghost)
    sorted_new_lat = new_lat[np.argsort(new_times)]

    new_lon = np.array(lon_compl + lon_ghost)
    sorted_new_lon = new_lon[np.argsort(new_times)]

    new_depth = np.array(depth_compl + z_ghost)
    sorted_new_depth = new_depth[np.argsort(new_times)]

    # flag: 0 for original events, 1 for simulated events
    flag = np.array([0] * len(new_times))
    for i in range(len(new_times)):
        if i < len(serial_times_compl):
            flag[i] = 0
        else:
            flag[i] = 1

    sorted_flag = flag[np.argsort(new_times)]

    sorted_new_year, sorted_new_month, sorted_new_day, sorted_new_hour, sorted_new_min, sorted_new_sec = \
        [], [], [], [], [], []

    for i in range(len(new_times)):
        date_out = timestamp_to_datetime(sorted_new_times[i])
        sorted_new_year.append(int(date_out.year))
        sorted_new_month.append(int(date_out.month))
        sorted_new_day.append(int(date_out.day))
        sorted_new_hour.append(int(date_out.hour))
        sorted_new_min.append(int(date_out.minute))
        sec = float(date_out.second) + (float(date_out.microsecond) / 1000000)
        sorted_new_sec.append(int(sec*100)/100.)

    nrows = len(new_times)
    ncols = 11
    replenished_catalog = np.ndarray((nrows, ncols), dtype=object)
    replenished_catalog[:, LonColumn] = sorted_new_lon[:]
    replenished_catalog[:, LatColumn] = sorted_new_lat[:]
    replenished_catalog[:, YearColumn] = sorted_new_year[:]
    replenished_catalog[:, MonthColumn] = sorted_new_month[:]
    replenished_catalog[:, DayColumn] = sorted_new_day[:]
    replenished_catalog[:, MagnColumn] = sorted_new_mag[:]
    replenished_catalog[:, DepthColumn] = sorted_new_depth[:]
    replenished_catalog[:, HourColumn] = sorted_new_hour[:]
    replenished_catalog[:, MinuteColumn] = sorted_new_min[:]
    replenished_catalog[:, SecondColumn] = sorted_new_sec[:]
    replenished_catalog[:, ncols - 1] = sorted_flag[:]

    rows = zip(sorted_new_lon, sorted_new_lat, sorted_new_year, sorted_new_month, sorted_new_day,
               sorted_new_mag, sorted_new_depth, sorted_new_hour, sorted_new_min, sorted_new_sec, sorted_flag)
    with open('Replenished_catalog.txt', 'w') as f:
        writer = csv.writer(f, delimiter=' ', skipinitialspace=False, quoting=csv.QUOTE_NONE, lineterminator='\n')
        for x in rows:
            writer.writerow(x)

    magnitudes_original, log_bin_counts_original, log_cum_counts_original = \
        magnitude_distribution(catalog[:, MagnColumn])
    magnitudes_replenished, log_bin_counts_replenished, log_cum_counts_replenished = \
        magnitude_distribution(replenished_catalog[:, MagnColumn])

    fig, ax = plt.subplots(figsize=(8, 10))
    ax.plot(magnitudes_original, log_bin_counts_original, 'o', color='0.5', markersize=3, label='Non-cumulative counts (Original)', alpha=0.7)
    ax.plot(magnitudes_original, log_cum_counts_original, '-', color='0.5', label='Cumulative counts (Original)', alpha=0.7)
    ax.plot(magnitudes_replenished, log_bin_counts_replenished, 'o', color='black', markersize=3, label='Non-cumulative counts (Replenished)', alpha=0.7)
    ax.plot(magnitudes_replenished, log_cum_counts_replenished, '-', color='black', label='Cumulative counts (Replenished)', alpha=0.7)
    ax.set_xlabel("Magnitude", labelpad=10, fontsize=14)
    ax.set_ylabel("Frequency", labelpad=10, fontsize=14)
    ax.legend(loc='upper right')
    plt.savefig('fig/Magnitude_Frequency_Distribution.pdf', format='pdf')
    plt.close()

    _write_run_summary(run_summary)
    
    return replenished_catalog
