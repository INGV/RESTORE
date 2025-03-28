"""

RESTORE: REal catalogs STOchastic REplenishment
Written by Angela Stallone with help from Giuseppe Falcone
For any comment, question or suggestion write to:
angela.stallone@ingv.it

This project has been founded by the Seismic Hazard Center
(Centro di Pericolosità Sismica, CPS, at the Istituto Nazionale di Geofisica e Vulcanologia, INGV)

|---------------------------|

Latest revision: Jan, 2024

|--------------------------------------------------------------------------|

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
from urllib.request import urlopen
import json
import dateparser
import csv
from numpy import flipud
import xmltodict
import datetime
import pandas as pd
import sys
import matplotlib.pyplot as plt
from math import sqrt, log, log10
import matplotlib.dates as mdates
from pandas.plotting import register_matplotlib_converters
from mpl_toolkits.basemap import Basemap
from scipy.stats import norm, lognorm, beta
from scipy.optimize import curve_fit
import os
import math
import decimal
decimal.getcontext().prec = 2

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


def read_catalog_zmap(name_catalog, delimiter):
    
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
            mag.append(round(float(row[5]), 1))
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


def acquisition_data(depthm, mmin, xmin, xmax, ymin, ymax, time_start, time_end):

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
    type_mag = []
    url = "http://webservices.rm.ingv.it/fdsnws/event/1/query?starttime=" + str(time_start) + "&endtime=" + str(
        time_end) + "&minlat=" + str(ymin) + "&maxlat=" + str(ymax) + "&minlon=" + str(xmin) + "&maxlon=" + str(
        xmax) + "&maxdepth=" + str(depthm * 1000) + "&minmag=" + str(mmin) + "&format=geojson"
    print(url)
    html = urlopen(url).read()
    data = json.loads(html)
    field_list = data['features']
    i = 0
    for field in field_list:
        mytime = data['features'][i]['properties']['time']
        time = dateparser.parse(mytime)
        year.append(int(time.year))
        month.append(int(time.month))
        day.append(int(time.day))
        hour.append(int(time.hour))
        minute.append(int(time.minute))
        second = float(time.second) + (float(time.microsecond) / 1000000)
        sec.append(round(second, 2))
        lon.append(float(data['features'][i]['geometry']['coordinates'][0]))
        lat.append(float(data['features'][i]['geometry']['coordinates'][1]))
        depth.append(float(data['features'][i]['geometry']['coordinates'][2]))
        mag.append(data['features'][i]['properties']['mag'])
        type_mag.append(data['features'][i]['properties']['magType'])
        i = i + 1
    nrows2 = len(lon)
    ncols = 11
    catalog = np.ndarray((nrows2, ncols), dtype=object)
    catalog[:, LonColumn] = flipud(lon[:])
    catalog[:, LatColumn] = flipud(lat[:])
    catalog[:, YearColumn] = flipud(year[:])
    catalog[:, MonthColumn] = flipud(month[:])
    catalog[:, DayColumn] = flipud(day[:])
    catalog[:, MagnColumn] = flipud(mag[:])
    catalog[:, DepthColumn] = flipud(depth[:])
    catalog[:, HourColumn] = flipud(hour[:])
    catalog[:, MinuteColumn] = flipud(minute[:])
    catalog[:, SecondColumn] = flipud(sec[:])
    rows = zip(flipud(lon), flipud(lat), flipud(year), flipud(month), flipud(day), flipud(mag), flipud(depth),
               flipud(hour), flipud(minute), flipud(sec))
    with open('Zmap_catalog', 'w', newline='') as f:
        writer = csv.writer(f, delimiter=' ', skipinitialspace=True, quoting=csv.QUOTE_NONE)
        for x in rows:
            writer.writerow(x)
    return catalog


def acquisition_xml(depth_m, mmin, xmin, xmax, ymin, ymax, time_start, time_end):

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
    type_mag = []
    url = "http://webservices.rm.ingv.it/fdsnws/event/1/query?starttime=" + str(time_start) + "&endtime=" + str(
        time_end) + "&minlat=" + str(ymin) + "&maxlat=" + str(ymax) + "&minlon=" + str(xmin) + "&maxlon=" + str(
        xmax) + "&maxdepth=" + str(depth_m * 1000) + "&minmag=" + str(mmin) + "&includeallmagnitudes"
    print(url)
    html = urlopen(url).read()
    doc = xmltodict.parse(html)
    root = doc["q:quakeml"]
    eventparameters = root["eventParameters"]
    i = 0
    for data in eventparameters["event"]:
        mytime = data['origin']['time']['value']
        time = dateparser.parse(mytime)
        year.append(int(time.year))
        month.append(int(time.month))
        day.append(int(time.day))
        hour.append(int(time.hour))
        minute.append(int(time.minute))
        second = float(time.second) + (float(time.microsecond) / 1000000)
        sec.append(round(second, 2))
        lon.append(float(data['origin']['longitude']['value']))
        lat.append(float(data['origin']['latitude']['value']))
        depth.append(float(data['origin']['depth']['value']) / 1000)
        mag.append(data['magnitude']['mag']['value'])
        type_mag.append(data['magnitude']['type'])
        i = i + 1
    nrows2 = len(lon)
    ncols = 11
    catalog = np.ndarray((nrows2, ncols), dtype=object)
    catalog[:, LonColumn] = flipud(lon[:])
    catalog[:, LatColumn] = flipud(lat[:])
    catalog[:, YearColumn] = flipud(year[:])
    catalog[:, MonthColumn] = flipud(month[:])
    catalog[:, DayColumn] = flipud(day[:])
    catalog[:, MagnColumn] = flipud(mag[:])
    catalog[:, DepthColumn] = flipud(depth[:])
    catalog[:, HourColumn] = flipud(hour[:])
    catalog[:, MinuteColumn] = flipud(minute[:])
    catalog[:, SecondColumn] = flipud(sec[:])
    rows = zip(flipud(lon), flipud(lat), flipud(year), flipud(month), flipud(day), flipud(mag), flipud(depth),
               flipud(hour), flipud(minute), flipud(sec))
    with open('Zmap_catalog', 'w', newline='') as f:
        writer = csv.writer(f, delimiter=' ', skipinitialspace=True, quoting=csv.QUOTE_NONE)
        for x in rows:
            writer.writerow(x)
    return catalog


def filter_catalog(catalog, t_start, serial_times, fault_length_multiplier):
    # Filter catalog to a meaningful region surrounding the large shock

    idx_big_shock = np.argmax(catalog[:, MagnColumn])
    mag_big_shock = catalog[idx_big_shock, MagnColumn]
    print("Mag_big_shock = ", mag_big_shock)
    lat_big_shock = catalog[idx_big_shock, LatColumn]
    lon_big_shock = catalog[idx_big_shock, LonColumn]

    moment_big_shock = 10 ** (3 / 2 * (mag_big_shock + 10.7)) * 10 ** (-7)
    l_big_shock = 10 ** (-5.20 + 0.35 * math.log10(moment_big_shock))  # rupture length (from Mai and Beroza (2000))

    # Max dist in lon and lat direction is fault_length_multiplier * rupture length 
    ymax, xmax = trainv(lat_big_shock, lon_big_shock, fault_length_multiplier * l_big_shock, fault_length_multiplier * l_big_shock)
    ymin, xmin = trainv(lat_big_shock, lon_big_shock, -fault_length_multiplier * l_big_shock, -fault_length_multiplier * l_big_shock)
    print(f"Subcatalog limits: "
      f"LAT_MIN = {round(ymin, 2)}, LAT_MAX = {round(ymax, 2)} "
      f"LON_MIN = {round(xmin, 2)}, LON_MAX = {round(xmax, 2)}")

    mask = (ymin <= catalog[:, LatColumn]) & (catalog[:, LatColumn] <= ymax) & \
           (xmin <= catalog[:, LonColumn]) & (catalog[:, LonColumn] <= xmax)

    filtered_catalog = catalog[mask, :]

    # Extract magnitudes pre-sequence
    mask_preseq = (np.array(serial_times) < t_start) & (ymin <= catalog[:, LatColumn]) & (catalog[:, LatColumn] <= ymax) & \
           (xmin <= catalog[:, LonColumn]) & (catalog[:, LonColumn] <= xmax)
    
    magcat_presequence = catalog[mask_preseq, MagnColumn]

    return filtered_catalog, mask, magcat_presequence


def write_filterd_catalog(filtered_catalog):
    # Write filtered catalog to file

    with open('Filtered_Catalog.txt', 'w') as txtfile:
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
    # Converts the input date (string format) or date vectors (year, month, day, hour, minute, second) into timestamps
    # Origin time: Year 0, Month 0, Day 0 (same as Matlab)

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


def dynamic_date_formatter(dates):
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


def equally_spaced_ticks(dates, num_ticks=6):
    min_date = min(dates)
    max_date = max(dates)
    return mdates.drange(min_date, max_date, (max_date - min_date) / (num_ticks - 1))


def timestamp_to_datetime(timestamp):
    # Converts custom timestamps to datetime objects

    days = timestamp % 1
    hours = days % 1 * 24
    minutes = hours % 1 * 60
    seconds = minutes % 1 * 60

    date_out = datetime.datetime.fromordinal(int(timestamp)) + datetime.timedelta(days=int(days)) + \
               datetime.timedelta(hours=int(hours)) + datetime.timedelta(minutes=int(minutes)) + \
               datetime.timedelta(seconds=seconds) - datetime.timedelta(days=366)

    return date_out


def lilliefors(magcat, alpha):
    # Implements the Python routine "mc_lilliefors" 
    # Source: Herrmann, M. and W. Marzocchi (2021). Inconsistencies and Lurking Pitfalls 
    # in the Magnitude–Frequency Distribution of High-Resolution Earthquake Catalogs. 
    # Seismological Research Letters 92(2A). doi: 10.1785/0220200337
    # Gitlab: https://gitlab.com/marcus.herrmann/mc-lilliefors

    try:
        import mc_lilliefors
    except ImportError:
        print("Ops! Module mc_lilliefors is missing! Download it here: https://gitlab.com/marcus.herrmann/mc-lilliefors")
        sys.exit(1)

    mag_series = pd.Series(magcat, dtype='float64')

    lill = mc_lilliefors.McLilliefors(
        mag_series,
        alpha
    )

    lill.calc_testdistr_mcutoff(
        n_repeats=50  # number of iterations for the random noise
       )

    mc = lill.estimate_Mc_expon_test()

    return mc


def b_value_zmap(magcat):
    # Calculate the b_value, tha a_value and the completeness magnitude using the ZMAP approach

    bin_m = 0.1
    magmin = 0.0
    magmax = max(magcat)
    x = np.arange(magmin, magmax + bin_m, bin_m)
    nevents = len(magcat)
    y = np.arange(0, nevents)
    rad_mag_cutoff = np.zeros(shape=(len(x), nevents))
    ntotmag = np.zeros(shape=(len(x)))
    difference = np.zeros(shape=(len(x)))
    sumup = np.zeros(shape=(len(x)))
    rad_mag_cutoff.fill(np.nan)
    idx = -1
    for i in x:
        idx = idx + 1
        for j in y:
            if magcat[j] >= i:
                ntotmag[idx] = ntotmag[idx] + 1
                rad_mag_cutoff[idx][j] = magcat[j]
    mean_mag_cutoff = np.nanmean(rad_mag_cutoff, axis=1)
    idx = -1
    for i in x:
        idx = idx + 1
        for j in y:
            rad_mag_cutoff[idx][j] = (rad_mag_cutoff[idx][j] - mean_mag_cutoff[idx]) ** 2
    sum_mag_cutoff = np.nansum(rad_mag_cutoff, axis=1)
    bin_centers = 0.5 * (x[1:] + x[:-1])
    no_cumulative, av, b = plt.hist(magcat, bin_centers)
    avm = np.delete(av, [len(av) - 1])
    # B_value from Marzocchi and Sandri (2003)
    b_cutoff = 1 / (log(10) * ((mean_mag_cutoff) - (x - bin_m / 2)))
    deltab_value = 2.3 * (b_cutoff ** 2) * np.sqrt(sum_mag_cutoff / (ntotmag * (ntotmag - 1)))
    d_m1 = 0.5
    dm = 0.1
    b_ave = [(a + b + c + d + e) * dm / d_m1 for a, b, c, d, e in
             zip(b_cutoff[0::1], b_cutoff[1::1], b_cutoff[2::1], b_cutoff[3::1], b_cutoff[4::1])]
    dbi_old = [(a - b) for a, b in zip(b_cutoff[0::1], b_cutoff[1::1])]
    idx = -1
    q = np.arange(0, len(b_ave))
    flag = 0
    for i in q:
        idx = idx + 1
        difference[idx] = abs(b_ave[idx] - b_cutoff[idx])
        sumup[idx] = b_ave[idx] + b_cutoff[idx]
        if difference[idx] <= deltab_value[idx] and abs(dbi_old[idx]) <= 0.03:
            sigmacomp = b_cutoff[idx] / sqrt(ntotmag[idx])
            avaluecomp = log10(ntotmag[idx]) + b_cutoff[idx] * x[idx]
            if flag == 0:
                print('The completeness magnitude is=', x[idx], 'b-value=', b_cutoff[idx], 'error b-value=', sigmacomp,
                      'a-value=', avaluecomp)
                x = np.linspace(x[idx], magmax, 1000)
                plt.semilogy(x, pow(10, avaluecomp - b_cutoff[idx] * x), linestyle='-', label='b-value')
            flag = 1
    # compare to normed hist output"
    plt.semilogy(avm, no_cumulative, '^r')
    fig21 = plt.gcf()
    fig21.set_size_inches(6.5, 4.0)
    cumulativerev = list(reversed(no_cumulative))
    d2 = list(reversed(np.cumsum(cumulativerev)))
    plt.semilogy(avm, d2, 'ok', linestyle='--')
    plt.xlabel('$Magnitude$')
    plt.ylabel('$Cumulative Number$')
    
    return b_cutoff, deltab_value


def b_value(magcat, mc):

    bin_m = 0.1
    magcat_new = ([mgev for mgev in magcat if mgev >= mc])
    n = len(magcat_new)
    diffmagcat = np.nansum(magcat_new - np.nanmean(magcat_new)) ** 2
    b = 1.0 / (log(10) * (np.nanmean(magcat_new) - (mc - bin_m / 2)))
    deltab = 2.3 * (b ** 2) * sqrt(diffmagcat / (n * (n - 1)))
    sigma = b / sqrt(n)
    avalue = log10(n) + b * mc
    return b, deltab, sigma, avalue, magcat_new


def bootstrap_mc(mag_data, alpha):
    # Estimates Mc uncertainty by bootstrap method

    print("Estimating Mc uncertainty")
    iterations = 50
    mc_bootstrap = []
    for i in range(iterations):
        boot = np.random.choice(mag_data, size=len(mag_data), replace=True)
        boot_mc = lilliefors(boot, alpha)
        mc_bootstrap.append(boot_mc)
    mc_sigma_bootstrap = np.std(mc_bootstrap)  # mc standard deviation

    return mc_sigma_bootstrap


def mc_vs_time(magcat, magcat_presequence, mc_ok, st_dev_multiplier, alpha, serial_times, size, step):
    # Analyses how the magnitude of completeness mc varies with time and identifies
    # critical regions ("STAI gaps"), where mc >= mc_ok + n * sigma, with n = st_dev_multiplier
    # mc_ok is the reference value for Mc estimated for the pre-sequence period

    register_matplotlib_converters()

    mc_ok_sigma = bootstrap_mc(magcat_presequence, alpha)
    print("Mc Sigma = ", mc_ok_sigma)
    upper_lim = mc_ok + st_dev_multiplier * mc_ok_sigma

    i = np.arange(0, len(magcat) - size, step)  # starting indexes of the windows

    # Moving window approach to estimate Mc vs time
    # Mc is estimated every "step" ( = 250 by default) events, at t_window (the end time of the moving window)
    print("----- Estimating Mc vs Time in the subcatalog -----")
    mc_time = []
    t_window = []  # time of the moving window (represented by the time of the last event within each window)
    for j in i:
        window = magcat[j: j + size]
        mc_new = lilliefors(window, alpha)
        mc_time.append(mc_new)
        t_window.append(serial_times[j + size])

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

    # Exclude too small gaps, which could simply arise from statistical fluctuations  of the magnitude of completeness
    hole_lower_lim , hole_upper_lim = [], []
    for j in range(len(tmp_hole_lower_lim)):
        index = [idx for idx, val in enumerate(serial_times) if
                 tmp_hole_lower_lim[j] <= val < tmp_hole_upper_lim[j]]
        if len(index) >= 2 * step:
            hole_lower_lim.append(tmp_hole_lower_lim[j])
            hole_upper_lim.append(tmp_hole_upper_lim[j])

    print('Found:', len(hole_lower_lim), 'STAI gaps')
    if len(hole_lower_lim) == 0:
        print("No STAI gaps found! Consider the following options:\n"
            "1. Decrease the fault_length_multiplier\n"
            "2. Decrease the st_dev_multiplier\n"
            "3. Check if the reference Mc is too high: either decrease the duration of the seismically quiescent period or set mc manually in the input file")
        sys.exit()

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
    ticks = equally_spaced_ticks(dates1)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(dynamic_date_formatter(dates1))
    ax.plot(datenums1, mc_time_plot, label='Mc(t)', ls='-')
    ax.fill_between(datenums1, upper_lim, mc_time_plot,
                    where=upper_lim <= mc_time_plot, color='C1', alpha=0.7,
                    label=f"$M_c \geq M^*_c$ + {st_dev_multiplier} $\sigma$")
    ax.set_xlabel("Time", labelpad=10, fontsize=14) 
    ax.set_ylabel("Mc", labelpad=10, fontsize=14)  
    ax.legend(loc='upper right')
    plt.savefig('fig/Mc_Time.pdf', format='pdf')
    
    return hole_lower_lim, hole_upper_lim, mc_times_hole, mc_hole


def map_plot(catalog, lon1, lat1, lon2, lat2):

    llcrnrlon, llcrnrlat = min(catalog[:, LonColumn]), min(catalog[:, LatColumn])
    urcrnrlon, urcrnrlat = max(catalog[:, LonColumn]), max(catalog[:, LatColumn])

    m = Basemap(projection='merc', resolution='i', llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)

    try:
        m.drawcoastlines(linewidth=0.5)
    except:
        pass
    m.drawcountries()
    m.drawmapboundary()

    lon_range = urcrnrlon - llcrnrlon
    lat_range = urcrnrlat - llcrnrlat
    lat_interval = round(0.2 * lat_range, 2)
    lon_interval = round(0.2 * lon_range, 2)
  
    parallels = np.arange(-90, 90, lat_interval)
    meridians = np.arange(-180, 180, lon_interval)
    m.drawparallels(parallels, labels=[1, 1, 1, 1], fontsize=10, dashes=[1, 5])
    m.drawmeridians(meridians, labels=[1, 1, 1, 1], fontsize=10, dashes=[1, 5])
    m.fillcontinents(color='#F0F0F0', lake_color='#F0F0F0')
    x1, y1 = m(lon1, lat1)
    x2, y2 = m(lon2, lat2)
    plt.plot(x1, y1, 'o', color='0.3', alpha=.7, markersize=5, markeredgecolor='w')
    plt.plot(x2, y2, 'o', color='steelblue', alpha=.7, markersize=5, markeredgecolor='w')


    return m

def spatial_smoothing(catalog, xmin, xmax, ymin, ymax, Sigma, sbin):
    # Smooth seismicity with Gaussian kernel

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
        R = np.arccos(A) * 6371
        W = (1 / (2 * pi * pow(Sigma, 2))) * (np.exp((-pow(R, 2)) / (2 * pow(Sigma, 2))))
        Smooth[:, 3] = Smooth[:, 3] + W

    TotalRate = sum(Smooth[:, 3])
    Smooth[:, 4] = Smooth[:, 3] / TotalRate

    lon_cell_grid = Smooth[:, 0]
    lan_cell_grid = Smooth[:, 1]
    smoothed_rate = Smooth[:, 4]

    return lon_cell_grid, lan_cell_grid, smoothed_rate


def trainv(lat1, long1, distx, disty):

    r0 = 6367.
    alfa = 57.29578
    fi1 = lat1 + alfa * disty / r0
    fr1 = fi1 / alfa
    lat2 = fi1 - alfa * distx * distx / (2. * r0 * r0) * math.tan(lat1 / alfa)
    lon2 = long1 + alfa * distx / (r0 * math.cos(fr1)) - alfa * distx ** 3 / (3. * r0 ** 3) * (
        math.tan(fr1)) ** 2 / math.cos(fr1)

    return lat2, lon2


def magnitude_distribution(magcat):

    min_m = min(magcat)
    max_m = max(magcat)
    bin_m = 0.1

    bin_edges = np.arange(min_m, max_m + bin_m, bin_m)
    bin_counts = np.zeros(len(bin_edges), dtype=int)
    for mag in magcat:
        bin = int(round(mag / bin_m - min_m*10))
        bin_counts[bin] += 1

    idx_zeros = np.where(bin_counts == 0)[0]
    bin_counts = np.delete(bin_counts, idx_zeros)
    bin_edges = np.delete(bin_edges, idx_zeros)

    log_bin_counts = [math.log10(i) for i in bin_counts]  # non-cumulative counts

    cum_counts = np.cumsum(bin_counts[::-1])[::-1]
    log_cum_counts = [math.log10(i) for i in cum_counts]  # cumulative counts

    return bin_edges, log_bin_counts, log_cum_counts


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
    
def select_distribution(user_input):
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
                  depth_distribution, p0):
    # Simulates magnitude, time, lat, lon and depth of ghost events

    magcat = catalog[:, MagnColumn]
    filtered_catalog, mask, magcat_presequence = filter_catalog(catalog, t_start, serial_times, fault_length_multiplier)

    # Save meaningful filtered dataset to file 
    write_filterd_catalog(filtered_catalog)

    if str(mc_user) == 'None':
        ref_mc = lilliefors(magcat_presequence, alpha)
    else:
        ref_mc = mc_user
    print("Reference Mc  = ", ref_mc)  
    if str(b_user) == 'None':  
        b, _, _, _, _ = b_value(catalog[:, MagnColumn], ref_mc)
    else:
        b = b_user
    print("b-value = ", b)

    magcat_filtered = filtered_catalog[:, MagnColumn]
    serial_times_filtered = np.array(serial_times)[mask]

    bin_m = 0.1

    # "hole" --> STAI gap
    hole_lower_lim, hole_upper_lim, mc_times_hole, mc_hole = mc_vs_time(magcat_filtered, magcat_presequence, ref_mc,
                                                                        st_dev_multiplier, alpha, serial_times_filtered, size, step)
    
    # Fit hypo depth distribution
    z_data = filtered_catalog[:, DepthColumn]
    bins = np.linspace(min(z_data), max(z_data), 101) 
    if depth_distribution == 'bimodal':
        y, x = np.histogram(z_data, bins=bins)
    else:
        y, x = np.histogram(z_data, density=1., bins=bins)
    x = (x[1:]+x[:-1])/2 # bin centers

    distribution_function = select_distribution(depth_distribution)
    params = fit_distribution(distribution_function, x, y, p0=p0)
    print('Fitted parameters = ', params)
    y, x, _ = plt.hist(z_data, bins)
    x_fit = np.linspace(x.min(), x.max(), 500)
    y_fit = distribution_function(x_fit, *params)

    print("----- Starting catalog filling -----")
    m_ghost, t_ghost, n_ghost_hole = [], [], []
    for i in range(len(hole_lower_lim)):  # iterate over all the STAI gaps

        n_ghost_step_list = []
        for j in range(len(mc_hole[i])):  # iterate over all the steps in the STAI gap

            mc_step = mc_hole[i][j]  # Mc value for a given step

            if j == 0:

                index = [idx for idx, val in enumerate(serial_times) if hole_lower_lim[i] < val <= mc_times_hole[i][j]]

            else:

                index = [idx for idx, val in enumerate(serial_times) if
                         mc_times_hole[i][j - 1] < val <= mc_times_hole[i][j]]

            nbin = int(mc_step * 10 - ref_mc * 10)  # total number of bins in the step

            # Compare expected and observed number of events for each bin in the step
            count = 0
            m2 = mc_step  # current upper bound of magnitude bin
            m1 = float(decimal.Decimal(mc_step) - decimal.Decimal(bin_m))  # current lower bound of magnitude bin

            # Expected number of magnitudes m >= mc_step (since magnitudes in the step are complete for m >= mc_step,
            # the EXPECTED number of magnitudes with m >= mc_step equals the OBSERVED number of magnitudes with m >= mc_step)
            n_m2 = len([magcat[i] for i in index if magcat[i] >= m2])

            n_ghost_bin_list = []
            while count < nbin:  # iterate over all the bins in the step

                # Expected number of magnitude m >= m1
                n_m1 = int(n_m2 * 10 ** (b * bin_m))
                # Expected number of magnitudes in the bin (i.e. with m == m1)
                n_exp_bin = n_m1 - n_m2

                # Observed number of magnitudes in the bin (i.e. with m == m1)
                n_obs_bin = len([magcat[i] for i in index if magcat[i] == m1])

                # Estimate the number of ghost magnitudes in the bin
                n_ghost_bin = n_exp_bin - n_obs_bin

                if n_ghost_bin <= 0:  # due to the inherent variability, bins can have more magnitude values than those expected!
                    # Update variables
                    m2 = m1
                    m1 = float(decimal.Decimal(m2) - decimal.Decimal(bin_m))
                    n_m2 = n_m1
                    count += 1
                    continue

                else:
                    # List of missing magnitudes
                    m_ghost_bin = [m1] * n_ghost_bin

                    # Update variables
                    m2 = m1
                    m1 = float(decimal.Decimal(m2) - decimal.Decimal(bin_m))
                    n_m2 = n_m1
                    count += 1

                    m_ghost.extend(m_ghost_bin)  # append all magnitudes to one list
                    n_ghost_bin_list.append(n_ghost_bin)  # from bin to step scale

            n_ghost_step = sum(n_ghost_bin_list)  # number of ghost events in the current step

            # Times random sampling
            for _ in range(n_ghost_step):

                u = np.random.random()

                # times are sampled from an uniform distribution whose support are the limits of the interval
                if j == 1:
                    random_t = hole_lower_lim[i] + u * (mc_times_hole[i][1] - hole_lower_lim[i])
                    # Add the events occurred after the large event to the first interval, i.e.
                    # the events occurred in the range [hole_lower_lim[i], mc_times_hole[i][0]]
                    # to those occurred in the first interval ( [mc_times_hole[i][0], mc_times_hole[i][1]] )
                else:
                    random_t = mc_times_hole[i][j - 1] + u * (mc_times_hole[i][j] - mc_times_hole[i][j - 1])
                    # Eqs occurred in the range [mc_times_hole[i][j-1], mc_times_hole[i][j]]

                t_ghost.append(random_t)  # append all occurrence times to one list

            n_ghost_step_list.append(n_ghost_step)  # from step to STAI gap scale
        print('Filled', i + 1, 'over a total of', len(hole_lower_lim), 'STAI gaps with', sum(n_ghost_step_list),
                      'events')

        n_ghost_hole.append(sum(n_ghost_step_list))  # append the number of ghost events in the current STAI gap

    # Simulate lat, lon and depth of ghost events
        
    xmin, ymin = min(filtered_catalog[:, LonColumn]), min(filtered_catalog[:, LatColumn])
    xmax, ymax = max(filtered_catalog[:, LonColumn]), max(filtered_catalog[:, LatColumn])
    lon_ghost, lat_ghost, z_ghost = [], [], []
    for i in range(len(hole_lower_lim)):  # iterate over all the STAI gaps
        lon_ghost_hole = []
        lat_ghost_hole = []

        index = [idx for idx, val in enumerate(serial_times_filtered) if hole_lower_lim[i] <= val <= hole_upper_lim[i]]
        subcatalog = filtered_catalog[index]
        lon_cell_grid, lat_cell_grid, smoothed_rate = spatial_smoothing(subcatalog, xmin, xmax, ymin, ymax, Sigma, sbin)

        # index = [idx for idx, val in enumerate(serial_times) if hole_lower_lim[i] <= val <= hole_upper_lim[i]]
        # subcatalog_tmp = catalog[index]
        # mag_big_shock = max(subcatalog_tmp[:, MagnColumn])
        # lat_big_shock = subcatalog_tmp[(subcatalog_tmp[:, MagnColumn] == mag_big_shock), LatColumn][0]
        # lon_big_shock = subcatalog_tmp[(subcatalog_tmp[:, MagnColumn] == mag_big_shock), LonColumn][0]

        # moment_big_shock = 10 ** (3 / 2 * (mag_big_shock + 10.7)) * 10 ** (-7)
        # l_big_shock = 10 ** (-5.20 + 0.35 * math.log10(moment_big_shock))  # rupture length (from Mai and Beroza (2000))

        # ymax, xmax = trainv(lat_big_shock, lon_big_shock, l_big_shock, l_big_shock)
        # ymin, xmin = trainv(lat_big_shock, lon_big_shock, -l_big_shock, -l_big_shock)

        # subcatalog = subcatalog_tmp[
        #              (subcatalog_tmp[:, LatColumn] >= ymin) & (subcatalog_tmp[:, LatColumn] <= ymax) &
        #              (subcatalog_tmp[:, LonColumn] >= xmin) & (subcatalog_tmp[:, LonColumn] <= xmax), :]

        # lon_cell_grid, lat_cell_grid, smoothed_rate = spatial_smoothing(subcatalog, xmin, xmax, ymin, ymax, Sigma, sbin)

        # Simulate lat and lon of ghost eqs
        idx = np.argsort(smoothed_rate)
        cumulative_sum = np.cumsum(smoothed_rate[idx])

        # Number of eqs to simulate
        n_to_sample = n_ghost_hole[i]

        # Depth
        samples = np.random.choice(x_fit, size=n_to_sample, p=y_fit / np.sum(y_fit))
        z_ghost.extend([round(value, 1) for value in samples])

        # Generate n=n_to_sample random numbers btw 0 and 1
        u = np.random.random(n_to_sample)

        for j in range(len(u)):
            if u[j] < cumulative_sum[0]:
                rand_lon = lon_cell_grid[idx[0]] - np.random.random() * sbin
                rand_lat = lat_cell_grid[idx[0]] - np.random.random() * sbin
            else:
                for k in range(1, len(smoothed_rate)):
                    if (cumulative_sum[k - 1] <= u[j]) and (u[j] < cumulative_sum[k]):
                        rand_lon = lon_cell_grid[idx[k - 1]] + np.random.random() * sbin
                        rand_lat = lat_cell_grid[idx[k - 1]] + np.random.random() * sbin
            lon_ghost_hole.append(round(rand_lon, 4))
            lat_ghost_hole.append(round(rand_lat, 4))

        lon_ghost.extend(lon_ghost_hole)
        lat_ghost.extend(lat_ghost_hole)

    # Plot hypo depth distribution

    data_min = np.min(z_ghost)
    data_max = np.max(z_ghost)
    margin = (data_max - data_min) * 0.95  
    depth_min, depth_max = data_min - margin, data_max + margin
    
    fig, ax = plt.subplots(figsize=(8, 6))
    if depth_distribution == 'bimodal':
        ax.hist(z_data, bins=bins, alpha=0.3, color='0.3', label='Hypo depths (filtered) original')
        ax.plot(x_fit, y_fit, color='red', label=f"Fitted distribution: {depth_distribution}")
        ax.hist(z_ghost, bins=bins, color='steelblue', label='Hypo depths replenished')
    else:
        ax.hist(z_data, bins=bins, density=True, alpha=0.3, color='0.3', label='Hypo depths (filtered) original')
        ax.plot(x_fit, y_fit, color='red', label=f"Fitted distribution: {depth_distribution}")
        ax.hist(z_ghost, bins=bins, density=True, color='steelblue', label='Hypo depths replenished')
    ax.set_xlim(depth_min, depth_max)
    ax.set_xlabel('Hypocenter Depth [km]', labelpad=10, fontsize=14)
    ax.legend()
    plt.savefig('fig/Hypo_Depth_Distribution.pdf', format='pdf')
    
    ####    

    # idx_magcatok = [idx for idx, val in enumerate(magcat) if val >= ref_mc]
    # magcat_ok = [magcat[i] for i in idx_magcatok]
    # serial_times_ok = [serial_times[i] for i in idx_magcatok]

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
    

    fig, ax = plt.subplots()
    ax.set_xlabel("Sequential number")
    ax.set_ylabel("Magnitude")
    ax.set_title('Replenished (filtered) catalog')
    plt.plot(x2, sorted_new_mag, 'o', markersize=3, markerfacecolor='0.3', markeredgecolor='None')
    plt.savefig('fig/Magnitude_SeqNumbers_Replenished.pdf', format='pdf')
    
    # Map
    # Plot original (m >= Mc) + simulated eqs

    lat = filtered_catalog[:, LatColumn]
    lon = filtered_catalog[:, LonColumn]
    lat_compl = [lat[i] for i in idx_magcatok_filtered]
    lon_compl = [lon[i] for i in idx_magcatok_filtered]

    fig, ax = plt.subplots()
    m = map_plot(filtered_catalog, lon_compl, lat_compl, lon_ghost, lat_ghost)
    plt.savefig('fig/Spatial_map.pdf', format='pdf')
    
    return ref_mc, b, m_ghost, t_ghost, lon_ghost, lat_ghost, z_ghost


def replenished_catalog(catalog, t_start, mc_user, b_user, size, step, st_dev_multiplier,
                        alpha, serial_times, Sigma, sbin, fault_length_multiplier,
                        depth_distribution, p0):
    # Returns the replenished catalog (original + ghost events)

    path = "fig"
    if not os.path.exists(path):
        os.mkdir(path)

    ref_mc, b, m_ghost, t_ghost, lon_ghost, lat_ghost, z_ghost = ghost_element(catalog, t_start, mc_user, b_user,
                                                           size, step, st_dev_multiplier, alpha, serial_times, Sigma, sbin, fault_length_multiplier,
                                                           depth_distribution, p0)

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

    depth = np.zeros(len(new_times))

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
    replenished_catalog[:, DayColumn] = sorted_new_month[:]
    replenished_catalog[:, MagnColumn] = sorted_new_mag[:]
    replenished_catalog[:, DepthColumn] = np.zeros(len(new_times))
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

    bin_edges_original, log_bin_counts_original, log_cum_counts_original = \
        magnitude_distribution(catalog[:, MagnColumn])
    bin_edges_replenished, log_bin_counts_replenished, log_cum_counts_replenished = \
        magnitude_distribution(replenished_catalog[:, MagnColumn])

    fig, ax = plt.subplots(figsize=(8, 10))
    ax.plot(bin_edges_original, log_bin_counts_original, 'o', color='0.5', markersize=3, label='Non-cumulative counts (Original)', alpha=0.7)
    ax.plot(bin_edges_original, log_cum_counts_original, '-', color='0.5', label='Cumulative counts (Original)', alpha=0.7)
    ax.plot(bin_edges_replenished, log_bin_counts_replenished, 'o', color='black', markersize=3, label='Non-cumulative counts (Replenished)', alpha=0.7)
    ax.plot(bin_edges_replenished, log_cum_counts_replenished, '-', color='black', label='Cumulative counts (Replenished)', alpha=0.7)
    ax.set_xlabel("Magnitude", labelpad=10, fontsize=14)
    ax.set_ylabel("Frequency", labelpad=10, fontsize=14)
    ax.legend(loc='upper right')
    plt.savefig('fig/Magnitude_Frequency_Distribution.pdf', format='pdf')
    
    return replenished_catalog
