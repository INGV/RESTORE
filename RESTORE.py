"""

RESTORE: REal catalogs STOchastic REplenishment
Written by Angela Stallone with help from Giuseppe Falcone
For any comment, question or suggestion write to:
angela.stallone@ingv.it

This project has been founded by the Seismic Hazard Center
(Centro di Pericolosità Sismica, CPS, at the Istituto Nazionale di Geofisica e Vulcanologia, INGV)

|---------------------------|

Latest revision: December, 2020

|---------------------------|

To cite:
Stallone A., Falcone G. 2020. Missing earthquake data reconstruction in the space-time-magnitude domain.
Preprint on https://essoar.org (2020) DOI: 10.1002/essoar.10504916.2

|--------------------------------------------------------------------------|

ABOUT

RESTORE simulates time, location and magnitude of those earthquakes
that have not been detected by the seismic network due to the overlap of earthquake signals in seismic records
after the occurrence of a large earthquake (short term aftershock incompleteness - STAI).
The algorithm assesses the temporal variability of the magnitude of completeness Mc
by means of a sliding overlapping windows approach, which collects estimates of
Mc at the end of each window. Since the window has a fixed number of events k and its
shift dk is constant, estimates of Mc are elapsed by dk events.
A statistic-based approach is implemented to pinpoint those time intervals where a threshold value for the
magnitude of completeness is significantly exceeded ("STAI gaps").
The number of missing events within each dk-step is estimated by calculating the difference between the
observed counts and the counts predicted by the Gutenberg-Richter relationship. Magnitude,
occurrence time and location of the simulated events are reconstructed implementing Monte
Carlo sampling techniques.

"""

# Column numbers in the saved catalog

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
    import csv
    import numpy as np
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
    from urllib.request import urlopen
    import json
    import dateparser
    import csv
    import numpy as np
    from numpy import flipud
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
    import xmltodict
    from urllib.request import urlopen
    import dateparser
    import csv
    import numpy as np
    from numpy import flipud
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


def serial_time(input_):
    # Converts the input date (string format) or date vectors (year, month, day, hour, minute, second) into timestamps
    # Origin time: Year 0, Month 0, Day 0 (same as Matlab)

    import datetime
    import numpy as np

    secs_per_day = 24.0 * 60.0 * 60.0

    # input: catalog
    if isinstance(input_, np.ndarray):
        catalog_sel = input_
        out1 = []
        for i in range(len(catalog_sel[:, YearColumn])):
            microseconds = round((catalog_sel[i, SecondColumn] - int(catalog_sel[i, SecondColumn])) * 1000000.0)
            dates = datetime.datetime(int(catalog_sel[i, YearColumn]), int(catalog_sel[i, MonthColumn]),
                                      int(catalog_sel[i, DayColumn]), int(catalog_sel[i, HourColumn]),
                                      int(catalog_sel[i, MinuteColumn]), int(catalog_sel[i, SecondColumn]),
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


def timestamp_to_datetime(timestamp):
    # Converts custom timestamps to datetime objects

    import datetime

    days = timestamp % 1
    hours = days % 1 * 24
    minutes = hours % 1 * 60
    seconds = minutes % 1 * 60

    date_out = datetime.datetime.fromordinal(int(timestamp)) + datetime.timedelta(days=int(days)) + \
               datetime.timedelta(hours=int(hours)) + datetime.timedelta(minutes=int(minutes)) + \
               datetime.timedelta(seconds=seconds) - datetime.timedelta(days=366)

    return date_out


def lilliefors(magcat, mmin, alpha):
    # Estimates magnitude of completeness with Lilliefors test (Clauset et al., 2009)
    from statsmodels.stats.diagnostic import lilliefors as lillie
    from numpy import random

    bin_m = 0.1
    incert = ((random.randint(0, 10001, size=len(magcat))) - 5000) / 100000
    mag = magcat[:] + incert[:]
    upperlim = 40

    for i in range(int(mmin * 10), upperlim, 1):
        lowlim = float(i / 10) - bin_m / 2
        magsel = ([mgev for mgev in mag if mgev > lowlim])
        magselmin = [x - lowlim for x in magsel]
        kstest, pvalue_lilli = lillie(magselmin, dist='exp')
        if pvalue_lilli < alpha:
            continue
        else:
            break
    mc = float(i / 10)
    pval = pvalue_lilli

    return mc, pval


def b_value_zmap(magcat):
    # Calculate the b_value, tha a_value and the completeness magnitude using the ZMAP approach

    import numpy as np
    import matplotlib.pyplot as plt
    from math import sqrt, log, log10

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
    plt.show()
    return b_cutoff, deltab_value


def b_value(magcat, mc):
    from math import sqrt, log, log10
    import numpy as np
    bin_m = 0.1
    magcat_new = ([mgev for mgev in magcat if mgev >= mc])
    n = len(magcat_new)
    diffmagcat = np.nansum(magcat_new - np.nanmean(magcat_new)) ** 2
    b = 1.0 / (log(10) * (np.nanmean(magcat_new) - (mc - bin_m / 2)))
    deltab = 2.3 * (b ** 2) * sqrt(diffmagcat / (n * (n - 1)))
    sigma = b / sqrt(n)
    avalue = log10(n) + b * mc
    return b, deltab, sigma, avalue, magcat_new


def bootstrap_mc(mag_data, mmin, alpha):
    # Estimates uncertainty on magnitude of completeness mc by bootstrap method

    import numpy as np

    iterations = 200
    mc_bootstrap = []
    for _ in range(iterations):
        boot = np.random.choice(mag_data, size=len(mag_data), replace=True)
        boot_mc, _ = lilliefors(boot, mmin, alpha)
        mc_bootstrap.append(boot_mc)
    mc_sigma_bootstrap = np.std(mc_bootstrap)  # mc standard deviation

    return mc_sigma_bootstrap


def mc_vs_time(magcat, magcat_presequence, mmin, mc_ok, alpha, serial_times, size, step):
    # Analyses how the magnitude of completeness mc varies with time and identifies
    # critical regions ("STAI gaps"), where mc >= mc_ok + 2 * sigma
    # mc_ok is the reference value for mc estimated for the pre-sequence period

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from pandas.plotting import register_matplotlib_converters
    register_matplotlib_converters()

    mc_ok_sigma = bootstrap_mc(magcat_presequence, mmin, alpha)
    upper_lim = mc_ok + 2 * mc_ok_sigma

    i = np.arange(0, len(magcat) - size, step)  # starting indexes of the windows

    # Moving window approach to estimate mc vs time
    # mc is estimated every "step" ( = 250 by default) events, at t_window (the end time of the moving window)
    mc_time = []
    t_window = []  # time of the moving window (represented by the time of the last event within each window)
    for j in i:
        window = magcat[j: j + size]
        mc_new, _ = lilliefors(window, mmin, alpha)
        mc_time.append(mc_new)
        t_window.append(serial_times[j + size])

    # Find the temporal bounds of the STAI gaps (where mc >= mc_ok + 2 * sigma)
    tmp_hole_lower_lim = []
    tmp_hole_upper_lim = []
    for i in range(len(t_window) - 1):

        # If these conditions are met, mc is exceeding the upper limit (mc >= mc_ok + 2 * sigma)
        # --> STAI gap starts
        if i == 0:
            if mc_time[i] >= upper_lim:
                # Set the STAI gap starting time to the time of the first event
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

    # Extract mc and relative times within each STAI gap
    mc_times_hole, mc_hole = [], []
    for j in range(len(hole_lower_lim)):
        index = [idx for idx, val in enumerate(t_window) if hole_lower_lim[j] <= val < hole_upper_lim[j]]
        tmp_mc_times_hole = [t_window[i] for i in index]
        tmp_mc_hole = [mc_time[i] for i in index]
        mc_times_hole.append(tmp_mc_times_hole)
        mc_hole.append(tmp_mc_hole)

    # plot mc vs time
    t_window_plot = [val for idx, val in enumerate(t_window) if val >= hole_lower_lim[0]]
    idx_t_window_plot = [idx for idx, val in enumerate(t_window) if val >= hole_lower_lim[0]]
    dates1 = [timestamp_to_datetime(t_window_plot[i]) for i in range(len(t_window_plot))]
    datenums1 = mdates.date2num(dates1)
    mc_time_plot = [mc_time[i] for i in idx_t_window_plot]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.xaxis_date()
    fmt = mdates.DateFormatter('%b %Y')
    ax.xaxis.set_major_formatter(fmt)
    ax.plot(datenums1, mc_time_plot, label='Mc(t)', ls='-')
    ax.fill_between(datenums1, upper_lim, mc_time_plot,
                    where=upper_lim <= mc_time_plot, color='C1', alpha=0.7, label='$ M_{c} \geq M^{*}_{c} + 2 \sigma$')
    ax.set_xlabel("Time")
    ax.set_ylabel("Mc")
    ax.legend(loc='upper right')
    plt.savefig('fig/Mc_Time.pdf', format='pdf')

    return hole_lower_lim, hole_upper_lim, mc_times_hole, mc_hole


def map_plot(lon1, lat1, lon2, lat2, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat):
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import numpy as np

    m = Basemap(projection='merc', resolution='i', llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)

    try:
        m.drawcoastlines(linewidth=0.5)
    except:
        pass
    m.drawcountries()
    m.drawmapboundary()
    parallels = np.arange(-90, 90, 0.5)
    meridians = np.arange(-180, 180, 0.5)
    m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=10, dashes=[1, 5])
    m.drawmeridians(meridians, labels=[0, 0, 1, 0], fontsize=10, dashes=[1, 5], rotation=45)

    x1, y1 = m(lon1, lat1)
    x2, y2 = m(lon2, lat2)

    plt.plot(x1, y1, 'o', color='0.3', alpha=.7, markersize=5, markeredgecolor='w')
    plt.plot(x2, y2, 'o', color='steelblue', alpha=.7, markersize=5, markeredgecolor='w')
    plt.xlabel('Longitude', fontsize=12, labelpad=15)
    plt.ylabel('Latitude', fontsize=12, labelpad=40)

    return m


def spatial_smoothing(subcatalog, xmin, xmax, ymin, ymax):
    # Smooth seismicity with Gaussian kernel

    import numpy as np
    pi = np.pi
    sbin = 0.01
    sbinx = sbin
    sbiny = sbin
    Sigma = 1

    rows_cat, columns_cat = np.shape(subcatalog)
    X = np.arange(xmin, xmax + sbinx, sbinx)
    Y = np.arange(ymin, ymax + sbiny, sbiny)
    xi, yi = np.meshgrid(X, Y)
    Smooth = np.zeros((len(X) * len(Y), 5))
    Smooth[:, 0] = np.reshape(xi, (len(X) * len(Y)))
    Smooth[:, 1] = np.reshape(yi, (len(X) * len(Y)))
    onetoten = range(0, rows_cat)
    for i in onetoten:
        A = (np.sin(Smooth[:, 1] * pi / 180) * np.sin(float(subcatalog[i, 1]) * pi / 180) + np.cos(
            Smooth[:, 1] * pi / 180) * np.cos(float(subcatalog[i, 1]) * pi / 180) * np.cos(
            ((float(subcatalog[i, 0])) - Smooth[:, 0]) * pi / 180))
        R = np.arccos(A) * 6371
        W = (1 / (2 * pi * pow(Sigma, 2))) * (np.exp((-pow(R, 2)) / (2 * pow(Sigma, 2))))
        Smooth[:, 3] = Smooth[:, 3] + W

    TotalRate = sum(Smooth[:, 3])
    Smooth[:, 4] = Smooth[:, 3] / TotalRate

    lon_cell_grid = Smooth[:, 0]
    lan_cell_grid = Smooth[:, 1]
    smoothed_rate = Smooth[:, 4]

    return lon_cell_grid, lan_cell_grid, smoothed_rate, sbin


def trainv(lat1, long1, distx, disty):
    import math

    r0 = 6367.
    alfa = 57.29578
    fi1 = lat1 + alfa * disty / r0
    fr1 = fi1 / alfa
    lat2 = fi1 - alfa * distx * distx / (2. * r0 * r0) * math.tan(lat1 / alfa)
    lon2 = long1 + alfa * distx / (r0 * math.cos(fr1)) - alfa * distx ** 3 / (3. * r0 ** 3) * (
        math.tan(fr1)) ** 2 / math.cos(fr1)

    return lat2, lon2


def magnitude_distribution(magcat):
    import numpy as np
    import math

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


def ghost_element(catalog_sel, mmin, mc_ok, magcat, magcat_presequence, b,
                  size, step, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
                  alpha, serial_times):
    # Simulates magnitude, time, lat and lon of ghost events

    from numpy import random
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import math
    import numpy as np

    bin_m = 0.1

    # "hole" --> STAI gap
    hole_lower_lim, hole_upper_lim, mc_times_hole, mc_hole = mc_vs_time(magcat, magcat_presequence, mmin, mc_ok,
                                                                        alpha, serial_times, size, step)

    m_ghost, t_ghost, n_ghost_hole = [], [], []
    for i in range(len(hole_lower_lim)):  # iterate over all the STAI gaps
        n_ghost_tmp = []
        for j in range(1, len(mc_hole[i])):  # iterate over all the steps in the STAI gap

            mc_step = mc_hole[i][j]  # mc value for a given step

            if j == 1:

                index = [idx for idx, val in enumerate(serial_times) if hole_lower_lim[i] < val <= mc_times_hole[i][1]]

            else:

                index = [idx for idx, val in enumerate(serial_times) if
                         mc_times_hole[i][j - 1] < val <= mc_times_hole[i][j]]

            n_obs = len([magcat[i] for i in index if magcat[i] >= mc_ok])

            n_mc_step = len([magcat[i] for i in index if magcat[i] >= mc_step])
            nbin = int(mc_step * 10 - mc_ok * 10)
            n_exp = int(n_mc_step * 10 ** (b * nbin * bin_m))

            n_ghost = n_exp - n_obs

            if n_ghost <= 0:
                continue

            else:
                n_ghost_tmp.append(n_ghost)

                # Magnitudes random sampling (from Gutenberg-Richter law)
                count = 1
                while count <= n_ghost:

                    u = random.random()

                    random_m = round((- 1/b * math.log10(u) + mc_ok - bin_m/2) * 10)/10

                    if random_m <= mc_step:

                        m_ghost.append(random_m)

                        count += 1

                    else:
                        count += 0

                # Times random sampling
                for _ in range(n_ghost):

                    u = random.random()

                    # times are sampled from an uniform distribution whose support are the limits of the interval
                    if j == 1:
                        random_t = hole_lower_lim[i] + u * (mc_times_hole[i][1] - hole_lower_lim[i])
                        # Add the events occurred after the large event to the first interval, i.e.
                        # the events occurred in the range [hole_lower_lim[i], mc_times_hole[i][0]]
                        # to those occurred in the first interval ( [mc_times_hole[i][0], mc_times_hole[i][1]] )
                    else:
                        random_t = mc_times_hole[i][j - 1] + u * (mc_times_hole[i][j] - mc_times_hole[i][j - 1])
                        # Eqs occurred in the range [mc_times_hole[i][j-1], mc_times_hole[i][j]]
                    t_ghost.append(random_t)

        n_ghost_hole.append(sum(n_ghost_tmp))

        print('Filled', i+1, 'over a total of', len(hole_lower_lim), 'STAI gaps with', sum(n_ghost_tmp), 'events')

    # Simulate lat e lon of ghost events
    lon_ghost, lat_ghost = [], []
    for i in range(len(hole_lower_lim)):  # iterate over all the STAI gaps
        lon_ghost_hole = []
        lat_ghost_hole = []
        index = [idx for idx, val in enumerate(serial_times) if hole_lower_lim[i] <= val <= hole_upper_lim[i]]
        subcatalog_tmp = catalog_sel[index]
        mag_big_shock = max(subcatalog_tmp[:, MagnColumn])
        lat_big_shock = subcatalog_tmp[(subcatalog_tmp[:, MagnColumn] == mag_big_shock), LatColumn][0]
        lon_big_shock = subcatalog_tmp[(subcatalog_tmp[:, MagnColumn] == mag_big_shock), LonColumn][0]

        moment_big_shock = 10 ** (3 / 2 * (mag_big_shock + 10.7)) * 10 ** (-7)
        l_big_shock = 10 ** (-5.20 + 0.35 * math.log10(moment_big_shock))  # rupture length (from Mai and Beroza (2000))

        ymax, xmax = trainv(lat_big_shock, lon_big_shock, l_big_shock, l_big_shock)
        ymin, xmin = trainv(lat_big_shock, lon_big_shock, -l_big_shock, -l_big_shock)

        subcatalog = subcatalog_tmp[
                     (subcatalog_tmp[:, LatColumn] >= ymin) & (subcatalog_tmp[:, LatColumn] <= ymax) &
                     (subcatalog_tmp[:, LonColumn] >= xmin) & (subcatalog_tmp[:, LonColumn] <= xmax), :]

        lon_cell_grid, lat_cell_grid, smoothed_rate, sbin = spatial_smoothing(subcatalog, xmin, xmax, ymin, ymax)

        # Simulate lat and lon of ghost eqs
        idx = np.argsort(smoothed_rate)
        cumulative_sum = np.cumsum(smoothed_rate[idx])

        # Generate n=n_ghost_hole[i] random numbers btw 0 and 1
        u = random.random(n_ghost_hole[i])

        for j in range(len(u)):
            if u[j] < cumulative_sum[0]:
                rand_lon = lon_cell_grid[idx[0]] - random.random() * sbin
                rand_lat = lat_cell_grid[idx[0]] - random.random() * sbin
            else:
                for k in range(1, len(smoothed_rate)):
                    if (cumulative_sum[k - 1] <= u[j]) and (u[j] < cumulative_sum[k]):
                        rand_lon = lon_cell_grid[idx[k - 1]] + random.random() * sbin
                        rand_lat = lat_cell_grid[idx[k - 1]] + random.random() * sbin

            lon_ghost_hole.append(round(rand_lon, 2))
            lat_ghost_hole.append(round(rand_lat, 2))

        lon_ghost.extend(lon_ghost_hole)
        lat_ghost.extend(lat_ghost_hole)

    idx_magcatok = [idx for idx, val in enumerate(magcat) if val >= mc_ok]
    magcat_ok = [magcat[i] for i in idx_magcatok]
    serial_times_ok = [serial_times[i] for i in idx_magcatok]

    # plots
    # Plot original (m >= mc) + simulated eqs

    dates1 = [timestamp_to_datetime(serial_times_ok[i]) for i in range(len(serial_times_ok))]
    dates2 = [timestamp_to_datetime(t_ghost[i]) for i in range(len(t_ghost))]
    datenums1 = mdates.date2num(dates1)
    datenums2 = mdates.date2num(dates2)

    fig, ax = plt.subplots(2, figsize=(10, 6))
    fmt = mdates.DateFormatter('%b %Y')
    ax[0].plot(datenums1, magcat_ok, 'o', markersize=3, markerfacecolor='0.3', markeredgecolor='None',
               label='Original data')
    ax[0].plot(datenums2, m_ghost, 'o', markersize=3, markerfacecolor='steelblue', markeredgecolor='None',
               label='Replenished data')
    ax[0].xaxis_date()
    ax[0].xaxis.set_major_formatter(fmt)
    ax[0].set_xlabel("Time")
    ax[0].set_ylabel("Magnitude")
    ax[0].legend(loc='upper left')
    ax[1].plot(datenums1, magcat_ok, 'o', markersize=3, markerfacecolor='0.3', markeredgecolor='None',
               label='Original data')
    ax[1].xaxis_date()
    ax[1].xaxis.set_major_formatter(fmt)
    ax[1].set_xlabel("Time")
    ax[1].set_ylabel("Magnitude")
    ax[1].legend(loc='upper left')
    plt.subplots_adjust(hspace=0.5)
    plt.savefig('fig/Magnitude_Time.pdf', format='pdf')

    # Map
    # Plot original (m >= mc) + simulated eqs

    lat = catalog_sel[:, LatColumn]
    lon = catalog_sel[:, LonColumn]
    lat_compl = [lat[i] for i in idx_magcatok]
    lon_compl = [lon[i] for i in idx_magcatok]

    fig, ax = plt.subplots()
    m = map_plot(lon_compl, lat_compl, lon_ghost, lat_ghost, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat)
    plt.savefig('fig/Spatial_map.pdf', format='pdf')

    return m_ghost, t_ghost, lon_ghost, lat_ghost


def replenished_catalog(catalog_sel, magcat, magcat_presequence, mmin, mc_ok, b, size, step,
                       llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, alpha, serial_times):
    # Returns the replenished catalog (original + ghost events)

    import numpy as np
    import csv
    import matplotlib.pyplot as plt
    import os

    path = "fig"
    if not os.path.exists(path):
        os.mkdir(path)

    m_ghost, t_ghost, lon_ghost, lat_ghost = ghost_element(catalog_sel, mmin, mc_ok, magcat, magcat_presequence, b,
                                                           size, step, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
                                                           alpha, serial_times)

    magcat_compl_idx = [idx for idx, val in enumerate(magcat) if val >= mc_ok]

    magcat_compl = [magcat[i] for i in magcat_compl_idx]
    serial_times_compl = [serial_times[i] for i in magcat_compl_idx]

    new_times = np.array(serial_times_compl + t_ghost)
    sorted_new_times = new_times[np.argsort(new_times)]

    new_mag = np.array(magcat_compl + m_ghost)
    sorted_new_mag = new_mag[np.argsort(new_times)]

    lat = catalog_sel[:, LatColumn]
    lon = catalog_sel[:, LonColumn]
    lat_compl = [lat[i] for i in magcat_compl_idx]
    lon_compl = [lon[i] for i in magcat_compl_idx]

    new_lat = np.array(lat_compl + lat_ghost)
    sorted_new_lat = new_lat[np.argsort(new_times)]

    new_lon = np.array(lon_compl + lon_ghost)
    sorted_new_lon = new_lon[np.argsort(new_times)]

    # flag: 0 for original events, 1 for simulated events
    flag = np.array([0] * len(new_times))
    for i in range(len(new_times)):
        if i <= len(serial_times_compl):
            flag[i] = 0
        else:
            flag[i] = 1

    sorted_flag = flag[np.argsort(new_times)]

    # Magnitude-Time plots

    x1 = np.arange(1, len(magcat_compl) + 1)
    x2 = np.arange(1, len(sorted_new_mag) + 1)

    fig, ax = plt.subplots()
    ax.set_xlabel("Sequential number")
    ax.set_ylabel("Magnitude")
    ax.set_title('Original catalog')
    plt.plot(x1, magcat_compl, 'o', markersize=3, markerfacecolor='0.3', markeredgecolor='None')
    plt.savefig('fig/Magnitude_SeqNumbers_Original.pdf', format='pdf')

    fig, ax = plt.subplots()
    ax.set_xlabel("Sequential number")
    ax.set_ylabel("Magnitude")
    ax.set_title('Replenished catalog')
    plt.plot(x2, sorted_new_mag, 'o', markersize=3, markerfacecolor='0.3', markeredgecolor='None')
    plt.savefig('fig/Magnitude_SeqNumbers_Replenished.pdf', format='pdf')

    depth = np.zeros(len(new_times))

    sorted_new_year, sorted_new_month, sorted_new_day, sorted_new_hour, sorted_new_min, sorted_new_sec = \
        [], [], [], [], [], [],

    for i in range(len(new_times)):
        date_out = timestamp_to_datetime(sorted_new_times[i])
        sorted_new_year.append(int(date_out.year))
        sorted_new_month.append(int(date_out.month))
        sorted_new_day.append(int(date_out.day))
        sorted_new_hour.append(int(date_out.hour))
        sorted_new_min.append(int(date_out.minute))
        sec = float(date_out.second) + (float(date_out.microsecond) / 1000000)
        sorted_new_sec.append(round(sec, 2))

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
               sorted_new_mag, depth, sorted_new_hour, sorted_new_min, sorted_new_sec, sorted_flag)
    with open('Replenished_catalog.txt', 'w') as f:
        writer = csv.writer(f, delimiter=' ', skipinitialspace=False, quoting=csv.QUOTE_NONE, lineterminator='\n')
        for x in rows:
            writer.writerow(x)

    bin_edges_original, log_bin_counts_original, log_cum_counts_original = \
        magnitude_distribution(catalog_sel[:, MagnColumn])
    bin_edges_replenished, log_bin_counts_replenished, log_cum_counts_replenished = \
        magnitude_distribution(replenished_catalog[:, MagnColumn])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 6))
    ax1.plot(bin_edges_original, log_bin_counts_original, 'o', color='0.5', markersize=3, label='Non-cumulative counts')
    ax1.plot(bin_edges_original, log_cum_counts_original, 'ok', markersize=3, label='Cumulative counts')
    ax1.legend(loc='upper right')
    ax1.set_xlabel("Magnitude")
    ax1.set_ylabel("Frequency")
    ax1.set_title("Original Catalog")
    ax2.plot(bin_edges_replenished, log_bin_counts_replenished, 'o', color='0.5', markersize=3, label='Non-cumulative counts')
    ax2.plot(bin_edges_replenished, log_cum_counts_replenished, 'ok', markersize=3, label='Cumulative counts')
    ax2.legend(loc='upper right')
    ax2.set_xlabel("Magnitude")
    ax2.set_ylabel("Frequency")
    ax2.set_title("Replenished Catalog")
    plt.savefig('fig/Magnitude_Frequency_Distribution.pdf', format='pdf')

    return replenished_catalog
