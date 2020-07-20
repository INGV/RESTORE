"""
This function simulates time, location and magnitude of those earthquakes
that have not been detected by the seismic network due to the decrease of signal-to-noise ratio
after the occurrence of a (relatively) large earthquake (STAI issue).

The algorithm assesses the temporal variability of the magnitude of completeness "mc"
by means of a sliding overlapping windows approach, which collects estimates of
"mc" at the end of each window. Since the window has a fixed number of events k and its
shift dk is constant, estimates of mc are elapsed by dk events.

A statistic-based approach is implemented to pinpoint those time intervals where a threshold value for the
magnitude of completeness is significantly exceeded ("STAI gaps").

The number of missing events within each dk-step is estimated by calculating the difference between the
observed counts and the counts predicted by the Gutenberg-Richter relationship. Magnitude,
occurrence time and location of the simulated events are reconstructed implementing Monte
Carlo sampling techniques.

|---------------------------|
 Latest revision: July, 2020
|---------------------------|


|--------------------------------------------|
 AUTHORS:
 Angela Stallone  (angela.stallone@ingv.it)
 Giuseppe Falcone (giuseppe.falcone@ingv.it)
|--------------------------------------------|


COPYRIGHT:
Stallone, A.; Falcone, G. (2020) "Missing earthquake data reconstruction in the space-time-magnitude domain."
Submitted to JGR Solid Earth

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


def read_catalog_cnt(name_catalog, depthm, magm, xmin, xmax, ymin, ymax):
    import csv
    import math
    import numpy as np
    from numpy import flipud
    import dateparser
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
    with open(name_catalog) as csvfile:
        read_csv = csv.reader(csvfile, skipinitialspace=True, delimiter='|')
        next(read_csv)
        i = 0
        for row in read_csv:
            if float(row[4]) <= depthm and xmin <= float(row[3]) <= xmax and ymin <= float(row[2]) <= ymax:
                if row[9] == 'Md':
                    magnitude = 1.612 * float(row[10]) - 1.633
                if row[9] == 'Mw':
                    magnitude = 0.938 * float(row[10]) + 0.154
                if row[9] == 'ML':
                    magnitude = float(row[10])
                if magnitude >= magm:
                    i = i + 1
                    time = dateparser.parse((row[1]))
                    lon.append(float(row[3]))
                    lat.append(float(row[2]))
                    year.append(int(time.year))
                    month.append(int(time.month))
                    day.append(int(time.day))
                    mag.append(round(magnitude, 1))
                    depth.append(float(row[4]))
                    hour.append(int(time.hour))
                    minute.append(int(time.minute))
                    second = float(time.second) + (float(time.microsecond) / 1000000)
                    sec.append(round(second, 2))
                    type_mag.append((row[9]))
            sec = [0 if math.isnan(x) else x for x in sec]
            minute = [0 if math.isnan(x) else x for x in minute]
            hour = [0 if math.isnan(x) else x for x in hour]
            day = [1 if math.isnan(x) else x for x in day]
            month = [1 if math.isnan(x) else x for x in month]
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
    with open('Zmap_catalog', 'w') as f:
        writer = csv.writer(f, delimiter=' ', skipinitialspace=True, quoting=csv.QUOTE_NONE, lineterminator='\n')
        for x in rows:
            writer.writerow(x)
    return catalog


def read_catalog_zmap(name_catalog):
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
        read_csv = csv.reader(csvfile, skipinitialspace=False, delimiter=',')
        for row in read_csv:
            lon.append(float(row[0]))
            lat.append(float(row[1]))
            year.append(int(row[2]))
            month.append(int(row[3]))
            day.append(int(row[4]))
            mag.append(round(float(row[5]), 1))
            depth.append(float(row[6]))
            hour.append(int(row[7]))
            minute.append(int(row[8]))
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
    rows = zip(lon, lat, year, month, day, mag, depth, hour, minute, sec)
    with open('Zmap_catalog', 'w') as f:
        writer = csv.writer(f, delimiter=' ', skipinitialspace=False, quoting=csv.QUOTE_NONE, lineterminator='\n')
        for x in rows:
            writer.writerow(x)
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
    # Converts the input date (string format) or date vectors
    # (year, month, day, hour, minute, second) into timestamps

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


def lilliefors(magcat, mmin):
    # Estimates magnitude of completeness with Lilliefors test (Clauset et al., 2009)

    from statsmodels.stats.diagnostic import lilliefors as lillie
    from numpy import random

    incert = ((random.randint(0, 10000, size=len(magcat))) - 5000) / 100000
    mag = magcat[:] + incert[:]
    bin_m = 0.1
    upperlim = max(magcat)

    m, pval, h = [], [], []
    count = 0
    for i in range(int(mmin * 10), int(upperlim * 10), 1):
        lowlim = float(i / 10) - bin_m / 2
        magsel = ([mgev for mgev in mag if mgev >= lowlim])
        magselmin = [x - lowlim for x in magsel]
        kstest, pvalue_lilli = lillie(magselmin, dist='exp')
        if pvalue_lilli >= 0.05:
            h.append(1)
        else:
            h.append(0)
        m.append(float(i / 10))
        pval.append(pvalue_lilli)
        count += 1

        if i == int(mmin * 10):
            continue
        else:
            if h[count - 2] == 1:
                break
            else:
                continue
    mc = m[count - 2]
    pvalue = pval[count - 2]
    return mc, pvalue


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
    differenza = np.zeros(shape=(len(x)))
    somma = np.zeros(shape=(len(x)))
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
        differenza[idx] = abs(b_ave[idx] - b_cutoff[idx])
        somma[idx] = b_ave[idx] + b_cutoff[idx]
        if differenza[idx] <= deltab_value[idx] and abs(dbi_old[idx]) <= 0.03:
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


def bootstrap_mc(mag_data, mmin):
    # Estimates uncertainty on magnitude of completeness mc by bootstrap method

    import numpy as np

    iterations = 100
    mc_bootstrap = []
    for _ in range(iterations):
        boot = np.random.choice(mag_data, size=len(mag_data), replace=True)
        boot_mc, _ = lilliefors(boot, mmin)
        mc_bootstrap.append(boot_mc)
    mc_sigma_bootstrap = np.std(mc_bootstrap)  # mc standard deviation

    return mc_sigma_bootstrap


def mc_vs_time(magcat, mmin, mc_ok, serial_times, size, step):
    # Analyses how the magnitude of completness mc varies with time and identifies
    # critical regions ("STAI gaps"), where mc > mc_ok + sigma

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from pandas.plotting import register_matplotlib_converters
    register_matplotlib_converters()

    mc_ok_sigma = bootstrap_mc(magcat, mmin)
    upper_lim = mc_ok + mc_ok_sigma

    i = np.arange(0, len(magcat) - size, step)  # starting indexes of the windows

    # Moving window approach to estimate mc vs time
    # mc is estimated every "step" ( = 250 by default) events, at t_window (the end time of the moving window)
    mc_time = []
    t_window = []  # time of the moving window (represented by the time of the last event within each window)
    for j in i:
        window = magcat[j: j + size]
        mc_new, _ = lilliefors(window, mmin)
        mc_time.append(mc_new)
        t_window.append(serial_times[j + size])

    # Temporal bounds of the STAI gaps (where mc > mc_ok + sigma)
    tmp_hole_lower_lim = []
    tmp_hole_upper_lim = []
    for i in range(1, len(t_window) - 1):
        if (i == 1 and mc_time[i] > upper_lim) or (mc_time[i] > upper_lim and mc_time[i - 1] <= upper_lim):
            # if these conditions are met, mc is exiting the upper limit (mc > mc_ok + sigma)
            # --> STAI gap starts
            tmp_hole_lower_lim.append(t_window[i])
        if mc_time[i] > upper_lim and mc_time[i + 1] <= upper_lim:
            # if these conditions are met, mc is re-entering in the range mc_ok + sigma
            # --> STAI gap ends
            tmp_hole_upper_lim.append(t_window[i])

    # Refine STAI gap start time by setting it to the time of the largest shock in the interval,
    # which has caused the raise of the magnitude of completeness
    for j in range(len(tmp_hole_lower_lim)):
        idx_critical = [idx for idx, val in enumerate(serial_times) if
                        (tmp_hole_lower_lim[j] - step) <= val <= tmp_hole_lower_lim[j]]
        val = magcat[idx_critical]
        iidx_main = [j for j, v in enumerate(val) if v == max(val)]
        idx_main = idx_critical[iidx_main[0]]
        tmp_hole_lower_lim[j] = serial_times[idx_main]

    # Exclude too small regions, which could simply arise from statistical fluctuations  of the magnitude of completeness
    a, b = [], []
    for j in range(len(tmp_hole_lower_lim)):
        index = [idx for idx, val in enumerate(serial_times) if
                 tmp_hole_lower_lim[j] <= val < tmp_hole_upper_lim[j]]
        if len(index) > 2 * step:
            a.append(tmp_hole_lower_lim[j])
            b.append(tmp_hole_upper_lim[j])
    hole_lower_lim = a
    hole_upper_lim = b

    # Extract mc and relative times within each STAI gap
    mc_times_hole, mc_hole = [], []
    for j in range(len(hole_lower_lim)):
        index = [idx for idx, val in enumerate(t_window) if hole_lower_lim[j] <= val < hole_upper_lim[j]]
        tmp_mc_times_hole = [t_window[i] for i in index]
        tmp_mc_hole = [mc_time[i] for i in index]
        mc_times_hole.append(tmp_mc_times_hole)
        mc_hole.append(tmp_mc_hole)

    # plot mc vs time
    fig, ax = plt.subplots()
    ax.xaxis_date()
    locator = mdates.MonthLocator()  # every month
    fmt = mdates.DateFormatter('%b')
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(fmt)
    ax.plot(t_window, mc_time, label='Mc vs Time', ls='-')
    ax.fill_between(t_window, upper_lim, mc_time,
                    where=upper_lim < mc_time, color='C1', alpha=0.7, label='Mc > sigma')
    x_holes = hole_lower_lim + hole_upper_lim
    for p in x_holes:
        plt.axvline(x=p, ls='--', linewidth=1, color='C2')
    ax.set_xlabel("Time")
    ax.set_ylabel("Mc")
    ax.legend(loc='upper right')
    plt.savefig('fig/Mc_Time.pdf', format='pdf')

    return hole_lower_lim, hole_upper_lim, mc_times_hole, mc_hole


def map_plot(lon1, lat1, lon2, lat2, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat):
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import numpy as np

    m = Basemap(projection='merc', lat_0=5, lon_0=5,
                resolution='i', area_thresh=0.1,
                llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)

    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary()
    parallels = np.arange(35, 48, 2)
    m.drawparallels(parallels, labels=[True, False, False, True], fontsize=10)
    meridians = np.arange(5, 21, 2)
    m.drawmeridians(meridians, labels=[True, False, False, True], fontsize=10)

    x1, y1 = m(lon1, lat1)
    x2, y2 = m(lon2, lat2)

    plt.plot(x1, y1, 'ko', markersize=0.5, label="filled")
    plt.plot(x2, y2, 'ro', markersize=0.5, label="original")
    plt.xlabel('Longitude', fontsize=12, labelpad=15)
    plt.ylabel('Latitude', fontsize=12, labelpad=30)


def spatial_smoothing(mc_ok, subcatalog, xmin, xmax, ymin, ymax):
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
    Smooth[:, 2] = mc_ok * np.ones((len(X) * len(Y)))
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


def ghost_element(catalog_sel, mmin, mc_ok, magcat, b, size, step, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat):
    # Simulates magnitude, time, lat and lon of ghost events

    import random
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import math
    import numpy as np

    serial_times = serial_time(catalog_sel)

    bin_m = 0.1

    # "hole" --> STAI gap
    hole_lower_lim, hole_upper_lim, mc_times_hole, mc_hole = mc_vs_time(magcat, mmin, mc_ok, serial_times, size, step)

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
                    random_m = round((- 1 / b * math.log10(u) + mc_ok - bin_m / 2) * 10) / 10

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

    # Simulate lat e lon of ghost events
    lon_ghost, lat_ghost = [], []
    for i in range(len(hole_lower_lim)):  # iterate over all the STAI gaps
        index = [idx for idx, val in enumerate(serial_times) if hole_lower_lim[i] <= val <= hole_upper_lim[i]]
        subcatalog_tmp = catalog_sel[index]
        mag_big_shock = max(subcatalog_tmp[:, MagnColumn])
        lat_big_shock = subcatalog_tmp[(subcatalog_tmp[:, MagnColumn] == mag_big_shock), LatColumn]
        lon_big_shock = subcatalog_tmp[(subcatalog_tmp[:, MagnColumn] == mag_big_shock), LonColumn]

        moment_big_shock = 10 ** (3 / 2 * (mag_big_shock + 10.7)) * 10 ** (-7)
        l_big_shock = 10 ** (-5.20 + 0.35 * math.log10(moment_big_shock))  # rupture length (from Mai and Beroza (2000))

        ymax, xmax = trainv(lat_big_shock, lon_big_shock, l_big_shock, l_big_shock)
        ymin, xmin = trainv(lat_big_shock, lon_big_shock, -l_big_shock, -l_big_shock)

        subcatalog = subcatalog_tmp[
                     (subcatalog_tmp[:, LatColumn] >= ymin) & (subcatalog_tmp[:, LatColumn] <= ymax) &
                     (subcatalog_tmp[:, LonColumn] >= xmin) & (subcatalog_tmp[:, LonColumn] <= xmax), :]

        lon_cell_grid, lat_cell_grid, smoothed_rate, sbin = spatial_smoothing(mc_ok, subcatalog, xmin, xmax, ymin,
                                                                              ymax)

        # Simulate lat and lon of ghost eqs
        idx = np.argsort(smoothed_rate)
        cumulative_sum = np.cumsum(smoothed_rate[idx])

        u = []
        for j in range(n_ghost_hole[i]):  # Generate n=n_ghost_hole[i] random numbers btw 0 and 1
            rand = random.random()
            u.append(rand)

        for j in range(len(u)):
            if u[j] < cumulative_sum[0]:
                rand_lon = lon_cell_grid[idx[0]] - random.random() * sbin
                rand_lat = lat_cell_grid[idx[0]] - random.random() * sbin
            else:
                for k in range(1, len(smoothed_rate)):
                    if (cumulative_sum[k - 1] <= u[j]) and (u[j] < cumulative_sum[k]):
                        rand_lon = lon_cell_grid[idx[k - 1]] + random.random() * sbin
                        rand_lat = lat_cell_grid[idx[k - 1]] + random.random() * sbin

            lon_ghost.append(round(rand_lon, 2))
            lat_ghost.append(round(rand_lat, 2))

    # plots
    # Plot original (m >= mc) + simulated eqs
    fig, ax = plt.subplots(2)
    formatter = mdates.DateFormatter("%m-%d")
    idx_magcatok = [idx for idx, val in enumerate(magcat) if val >= mc_ok]
    magcat_ok = [magcat[i] for i in idx_magcatok]
    serial_times_ok = [serial_times[i] for i in idx_magcatok]
    ax[0].xaxis_date()
    ax[0].xaxis.set_major_formatter(formatter)
    ax[0].scatter(serial_times_ok, magcat_ok, color='C1', s=1, label='Original data')
    ax[0].scatter(t_ghost, m_ghost, color='C4', s=1, label='Replenished data')
    ax[0].set_xlabel("Time")
    ax[0].set_ylabel("Magnitude")
    ax[0].legend(loc='upper left')
    ax[1].xaxis_date()
    ax[1].xaxis.set_major_formatter(formatter)
    ax[1].scatter(serial_times_ok, magcat_ok, color='C1', s=1, label='Original data')
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
    map_plot(lon_compl, lat_compl, lon_ghost, lat_ghost, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat)
    plt.savefig('fig/Spatial_map.pdf', format='pdf')

    return m_ghost, t_ghost, lon_ghost, lat_ghost


def timestamp_to_datetime(timestamp):
    import datetime

    days = timestamp % 1
    hours = days % 1 * 24
    minutes = hours % 1 * 60
    seconds = minutes % 1 * 60

    date_out = datetime.datetime.fromordinal(int(timestamp)) + datetime.timedelta(days=int(days)) + \
               datetime.timedelta(hours=int(hours)) + datetime.timedelta(minutes=int(minutes)) + \
               datetime.timedelta(seconds=seconds) - datetime.timedelta(days=366)

    return date_out


def replenished_catalog(catalog_sel, magcat, mmin, mc_ok, b, size, step, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat):
    # Returns the replenished catalog (original + ghost events)

    import numpy as np
    import csv
    import matplotlib.pyplot as plt

    serial_times = serial_time(catalog_sel)

    m_ghost, t_ghost, lon_ghost, lat_ghost = ghost_element(catalog_sel, mmin, mc_ok, magcat, b, size, step,
                                                           llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat)

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
    plt.scatter(x1, magcat_compl, 1)
    plt.savefig('fig/Magnitude_SeqNumbers_Original.pdf', format='pdf')

    fig, ax = plt.subplots()
    ax.set_xlabel("Sequential number")
    ax.set_ylabel("Magnitude")
    ax.set_title('Replenished catalog')
    plt.scatter(x2, sorted_new_mag, 1)
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
    return replenished_catalog





