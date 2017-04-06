
# coding: utf-8

# # Wrangling libraries

# In[1]:

# Import python modules
import os
import os.path
import sys
import csv
import ftplib
import bz2
import pandas
import wget

print 'Version1 4/4/17' 

# ## Define internal operations for downloading data
# In[1]:

# read in reference locations
def read_in_longlats(mappingfile):
    maptable=[]
    #pandas.read_table(mappingfile, sep=',')
    with open(mappingfile, 'rb') as csvfile:
        longlat = csv.reader(csvfile, delimiter=',')
        for row in longlat:
            maptable.append(row)
    csvfile.close()
    return(maptable)


# ### CIG (DHSVM)-oriented functions

# In[2]:

# index and extract longitude and latitude points
def compile_bc_Livneh2013_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')

    # compile a list of file urls for Livneh et al. 2013 (CIG)
    locations2013=[]
    for row in maptable:
        if maptable.index(row)!=0:
            basename='_'.join(['data',row[latitude], row[longitude]])
            url=['http://cses.washington.edu/rocinante/Livneh/bcLivneh_WWA_2013/forcs_dhsvm/',basename]
            locations2013.append(''.join(url))
    return(locations2013)

# index and extract longitude and latitude points
def compile_Livneh2013_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')

    # compile a list of lats and longs for Livneh 2013
    locations2013=[]
    for row in maptable:
        if maptable.index(row)!=0:
            basename='_'.join(['data',row[latitude], row[longitude]])
            url=['http://www.cses.washington.edu/rocinante/Livneh/Livneh_WWA_2013/forcs_dhsvm/',basename]
            locations2013.append(''.join(url))
    return(locations2013)


# ### VIC-oriented functions

# In[3]:

# compile file URLs
def compile_VICASCII_Livneh2016_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')
    
    # compile the VIC.ASCII data from Livneh 2016
    locations2016=[]
    for row in maptable:
        if maptable.index(row)!=0:
            loci='_'.join(['Fluxes_Livneh_NAmerExt_15Oct2014',row[latitude], row[longitude]])
            url=["ftp://192.12.137.7/pub/dcp/archive/OBS/livneh2014.1_16deg/VIC.ASCII/latitude.",row[latitude],'/',loci,'.bz2']
            locations2016.append(''.join(url))
    return(locations2016)

# compile file URLs
def compile_VICASCII_Livneh2013_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')
    
    # compile the VIC.ASCII data from Livneh 2013
    locations2013=[]
    for row in maptable:
        if maptable.index(row)!=0:
            loci='_'.join(['VIC_fluxes_Livneh_CONUSExt_v.1.2_2013',row[latitude], row[longitude]])
            url=["ftp://ftp.hydro.washington.edu/pub/blivneh/CONUS/Fluxes.asc.v.1.2.1915.2011.bz2/fluxes.125.120.37.49/",loci,".bz2"]
            #or
            #url=["ftp://ftp.hydro.washington.edu/pub/blivneh/CONUS/Fluxes.asc.v.1.2.1915.2011.bz2/fluxes.canada.columbia/",loci,".bz2"]
            
            locations2013.append(''.join(url))
    return(locations2013)

# compile file URLs
def compile_VICASCIIsubdaily_Livneh2013_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')
    
    # compile the VIC.ASCII data from Livneh 2013
    locations2013=[]
    for row in maptable:
        if maptable.index(row)!=0:
            loci='_'.join(['VIC_subdaily_fluxes_Livneh_CONUSExt_v.1.2_2013',row[latitude], row[longitude]])
            url=["ftp://ftp.hydro.washington.edu/pub/blivneh/CONUS/Derived.Subdaily.Outputs.asc.v.1.2.1915.2011.bz2/fluxes.125.120.37.49/",loci,".bz2"]
            locations2013.append(''.join(url))
    return(locations2013)


# ### Climate (Meteorological observations)-oriented functions

# In[ ]:

# index and extract longitude and latitude points for Livneh 2013
def compile_dailyMET_Livneh2013_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')
    
    # compile the daily MET data from Livneh 2013
    locations2013=[]
    for row in maptable:
        if maptable.index(row)!=0:
            loci='_'.join(['Meteorology_Livneh_CONUSExt_v.1.2_2013',row[latitude], row[longitude]])
            url=["ftp://ftp.hydro.washington.edu/pub/blivneh/CONUS/Meteorology.asc.v.1.2.1915.2011.bz2/data.125.120.37.49/",loci,".bz2"]
            locations2013.append(''.join(url))
    return(locations2013)

# index and extract longitude and latitude points for Livneh 2016
def compile_dailyMET_Livneh2016_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')

    # compile the daily data for Livneh 2016
    locations2016=[]
    for row in maptable:
        if maptable.index(row)!=0:
            loci='_'.join(['Meteorology_Livneh_NAmerExt_15Oct2014',row[latitude], row[longitude]])
            url=["ftp://192.12.137.7/pub/dcp/archive/OBS/livneh2014.1_16deg/ascii/daily/latitude.",row[latitude],"/",loci,".bz2"]
            locations2016.append(''.join(url))
    return(locations2016)


# ## Data file migration functions

# In[15]:

# check if the destination folder directory exists; if not, create it
def ensure_dir(f):
    if not os.path.exists(f):
        os.makedirs(f)
    os.chdir(f)

# Download the livneh 2013 files to the livneh2013 subdirectory
def wget_download(listofinterest, localfiledir):
           
    # check and download each location point, if it doesn't already exist in the download directory
    for fileurl in listofinterest:
        basename=os.path.basename(fileurl)
        filename=localfiledir + basename # file location on the local directory
        if os.path.isfile(filename):
            print('file already exists: ' + basename)
            continue
        wget.download(fileurl)
        wget.close()
        print(' downloaded: ' + basename)
        
# Download the livneh 2013 files to the livneh2013 subdirectory
def wget_download_one(fileurl):
    # check and download each location point, if it doesn't already exist in the download directory
    basename=os.path.basename(fileurl)
    filename=localfiledir + basename # file location on the local directory
    if os.path.isfile(filename):
        print('file already exists: ' + basename)
    else:
        wget.download(fileurl)
        print(' downloaded: ' + basename)
    
def wget_download_p(listofinterest):
    from multiprocessing import Pool
    pool = Pool(10) # submit 10 at once
    pool.map(wget_download_one, listofinterest)
    pool.close()
    pool.terminate()

# Download and decompress the livneh 2016 files
def ftp_download(listofinterest):
    for loci in listofinterest:
        
        # establish path info
        fileurl=loci.replace('ftp://','') # loci is already the url with the domain already appended
        ipaddress=fileurl.split('/',1)[0] # ip address
        path=os.path.dirname(fileurl.split('/',1)[1]) # folder path
        filename=os.path.basename(fileurl) # filename
        
        if os.path.isfile(filename):
            print('file already exists')
            continue
        
        # download the file from the ftp server
        ftp=ftplib.FTP(ipaddress)
        ftp.login()
        ftp.cwd(path)
        ftp.retrbinary("RETR " + filename ,open(filename, 'wb').write)
        ftp.close()
                
        # decompress the file
        decompbz2(filename)

# Download and decompress the livneh 2016 files
def ftp_download_one(loci):
    # establish path info
    fileurl=loci.replace('ftp://','') # loci is already the url with the domain already appended
    ipaddress=fileurl.split('/',1)[0] # ip address
    path=os.path.dirname(fileurl.split('/',1)[1]) # folder path
    filename=os.path.basename(fileurl) # filename
        
    if os.path.isfile(filename):
        print('file already exists')
    else:        
        # download the file from the ftp server
        ftp=ftplib.FTP(ipaddress)
        ftp.login()
        ftp.cwd(path)
        ftp.retrbinary("RETR " + filename ,open(filename, 'wb').write)
        ftp.quit()

        # decompress the file
        decompbz2(filename)

def ftp_download_p(listofinterest):
    from multiprocessing import Pool
    pool = Pool(10) # submit 10 at once
    pool.map(ftp_download_one, listofinterest)
    pool.close()
    pool.terminate()
    
# unzip the file
def decompbz2(filename):
    with open(filename.split(".bz2",1)[0], 'wb') as new_file, open(filename, 'rb') as zipfile:
        decompressor = bz2.BZ2Decompressor()
        for data in iter(lambda : zipfile.read(100 * 1024), b''):
            new_file.write(decompressor.decompress(data))
    os.remove(filename)
    zipfile.close()
    new_file.close()
    print(os.path.splitext(filename)[0] + ' unzipped')


# ## Wrapper scripts

# ### Get Daily Meteorological data from Livneh 2013
# formerly 'getClimateData_subdailyMET_Livneh2013.py'

# In[ ]:
# read in the longitude and latitude points from the reference mapping file
def getClimateData_DailyVIC_livneh2013(homedir, mappingfile):
    
    # generate table of lats and long coordinates
    maptable = read_in_longlats(mappingfile)
    
    # compile the longitude and latitude points
    dailyVIClocations2013 = compile_VICASCII_Livneh2013_locations(maptable)

    # check and generate VIC_ASCII Flux model livneh 2016 data directory
    filedir=homedir+'livneh2013/Daily_VIC_1915_2011/'
    ensure_dir(filedir)

    # Download the livneh 2016 VIC_ASCII Flux model data files
    ftp_download_p(dailyVIClocations2013)
    os.chdir(homedir)
    return(filedir)


# read in the longitude and latitude points from the reference mapping file
def getClimateData_DailyMET_livneh2013(homedir, mappingfile):
    
    # generate table of lats and long coordinates
    maptable = read_in_longlats(mappingfile)
    
    # compile the longitude and latitude points
    dailyMETlocations2013 = compile_dailyMET_Livneh2013_locations(maptable)

    # check and generate VIC_ASCII Flux model livneh 2016 data directory
    filedir=homedir+'livneh2013/Daily_MET_1915_2011/'
    ensure_dir(filedir)

    # Download the livneh 2016 VIC_ASCII Flux model data files
    ftp_download_p(dailyMETlocations2013)
    os.chdir(homedir)
    return(filedir)


# ### Get Daily Meteorological data from Livneh 2016
# formerly 'getClimateData_DailyMET_Livneh2016.py'

# In[ ]:

# read in the longitude and latitude points from the reference mapping file
def getClimateData_DailyMET_Livneh2016(homedir, mappingfile):
    
    # read in the longitude and latitude points from the reference mapping file
    maptable = read_in_longlats(mappingfile)
    
    # compile the longitude and latitude points
    dailyMETlocations2016 = compile_dailyMET_Livneh2016_locations(maptable)

    # check and generate baseline_corrected livneh2013 data directory
    filedir=homedir+'livneh2016/Daily_MET_1950_2013/'
    ensure_dir(filedir)

    # download the data
    ftp_download_p(dailyMETlocations2016)
    os.chdir(homedir)
    return(filedir)


# # Data Processing libraries

# In[ ]:

# import python modules
import pandas as pd
import datetime
import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
from collections import OrderedDict


# ## Define functions for reading-in downloaded files

# In[18]:

# create a list of files with their paths to be added to the HydroShare resource.
def compileContentfiles(directory):
    files = []
    filelist = os.listdir(directory) # the list of filenames
    for eachfile in filelist:
        files.append(directory + eachfile) # filepaths list
    return(files)


# Read in the mappingfile as a data frame
def mappingfileToDF(mappingfile):
    map_df= pandas.read_csv(mappingfile)

    # Select for station, latitude, longitude, and elevation columns
    if 'ELEV' in map_df.columns:
        map_df = map_df[['FID','LAT','LONG_','ELEV']]
    elif 'RASTERVALU' in map_df.columns:
        map_df = map_df[['FID','LAT','LONG_','RASTERVALU']]

    # Rename columns in climate locations dataframe
    map_df.columns = ['station','latitude','longitude','elevation']
    
    # compile summaries
    print map_df[0:4]
    print
    print 'Number of climate stations:', len(map_df)
    print 'Minimum elevation of climate stations:', np.min(map_df.elevation), 'm'
    print 'Mean elevation of climate stations:', int(np.mean(map_df.elevation)), 'm'
    print 'Maximum elevation of climate stations:', np.max(map_df.elevation), 'm'
    
    return(map_df, len(map_df))


# appends the folder path with the list of files within it
def filesWithPath(filedir):
    file_names=[]
    listOfFiles = os.listdir(filedir)
    for eachfile in listOfFiles:
        # logical rule to exclude hidden files
        if not eachfile.startswith('.'):
            file_names.extend([filedir+eachfile])
    return(file_names)

# define a function to read in the met files. The start and end date of the data must be input into the function since it is not 
# included in the raw data download.
def read_in_all_met_files(file_names, start_date, end_date):
    
    #initialize matrix and time sequence
    met_daily=[]
    met_daily_dates=pd.date_range(start_date, end_date, freq='D') # daily
    
    # Identify separator for files in the dataset
    sample = pd.read_table(file_names[0], header=None, sep='\t', nrows=1)
    if len(sample) != 1:
        separator = '\t'
    else:
        separator = '\s+'
    
    # import data for all climate stations
    for j in range(0, n_stations):
        met_daily.append(pd.read_table(file_names[j], header=None, sep=separator))
        met_daily[j].columns=['precip_mm','tmax_c', 'tmin_c', 'wind_m_s']
        met_daily[j].set_index(met_daily_dates, inplace=True) # set the index as the date

    # display first 10 lines of the first station's data frame
    print met_daily[0][0:10]
    return(met_daily)

# read in a daily streamflow data set
def read_daily_streamflow(file_name, drainage_area_m2, file_colnames=None, delimiter='\t', header='infer'):
    
    # if file_colnames are supplied, use header=None
    if file_colnames is not None:
        header=None
    
    # read in the data
    daily_data=pandas.read_table(file_name, delimiter=delimiter, header=header) 
    
    # set columns, if header=None
    if file_colnames is not None:
        daily_data.columns=file_colnames
    else:
        file_colnames=list(daily_data.columns)
        
    # calculate cfs to cms conversion, or vice versa
    if 'flow_cfs' in daily_data.columns:
        flow_cfs=daily_data['flow_cfs']
        flow_cms=flow_cfs/(3.28084**3)
        flow_mmday=flow_cms*1000*3600*24/drainage_area_m2
        
    elif 'flow_cms' in daily_data.columns:
        flow_cms=daily_data['flow_cms']
        flow_cfs=flow_cms*(3.28084**3)
        flow_mmday=flow_cms*1000*3600*24/drainage_area_m2
            
    # determine the datetime
    date_index=[file_colnames.index(each) for each in ['year','month','day']]
    row_dates=pandas.to_datetime(daily_data[date_index])
    
    # generate the daily_flow and set the datetime as row indices
    daily_flow=pandas.concat([flow_cfs, flow_cms, flow_mmday],axis=1)
    daily_flow.set_index(row_dates, inplace=True)
    daily_flow.columns=['flow_cfs', 'flow_cms', 'flow_mmday']
    return(daily_flow)

# read in a daily precipitation data set
def read_daily_precip(file_name, file_colnames=None, header='infer', delimiter='\s+'):
    
    # if file_colnames are supplied, use header=None
    if file_colnames is not None:
        header=None
    
    # read in the data
    daily_data=pandas.read_table(file_name, delimiter=delimiter, header=header) 
    
    # set columns, if header=None
    if file_colnames is not None:
        daily_data.columns=file_colnames
    else:
        file_colnames=list(daily_data.columns)
    
    # calculate cfs to cms conversion, or vice versa
    if 'precip_m' in daily_data.columns:
        precip_m=daily_data['precip_m']
        precip_mm=precip_m*1000
    
    # determine the datetime
    date_index=[file_colnames.index(each) for each in ['year','month','day']]
    row_dates=pandas.to_datetime(daily_data[date_index])
    
    # generate the daily_flow and set the datetime as row indices
    daily_precip=pandas.concat([precip_m, precip_mm],axis=1)
    daily_precip.set_index(row_dates, inplace=True)
    daily_precip.columns=['precip_m', 'precip_mm']
    return(daily_precip)

# read in a daily SNOTEL observation data set
def read_daily_snotel(file_name, file_colnames=None, usecols=None, delimiter=',', header='infer'):
    
    # if file_colnames are supplied, use header=None
    if file_colnames is not None:
        header=None
    
    # read in the data
    daily_data=pandas.read_table(file_name, usecols=usecols, header=header, delimiter=delimiter)
    
    # reset the colnames
    daily_data.columns=['Date', 'Tmax_C', 'Tmin_C', 'Tavg_C', 'Precip_mm']
    
    # transform the data
    daily_data['Tmax_C']=(daily_data['Tmax_C'] -32)/1.8
    daily_data['Tmin_C']=(daily_data['Tmin_C'] -32)/1.8
    daily_data['Tavg_C']=(daily_data['Tavg_C'] -32)/1.8
    daily_data['Precip_mm']=daily_data['Precip_mm'] *25.4
    
    # determine the datetime
    row_dates=pandas.to_datetime(daily_data.Date)
    
    # generate the daily_flow and set the datetime as row indices
    daily_snotel=daily_data[['Tmax_C', 'Tmin_C', 'Tavg_C', 'Precip_mm']]
    daily_snotel.set_index(row_dates, inplace=True)
    return(daily_snotel)


# ### Data Processing functions

# In[ ]:

# Determine dates (rows) and stations (columns). Number of stations 
# is the same for each dataset but number of dates may be different
def generateVarTables (listOfDates, listOfTables, n_stations):
    # NOTE: listOfTable must contain:
    # tmin_c
    # tmax_c
    # precip_mm
    # wind_m_s
    
    len_listOfDates=len(listOfDates) # number of dates
    station_list=range(0, n_stations) # list of stations number
    
    # Create arrays of for each variable of interest (Tmin, Tmax, Precip).
    # Rows are dates of analysis and columns are the station number
    temp_min_np=np.empty([len_listOfDates,n_stations])
    temp_max_np=np.empty([len_listOfDates,n_stations])
    precip_np=np.empty([len_listOfDates,n_stations])
    wind_np=np.empty([len_listOfDates,n_stations])
    
    # fill in each array with values from each station
    for i in station_list:
        temp_min_np[:,i]=listOfTables[i].tmin_c.values.astype(float)
        temp_max_np[:,i]=listOfTables[i].tmax_c.values.astype(float)
        precip_np[:,i]=listOfTables[i].precip_mm.values.astype(float)
        wind_np[:,i]=listOfTables[i].wind_m_s.values.astype(float)
        
    # generate each variable dataframe with rows as dates and columns as stations
    temp_min_df=pandas.DataFrame(temp_min_np, columns=station_list, index=listOfDates)    
    temp_max_df=pandas.DataFrame(temp_max_np, columns=station_list, index=listOfDates)    
    precip_df=pandas.DataFrame(precip_np, columns=station_list, index=listOfDates)    
    wind_df=pandas.DataFrame(wind_np, columns=station_list, index=listOfDates)
                                               
    # Create average temperature data frame as the average of Tmin and Tmax
    temp_avg_df=pandas.DataFrame((temp_min_np+temp_max_np)/2, columns=station_list, index=listOfDates)
    
    return(temp_min_df, temp_max_df, precip_df, wind_df, temp_avg_df)

# compare two date sets for the start and end of the overlapping dates
def overlappingDates(date_set1, date_set2):
    # find recent date
    if date_set1[0] > date_set2[0]:
        start_date = date_set1[0]
    else:
        start_date = date_set2[0]
    
    # find older date
    if date_set1[-1] < date_set2[-1]:
        end_date = date_set1[-1]
    else:
        end_date = date_set2[-1]
    return(start_date, end_date)

# Calculate means by 8 different methods
def multigroupMeans(VarTable, n_stations, start_date, end_date):
    Var_daily = VarTable.loc[start_date:end_date, range(0,n_stations)]
    
    # e.g., Mean monthly temperature at each station
    month_daily=Var_daily.groupby(Var_daily.index.month).mean() # average monthly minimum temperature at each station
    
    # e.g., Mean monthly temperature averaged for all stations in analysis
    meanmonth_daily=month_daily.mean(axis=1)
    
    # e.g., Mean monthly temperature for minimum and maximum elevation stations
    meanmonth_min_maxelev_daily=Var_daily.loc[:,analysis_elev_max_station].groupby(Var_daily.index.month).mean()
    meanmonth_min_minelev_daily=Var_daily.loc[:,analysis_elev_min_station].groupby(Var_daily.index.month).mean()
    
    # e.g., Mean annual temperature
    year_daily=Var_daily.groupby(Var_daily.index.year).mean()
    
    # e.g., mean annual temperature each year for all stations
    meanyear_daily=year_daily.mean(axis=1)
    
    # e.g., mean annual min temperature for all years, for all stations
    meanallyear_daily=np.nanmean(meanyear_daily)
    
    # e.g., anomoly per year compared to average
    anom_year_daily=meanyear_daily-meanallyear_daily
    
    return(month_daily, 
           meanmonth_daily, 
           meanmonth_min_maxelev_daily, 
           meanmonth_min_minelev_daily, 
           year_daily, 
           meanyear_daily, 
           meanallyear_daily,
           anom_year_daily)

def specialTavgMeans(VarTable):
    Var_daily = VarTable.loc[start_date:end_date, range(0,n_stations)]
    
    # Average temperature for each month at each station
    permonth_daily=Var_daily.groupby(pandas.TimeGrouper("M")).mean()
    
    # Average temperature each month averaged at all stations
    meanpermonth_daily=permonth_daily.mean(axis=1)
    
    # Average monthly temperature for all stations
    meanallpermonth_daily=meanpermonth_daily.mean(axis=0)
    
    # anomoly per year compared to average
    anom_month_daily=(meanpermonth_daily-meanallpermonth_daily)/1000
    
    return(permonth_daily,
          meanpermonth_daily,
          meanallpermonth_daily,
          anom_month_daily)

