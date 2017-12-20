"""
Module for general utilities

 - INCLUDE DROPOUT TEST
 - INCLUDE THIRD MOMENT TEST
 - CHANGE NOTATION IN QCONTROL'S SUMMARY
"""
from __future__ import absolute_import, print_function, division
from . import tests

def qc_replace(files, fileconfig,
             read_files_kw={'parse_dates':False, 'clean_dates':False, 'only_named_cols':False, 'return_units':False},
             begin_date=None, end_date=None, # date check
             file_lines=None, # lines test
             nans_test=True, # nans test
             lower_limits={}, upper_limits={}, # boundaries test
             spikes_test=True, # spikes test
             visualize_spikes=False,
             chunk_size=1200,
             spikes_vis_col='u',
             spikes_detrend = {'how':'linear'},
             spikes_func = lambda x: (abs(x - x.median()) > (7./0.6745)*x.mad()), 
             max_consec_spikes=10,
             max_replacement_count=180, # replacement count test
             replace_with='interpolation',
             passverbose=False, failverbose=True, # general
             passshow=False, failshow=False,
             passshow_vars=None,
             outdir='0_replaced',
             summary_file='control_replacement.csv',
             replaced_report='rreport.csv',
             full_report=None):

    """
    This function applies various quality control checks/tests to a set of datafiles and re-writes
    the successful files in another directory. This specific function focuses on point-analysis
    that can be fixed by replacements (removing spikes by interpolation, etc). A run fails this
    only if (1) it is beyond the accepted dates, (2) has a different number of lines than
    usual and (3) the number of points replaced is greater than max_replacement_count.

    A list of applied tests is found below in order of application. The only test available 
    by default is the spikes test. All others depend on their respective keywords.

    Tests are based on Vickers and Mahrt, Quality control and fux sampling problems for
    tower and aircraft data.

    Consistency checks
    ------------------

    - :date check:
        Files outside a date_range are left out.
        keywords: end_date, begin_date
    - :lines test:
        Checks each file to see if they have a certain number of lines.
        Files with a different number of lines fail this test. 
        keyworkds: file_lines

    Quality tests (in this order)
    -----------------------------
    - :NaN's filter:
        Checks for any NaN values. NaNs are replaced with interpolation or linear trend.
        Activate it by passing nans_test=True.
        - keywords: nans_test

    - :boundaries filter:
        Checks for values in any column lower than a pre-determined lower limit or higher
        than a upper limit. If found, these points are replacted (interpolated or otherwise).
        - keywords: lower_limits, upper_limits

    - :spikes filter:
        Search for spikes according to user definition. Spikes are replaced (interpolated or otherwise).
        - keywords: spikes_test, spikes_func, visualize_spikes, spikes_vis_col, 
                    max_consec_spikes and chunk_size keywords.
    - :replacement count test:
        Checks the total amount of points that were replaced (including NaN, boundaries and spikes test)
        against the max_replacement_count keyword. Fails if any columns has more replacements than that.
        - keywords: max_replacement_count


    Parameters
    ----------
    files: list
        list of filepaths
    fileconfig: pymicra.fileConfig object or str
        datalogger configuration object used for all files in the list of files or path to a dlc file.
    read_files_kw: dict
        keywords to pass to pymicra.timeSeries. Default is
        {'parse_dates':False} because parsing dates at every file is slow, so this
        makes the whole process faster.  However, {'parse_dates':True,
        'clean_dates':False} is recommended if time is not a problem because the window
        and chunk_size keywords may be used as, for example '2min', instead of 1200,
        which is the equivalent number of points.
    file_lines: int
        number of line a "good" file must have. Fails if the run has any other number of lines.
    begin_date: str
        dates before this automatically fail.
    end_date: str
        dates after this automatically fail.
    nans_test: bool
        whether or not to apply the nans test
    lower_limits: dict
        keys must be names of variables and values must be lower absolute limits for the values of each var.
    upper_limits: dict
        keys must be names of variables and values must be upper absolute limits for the values of each var.
    spikes_test: bool
        whether or not to check for spikes.
    spikes_detrend: dict
        keywords to pass to pymicra.detrend when detrending for spikes. Is it's empty, no detrending is done.
    visualize_spikes: bool
        whether or not to plot the spikes identification and interpolation (useful for calibration of spikes_func). Only
        one column is visualized at each time. This is set with the spikes_vis_col keyword.
    spikes_vis_col: str
        column to use to visualize spikes.
    spikes_func: function
        function used to look for spikes. Can be defined used numpy/pandas notation for methods with lambda functions.
        Default is: lambda x: (abs(x - x.mean()) > abs(x.std()*4.))
    replace_with: str
        method to use when replacing the spikes. Options are 'interpolation' and 'trend'.
    max_consec_spikes: int
        Limit of consecutive spike points to be interpolated. If the number of consecutive "spikes" is more than
        this, then we take all those points as not actually being spikes and no replacement is done. So if
        max_consec_spikes=0, no spike replacements is ever done.
    chunk_size: str
        string representing time length of chunks used in the spikes check.
        Default is "2Min".  Putting None will not separate in chunks. It's recommended
        to use rolling functions in this case (might be slow).
    max_replacement_count: int
        Maximum number of replaced point a variable can have in a run. If the
        replaced number of points is larger than this then the run fails and is
        discarded. Generally this should be about 1% of the file_lines.
    passverbose: bool
        whether or not to show details on the successful runs.
    failverbose: bool
        whether or not to show details on the failed runs.
    passshow: bool
        whether of not to plot the successful runs on screen.
    passshow_vars: list
        list of columns to plot if run is successfull.
    failshow: bool
        whether of not to plot the failed runs on screen.
    outdir: str
        name of directory in which to write the successful runs. Directory must already exist.
    summary_file: str
        path of file to be created with the summary of the runs. Will be overwriten if already exists.

    Returns
    -------
    ext_summary: pandas.DataFrame
        dict with the extended summary, which has the path of the files that got "stuck" in each test along with the successful ones
    """
    from . import timeSeries
    from os.path import basename, join
    from dateutil.parser import parse
    import pandas as pd
    import numpy as np
    from . import algs
    from .io import write_as_fconfig

    if passshow:
        import matplotlib.pyplot as plt

    if begin_date: begin_date=parse(begin_date)
    if end_date: end_date=parse(end_date)

    total_name='total'
    lines_name='failed lines test'
    replacement_name = 'failed replacement test'
    successful_name = 'passed all tests'
    replaced_nans_name = 'Runs with replaced nans'
    replaced_bound_name = 'Runs with replaced bound'
    replaced_spikes_name = 'Runs with replaced spikes'

    order = [total_name, lines_name, replacement_name,
                successful_name, replaced_nans_name, replaced_bound_name, replaced_spikes_name]

    #--------------
    # If the path to the dlc is provided, we read it as a dataloggerConfig object
    if isinstance(fileconfig, str):
        from . import fileConfig
        fileconfig = fileConfig(fileconfig)
    #--------------

    #-------------------------------------
    # Identifying columns that are not part of the datetime
    variables_list=fileconfig.variables
    if type(variables_list) == dict:
        usedvars=[ v for v in variables_list.values() if r'%' not in v ]
    else:
        raise TypeError('Check variables of the fileConfig object.')
    #-------------------------------------

    #--------------
    # We first create the dataframe to hols our limit values and the discarded dict, which
    # is what we use to produce our summary
    tables = pd.DataFrame(columns=usedvars)
    replaced = pd.DataFrame(columns=usedvars)
    control = pd.DataFrame({total_name: [], successful_name: []})
    #--------------

    #--------------
    # We update discarded and tables based on the tests we will perform.
    # If a test is not marked to be perform, it will not be on this list.
    if file_lines:
        control[ lines_name ] = []
    if nans_test:
        control[ replaced_nans_name ] = []

    if upper_limits:
        tables = tables.append( pd.DataFrame(upper_limits, index=['upper_limits']) )
        control[ replaced_bound_name ] = []
    if lower_limits:
        tables = tables.append( pd.DataFrame(lower_limits, index=['lower_limits']) )
        control[ replaced_bound_name ] = []

    if spikes_test:
        control[ replaced_spikes_name ] = []

    if max_replacement_count:
        control[ replacement_name ] = []

    tables = tables.fillna(value=np.nan)
    #--------------

    #-------------------------------------
    # BEGINNING OF MAIN PROGRAM
    #-------------------------------------

    for idx, filepath in enumerate(files):
        filename=basename(filepath)
        print(filename)

        #---------------
        # DATE CHECK
        if begin_date or end_date:
            cdate = algs.name2date(filename, fileconfig)
            if begin_date:
                if cdate<begin_date:
                    print('Skipped because of begin_date.\n')
                    continue
            if end_date: 
                if cdate>end_date:
                    print('Skipped because of end_date.\n')
                    continue
        #----------------

        #----------------
        # If the test passes the date check then we include it in the total amount
        control.loc[idx, total_name ] = filename
        #----------------
    
        #-------------------------------
        # LINE NUMBERS TEST
        if file_lines:
            valid = tests.check_numlines(filepath, numlines=file_lines, failverbose=failverbose)

            result, failed = algs.testValid(valid, testname=lines_name, passverbose=passverbose, filepath=filepath, failverbose=failverbose)
            if result == False:
                #discarded[ lines_name ].append(filename)
                control.loc[idx, lines_name ] = filename
                continue
        #-------------------------------
    
        #-------------------------------
        # OPENNING OF THE FILE HAPPENS HERE
        fin=timeSeries(filepath, fileconfig, **read_files_kw)
        #-------------------------------

        #-------------------------------
        # We save the full input for writting it later and exclude unnused variables
        fullfin=fin.copy()
        fin=fin[usedvars].copy()
        replaced.loc[filename] = 0
        #-------------------------------

        #-----------------
        # CHECK NANS
        if nans_test:
            valid, nans_replaced = tests.check_nans(fin, replace_with=replace_with)

#            result, failed = algs.testValid(valid, testname=nan_name, passverbose=passverbose, filepath=filepath, failverbose=failverbose)
#            control = algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=nan_name, filename=filename, failshow=failshow)

            #--------------
            # Add nans that were replaced to the full replaced list
            replaced.loc[ filename ] += nans_replaced
            control.loc[ idx, replaced_nans_name ] = nans_replaced.sum()
            #--------------

            if result==False: continue
        #-----------------

        #-------------------------------
        # BEGINNING OF LOWER AND UPPER VALUES CHECK (BOUNDARIES TEST)
        if lower_limits or upper_limits:
            fin, valid, limits_replaced = tests.check_limits(fin, tables, replace_with=replace_with)

#            result, failed = algs.testValid(valid, testname=bound_name, passverbose=passverbose, filepath=filepath, failverbose=failverbose)
#            control = algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=bound_name, filename=filename, failshow=failshow)

            #--------------
            # Add high/low values that were replaced to the full replaced list
            replaced.loc[ filename ] += limits_replaced
            control.loc[ idx, replaced_bound_name ] = limits_replaced.sum()
            #--------------

            if result==False: continue
        #----------------------------------
 
        #-----------------
        # BEGINNING OF SPIKES CHECK
        if spikes_test:
            fin, valid, spikes_replaced = tests.check_spikes(fin, detrend=spikes_detrend,
                            visualize=visualize_spikes, vis_col=spikes_vis_col, chunk_size=chunk_size, replace_with=replace_with,
                            cut_func=spikes_func, max_consec_spikes=max_consec_spikes)

#            result, failed = algs.testValid(valid, testname=spikes_name, passverbose=passverbose, filepath=filepath, failverbose=failverbose)
#            control = algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=spikes_name, filename=filename, failshow=failshow)

            #--------------
            # Add spikes that were replaced to the full replaced list
            replaced.loc[ filename ] += spikes_replaced
            control.loc[ idx, replaced_spikes_name ] = spikes_replaced.sum()
            #--------------

            if result==False: continue
        #-----------------
   
        #-----------------
        # REPLACEMENT COUNT TEST
        if max_replacement_count:
            valid = tests.check_replaced(replaced.iloc[-1], max_count=max_replacement_count)

            result, failed = algs.testValid(valid, testname=replacement_name, passverbose=passverbose, filepath=filepath, failverbose=failverbose)
            control = algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=replacement_name, filename=filename, failshow=failshow)

            if result==False: continue
        #-----------------

 
        #-----------------
        # END OF TESTS
        print('Passed all tests')
        if passshow:
            if passshow_vars:
                fin.loc[:, passshow_vars]
            else:
                fin.plot()
            plt.show()
        control.loc[idx, successful_name ] = filename
        #-----------------

        #-----------------
        # FINALLY, we write the result in the output directory in the same format
        if outdir:
            print('Re-writing',filepath)
            fullfin[usedvars] = fin[usedvars]       # This is because some spikes were removed during the process
            write_as_fconfig(fullfin, join(outdir, basename(filepath)), fileconfig)
        #-----------------
        print()
    
    if replaced_report:
        replaced.to_csv(replaced_report, na_rep='NaN')
    if full_report:
        control.to_csv(full_report, na_rep='NaN')
    #-----------------
    # We create the summary dataframe which is the count of all non-NaN values
    summary = pd.DataFrame({'control': control.count()})
    #-----------------
    
    #-----------------
    # But this doesnt work for the count of runs that had replacement, because zero is also a non-NaN value
    if nans_test:
        summary.loc[ replaced_nans_name ] = control[ replaced_nans_name ].astype(bool).sum()
    if spikes_test:
        summary.loc[ replaced_spikes_name ] = control[ replaced_spikes_name ].astype(bool).sum()
    if lower_limits or upper_limits:
        summary.loc[ replaced_bound_name ] = control[ replaced_bound_name ].astype(bool).sum()
    #-------------

    #-------------
    # We calculate the percentage and print to the user
    summary['percent']=100.*summary['control']/summary.loc['total', 'control']
    summary = summary.loc[ order ].dropna()
    summary[ 'control' ] = summary[ 'control' ].astype(int)
    print(summary)
    #-------------

    summary.to_csv(summary_file, na_rep='NaN')
    return control
 



def qc_discard(files, fileconfig,
             read_files_kw={'parse_dates':False, 'clean_dates':False, 'only_named_cols':False, 'return_units':False},
             std_limits={}, # min std test
             std_detrend=dict(how='linear'),
             dif_limits={}, # max. dif test
             maxdif_detrend=dict(),
             maxdif_trend=dict(how='movingmedian', window=600),
             chunk_size=1200,
             passverbose=False, failverbose=True,
             failshow=False, passshow=False,
             passshow_vars=None,
             outdir='1_filtered',
             summary_file='filter_summary.csv',
             full_report=None):

    """
    Function that applies various tests quality control to a set of datafiles and re-writes
    the successful files in another directory. A list of currently-applied tests is found 
    below in order of application. When some variable or set of points fails a test the
    whole file is discarded.

    Tests are based on Vickers and Mahrt, Quality control and fux sampling problems for
    tower and aircraft data.


    Quality tests (in this order)
    -----------------------------

    - :standard deviation (STD) check:
        runs with a standard deviation lower than a pre-determined value are left out.
        - keywords: std_limits, std_detrend

    - :maximum difference (stationarity) test:
        runs whose trend have a maximum difference greater than a certain value are left out.
        This excludes non-stationary runs. Activate it by passing a dif_limits keyword.
        - keywords: dif_limits, maxdif_detrend, maxdif_trend

    Parameters
    ----------
    files: list
        list of filepaths
    fileconfig: pymicra.fileConfig object or str
        datalogger configuration object used for all files in the list of files or path to a dlc file.
    read_files_kw: dict
        keywords to pass to pymicra.timeSeries. Default is
        {'parse_dates':False} because parsing dates at every file is slow, so this
        makes the whole process faster. However, {'parse_dates':True,
        'clean_dates':False} is recommended if time is not a problem because the window
        and chunk_size keywords may be used as, for example '2min', instead of 1200,
        which is the equivalent number of points.
    dif_limits: dict
        keys must be names of variables and values must be upper limits for the maximum difference
        of values that the linear trend of the run must have.
    maxdif_detrend: dict
        keywords to pass to pymicra.detrend when detrending for max difference test. If it's empty, no detrending is made.
    maxdif_trend: dict
        Keywords to pass to pymicra.detrend when trending for max difference
        test (passed to pymicra.data.trend).  If empty, no trending is used.  This is
        used in the max difference test, since the difference is taken between the max
        and min values of the trend, not of the raw timeSeries.
    chunk_size: str
        string representing time length of chunks used in the standard deviation check. Default is "2Min".
        Putting None will not separate in chunks. It's recommended to use rolling functions in this case (might be slow).
    passverbose: bool
        whether or not to show details on the successful runs.
    failverbose: bool
        whether or not to show details on the failed runs.
    passshow: bool
        whether of not to plot the successful runs on screen.
    passshow_vars: list
        list of columns to plot if run is successfull.
    failshow: bool
        whether of not to plot the failed runs on screen.
    outdir: str
        name of directory in which to write the successful runs. Directory must already exist.
    summary_file: str
        path of file to be created with the summary of the runs. Will be overwriten if already exists.

    Returns
    -------
    ext_summary: pandas.DataFrame
        dict with the extended summary, which has the path of the files that got "stuck" in each test along with the successful ones
    """
    from . import timeSeries
    from os.path import basename, join
    from dateutil.parser import parse
    import pandas as pd
    import numpy as np
    from . import algs
    from .io import write_as_fconfig

    if passshow:
        import matplotlib.pyplot as plt

    total_name='total'
    STD_name='failed STD test'
    maxdif_name='failed maxdif test'
    successful_name = 'passed all tests'

    order = [total_name, STD_name, maxdif_name, successful_name]

    #--------------
    # If the path to the dlc is provided, we read it as a dataloggerConfig object
    if isinstance(fileconfig, str):
        from . import fileConfig
        fileconfig = fileConfig(fileconfig)
    #--------------

    #-------------------------------------
    # Identifying columns that are not part of the datetime
    variables_list=fileconfig.variables
    if type(variables_list) == dict:
        usedvars=[ v for v in variables_list.values() if r'%' not in v ]
    else:
        raise TypeError('Check variables of the fileConfig object.')
    #-------------------------------------

    #--------------
    # We first create the dataframe to hols our limit values and the discarded dict, which
    # is what we use to produce our summary
    tables = pd.DataFrame(columns=usedvars)
    control = pd.DataFrame({total_name: [], successful_name: []})
    #--------------

    #--------------
    # We update discarded and tables based on the tests we will perform.
    # If a test is not marked to be perform, it will not be on this list.
    if std_limits:
        tables = tables.append( pd.DataFrame(std_limits, index=['std_limits']) )
        control[ STD_name ] = []

    if dif_limits:
        tables = tables.append( pd.DataFrame(dif_limits, index=['dif_limits']) )
        control[ maxdif_name ] = []

    tables = tables.fillna(value=np.nan)
    #--------------

    #-------------------------------------
    # BEGINNING OF MAIN PROGRAM
    #-------------------------------------

    for idx, filepath in enumerate(files):
        filename=basename(filepath)
        print(filename)

        #-------------------------------
        # OPENNING OF THE FILE HAPPENS HERE
        fin=timeSeries(filepath, fileconfig, **read_files_kw)
        #-------------------------------

        #-------------------------------
        # We save the full input for writting it later and exclude unnused variables
        fullfin=fin.copy()
        fin=fin[usedvars].copy()
        #-------------------------------

        #----------------
        # If the test passes the date check then we include it in the total amount
        control.loc[idx, total_name ] = filename
        #----------------
    
        #----------------------------------
        # STANDARD DEVIATION TEST
        if std_limits:
            valid = tests.check_std(fin, tables, detrend=std_detrend, chunk_size=chunk_size, failverbose=failverbose)

            result, failed = algs.testValid(valid, testname=STD_name, passverbose=passverbose, filepath=filepath, failverbose=failverbose)
            control = algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=STD_name, filename=filename, failshow=failshow)

            if result==False: continue
        #----------------------------------
    
        #-------------------------------
        # STATIONARITY TEST
        if dif_limits:
            valid = tests.check_stationarity(fin, tables, detrend=maxdif_detrend,
                                        trend=maxdif_trend)

            result, failed = algs.testValid(valid, testname=maxdif_name, passverbose=passverbose, filepath=filepath)
            control=algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=maxdif_name, filename=filename, failshow=failshow)

            if result==False: continue
        #-------------------------------
    
        #-----------------
        # END OF TESTS
        print('Passed all tests')
        if passshow:
            if passshow_vars:
                fin.loc[:, passshow_vars]
            else:
                fin.plot()
            plt.show()
        control.loc[idx, successful_name ] = filename
        #-----------------

        #-----------------
        # FINALLY, we write the result in the output directory in the same format
        if outdir:
            print('Re-writing',filepath)
            fullfin[usedvars] = fin[usedvars]       # This is because some spikes were removed during the process
            write_as_fconfig(fullfin, join(outdir, basename(filepath)), fileconfig)
        #-----------------
        print()
    
    if full_report:
        control.to_csv(full_report, na_rep='NaN')
    #-----------------
    # We create the summary dataframe which is the count of all non-NaN values
    summary = pd.DataFrame({'control': control.count()})
    #-----------------
    
    #-------------
    # We calculate the percentage and print to the user
    summary['percent']=100.*summary['control']/summary.loc['total', 'control']
    summary = summary.loc[ order ].dropna()
    summary[ 'control' ] = summary[ 'control' ].astype(int)
    print(summary)
    #-------------

    summary.to_csv(summary_file, na_rep='NaN')
    return control
 





def separateFiles(files, dlconfig, outformat='out_%Y-%m-%d_%H:%M.csv', outdir='',
                verbose=False, firstflag='.first', lastflag='.last', save_ram=False,
                frequency='30min', quoting=0, use_edges=False):
    """
    Separates files into (default) 30-minute smaller files. Useful for output files such
    as the ones by Campbell Sci, that can have days of data in one single file.

    Parameters
    -----------
    files: list
        list of file paths to be separated
    dlconfig: pymicra datalogger configuration file
        to tell how the dates are displayed inside the file
    outformat: str
        the format of the file names to output
    outdir: str
        the path to the directory in which to output the files
    verbose: bool
        whether to print to the screen
    firstflag: str
        flag to put after the name of the file for the first file to be created
    lastflag: str
        flag to put after the name of the fle for the last file to be created
    save_ram: bool
        if you have an amount of files that are to big for pandas to load on your ram this should be set to true
    frequency:
        the frequency in which to separate
    quoting: int
        for pandas (see read_csv documentation)
    use edges: bool
        use this carefully. This concatenates the last few lines of a file to the first few lines
        of the next file in case they don't finish on a nice round time with respect to the frequency

    Returns
    -------
    None
    """
    from os import path
    import io
    import pandas as pd
    from . import algs

    #-----------------------
    # First considerations
    outpath=path.join(outdir, outformat)
    firstpath= path.join(outdir, outformat+firstflag)
    lastpath = path.join(outdir, outformat+ lastflag)
    parser = lambda x: algs.line2date(x, dlconfig)
    badfiles = []
    #-----------------------

    #------------
    # if there's no need to save RAM, then I will use pandas, which is faster
    if save_ram == False:
    #------------
        for f in files:
            df=io.timeSeries(f, dlconfig, parse_dates=True, 
                parse_dates_kw={'clean' : False}, 
                read_data_kw={'quoting' : quoting})
            chunks, index = algs.splitData(df, frequency, return_index=True)
            for idx, chunk in zip(index,chunks):
                out = idx.strftime(outpath)
                if verbose: print('Writting chunk to ',out)
                chunk.to_csv(out, index=False, header=False)
        return
    
    #------------
    # If RAM consumption is an issue, then I won't use pandas. This saves RAM but takes a lot longer.
    else:
        header = dlconfig.header_lines
        if header==None: header=0
    #------------
        for fin in files:
            print('Now opening',fin)
            #------------
            # Creates a sequence of equaly-spaced dates based on the frequency and nicely rounded-up
            ft, lt = algs.first_last(fin)
            ft, lt = map(parser, [ft, lt])
            labeldates = pd.Series(index=pd.date_range(start=ft, end=lt, freq='min')).resample(frequency).index
            nfiles=len(labeldates)
            if nfiles == 0:
                badfiles.append(fin)
                continue
            #------------

            with open(fin, 'rt') as fin:
                for i in range(header):
                    fou.write(fin.readline())
                lines=[]
                pos=fin.tell()
                for i, line in enumerate(fin):
                    cdate=parser(line)
                    if (cdate>=labeldates[0]) and (cdate<labeldates[1]):
                        lines.append(line)
                    elif (cdate>=labeldates[1]):
                        if verbose:
                            print(labeldates[0])
                        #------------
                        # labels first and last runs separately, so they can be identified and put together later
                        if len(labeldates) == nfiles:
                            fou = open((labeldates[0]).strftime(firstpath), 'wt')
                        else:
                            fou = open((labeldates[0]).strftime(outpath), 'wt')
                        #------------
                        fou.writelines(lines)
                        fou.close()
                        labeldates=labeldates.drop(labeldates[0])
                        #---------
                        # this is in case the data has a big jump of more than the value of the frequency
                        while True:
                            if (len(labeldates)>1) and (cdate>=labeldates[1]):
                                labeldates=labeldates.drop(labeldates[0])
                            elif (cdate>=labeldates[0]) or (len(labeldates)==1):
                                break
                            else:
                                raise ValueError('SOMETHING WRONG\nCHECK ALGORITHM')
                        #---------

                        #---------
                        # Starts over the lines
                        lines = [line]
                        #---------
                #---------
                # This prints the last lines of the file
                        if len(labeldates)==1: break
                for line in fin:
                    lines.append(line)
                fou=open((labeldates[0]).strftime(lastpath), 'wt')
                fou.writelines(lines)
                fou.close()
                #---------

        #----------------
        # This feature avoids losing the ends of large files by concatenating the end of a file
        # with the beginning of the next one (only when either files do not end on a rounded time)
        if use_edges:
        #----------------
            from glob import glob
            ledges = sorted(glob(path.join(outdir,'*'+lastflag )))
            fedges = sorted(glob(path.join(outdir,'*'+firstflag)))
            for last in ledges:
                root = last[:-len(lastflag)]
                try:
                    idx = fedges.index(root+firstflag)
                except ValueError:
                    continue
                last_lines = open(last,'rt').readlines()
                first = fedges.pop(idx)
                first_lines = open(first,'rt').readlines()
                if verbose:
                    print('Concatenating ', last, ' and ', first)
                with open(root, 'wt') as fou:
                    fou.writelines(last_lines)
                    fou.writelines(first_lines)

        print('List of bad files:', badfiles)
        if verbose:
            print('Done!')
        return


def correctDrift(drifted, correct_drifted_vars=None, correct=None,
                get_fit=True, write_fit=True, fit_file='correctDrift_linfit.params',
                apply_fit=True, show_plot=False, return_plot=False, units={}, return_index=False):
    """
    Parameters
    -----------
    correct: pandas.DataFrame
        dataset with the correct averages
    drifted: pandas.DataFrame
        dataset with the averages that need to be corrected
    correct_drifted_vars: dict
        dictionary where every key is a var in the right dataset and 
        its value is its correspondent in the drifted dataset
    get_fit: bool
        whether ot not to fit a linear relation between both datasets. Generally slow. Should only be done once
    write_fit: bool
        if get_fit == True, whether or not to write the linear fit to a file (recommended)
    fit_file: string
        where to write the linear fit (if one is written) or from where to read the linear fit (if no fit is written)
    apply_fit: bool
        whether of not to apply the lineat fit and correct the data (at least get_fit and fit_file must be true)
    show_plot: bool
        whether or not to show drifted vs correct plot, to see if it's a good fit
    units: dict
        if given, it creates a {file_file}.units file, to tell write down in which units data has to be in
        order to be correctly corrected
    return_index: bool
        whether to return the indexes of the used points for the calculation. Serves to check the regression

    Returns
    -------
    outdf: pandas.DataFrame
        drifted dataset corrected with right dataset
    """
    from matplotlib import pyplot as plt
    import pandas as pd
    import numpy as np

    if correct_drifted_vars:
        rwvars = correct_drifted_vars
    else:
        if len(correct.columns)==1:
            rwvars = { cor : dft for cor, dft in zip(correct.columns, drifted.columns) }
        else:
            raise NameError('If correct is not provided or has more than one column, you should provide correct_drifted_vars.')

    cors=[]
    #----------------
    # This option is activated if we provide a correct dataset from which to withdraw the correction parameters
    if get_fit:
        for slw, fst in rwvars.items():
            slow=correct[slw]
            fast=drifted[fst]
            #----------------
            # Check to see if the frequency in both datasets are the same. Otherwise we are comparing different things
            try:
                if pd.infer_freq(correct.index) == pd.infer_freq(drifted.index):
                    slow, fast = map(np.array, [slow, fast] )
                else:
                    print('Frequencies must be the same, however, inferred frequencies appear to be different. Plese check.')
            except TypeError:
                print('Cannot determine if frequencies are the same. We will continue but you should check')
                slow, fast = map(np.array, [slow, fast] )
            #----------------
    
            #----------------
            # Does the 1D fitting filtering for NaN values (very important apparently)
            idx = np.isfinite(slow) & np.isfinite(fast)
            coefs, residuals, rank, singular_vals, rcond = np.polyfit(fast[idx], slow[idx], 1, full=True)
            #----------------
            if show_plot:
                plt.title('{} vs {}'.format(fst, slw))
                plt.plot(fast[idx], slow[idx], marker='o', linestyle='')
                plt.plot(fast[idx], np.poly1d(coefs)(fast[idx]), '-', linewidth=2)
                plt.xlabel(fst)
                plt.ylabel(slw)
                plt.grid(True)
                fig = plt.gcf()
                plt.show()
                if return_plot:
                    return fig
    
            correc=pd.DataFrame(columns=[ '{}_{}'.format(slw, fst) ], index=['angular', 'linear'], data=coefs).transpose()
            cors.append(correc)
        cors = pd.concat(cors, join='outer')
        print(cors)
    #----------------

        #----------------
        # Writes the fit parameters in a file to be used later
        if write_fit:
            cors.index.name='correct_drifted'
            cors.to_csv(fit_file, index=True)
            if units:
                with open(fit_file+'.units', 'wt') as fou:
                    for key, item in units.items():
                        fou.write('{"%s":"%s"}\n' % (key,item))
        #----------------

    #----------------
    # If you do not want to correct from an existing correct dataset. A file with the parameters must be read
    else:
        cors=pd.read_csv(fit_file, index_col=0, header=0)
    #----------------

    #------------
    # Applies the fit column by column
    if apply_fit:
        corrected=drifted.copy()
        for slw, fst in rwvars.items():
            coefs = np.array(cors.loc['{}_{}'.format(slw,fst), ['angular','linear']])
            corrected[ fst ] = np.poly1d(coefs)(drifted[ fst ])
    else:
        corrected=drifted.copy()
    #------------

    #----------------
    # The returning of the index idx is done mainly for checking purposes
    if return_index:
        return corrected, idx
    else:
        return corrected
    #----------------



def qcontrol(files, fileconfig,
             read_files_kw={'parse_dates':False, 'clean_dates':False, 'only_named_cols':False, 'return_units':False},
             begin_date=None, end_date=None, # date check
             file_lines=None, # lines test
             nans_test=True, # nans test
             accepted_nans_percent=1.,
             lower_limits={}, upper_limits={}, # boundaries test
             accepted_bound_percent=1.,
             spikes_test=True, # spikes test
             visualize_spikes=False,
             spikes_vis_col='u',
             spikes_detrend = {'how':'linear'},
             spikes_func = lambda x: (abs(x - x.mean()) > 4.*x.std()), 
             max_consec_spikes=3,
             accepted_spikes_percent=1.,
             max_replacement_count=180, # replacement count test
             std_limits={}, # min std test
             std_detrend=True, std_detrend_kw={'how':'movingmean', 'window':900},
             dif_limits={}, # max. dif test
             maxdif_detrend=True,
             maxdif_detrend_kw={'how':'movingmean', 'window':900},
             maxdif_trend=True, maxdif_trend_kw={'how':'movingmedian', 'window':600},
             RAT = False, # reverse arrangement test
             RAT_vars = None,
             RAT_detrend=True, RAT_detrend_kw={'how':'linear'},
             RAT_points = 50, RAT_significance = 0.05,
             chunk_size=1200,
             replace_with='interpolation',
             trueverbose=False, falseverbose=True,
             falseshow=False, trueshow=False,
             trueshow_vars=None,
             outdir='quality_controlled',
             summary_file='qcontrol_summary.csv',
             replaced_report=None,
             full_report=None):

    """
    Function that applies various tests quality control to a set of datafiles and re-writes
    the successful files in another directory. A list of currently-applied tests is found 
    below in order of application. The only test available by default is the spikes test.
    All others depend on their respective keywords.


    Consistency tests
    -----------------

    - :date check:
        files outside a date_range are left out.
        keywords: end_date, begin_date
    - :lines test:
        checks each file to see if they have a certain number of lines.
        Files with a different number of lines fail this test. 
        keyworkds: file_lines

    Quality tests (in this order)
    -----------------------------
    - :NaN's test:
        checks for any NaN values. NaNs are replaced with interpolation or linear trend. If the percentage
        of NaNs is greater than accepted_nans_percent, run is discarded. Activate it by passing nans_test=True.
        - keywords: accepted_nans_percent, nans_test

    - :boundaries test:
        runs with values in any column lower than a pre-determined lower limit or higher
        than a upper limits are left out.
        - keywords: lower_limits, upper_limits

    - :spikes test:
        replace for spikes and replace them according to some keywords.
        runs with more than a certain percetage of spikes are left out. 
        - keywords: spikes_test, spikes_func, visualize_spikes, spikes_vis_col, 
                    max_consec_spikes, accepted_spikes_percent and chunk_size keywords.
    - :replacement count test:
        checks the total amount of points that were replaced (including NaN, boundaries and spikes test)
        against the max_replacement_count keyword. Fails if any columns has more replacements than that.
        - keywords: max_replacement_count

    - :standard deviation (STD) check:
        runs with a standard deviation lower than a pre-determined value (generally close to the
        sensor precision) are left out.
        - keywords: std_limits, std_detrend, std_detrend_kw

    - :maximum difference (stationarity) test:
        runs whose trend have a maximum difference greater than a certain value are left out.
        This excludes non-stationary runs. Activate it by passing a dif_limits keyword.
        - keywords: dif_limits, maxdif_detrend, maxdif_detrend_kw, maxdif_trend, maxdif_trend_kw

    - :reverse arrangement test (RAT):
        runs that fail the reverse arrangement test for any variable are left out.
        - keywords: RAT, RAT_vars, RAT_detrend, RAT_detrend_kw, RAT_points, RAT_significance

    Parameters
    ----------
    files: list
        list of filepaths
    fileconfig: pymicra.fileConfig object or str
        datalogger configuration object used for all files in the list of files or path to a dlc file.
    read_files_kw: dict
        keywords to pass to pymicra.timeSeries. Default is {'parse_dates':False} because parsing dates
        at every file is slow, so this makes the whole process faster.
        However, {'parse_dates':True, 'clean_dates':False} is recommended if time is not a problem because
        the window and chunk_size keywords may be used as, for example '2min', instead of 1200, which is the
        equivalent number of points.
    file_lines: int
        number of line a "good" file must have. Fails if the run has any other number of lines.
    begin_date: str
        dates before this automatically fail.
    end_date: str
        dates after this automatically fail.
    nans_test: bool
        whether or not to apply the nans test
    accepted_nans_percent: float
        percentage of runs that is acceptable
    std_limits: dict
        keys must be names of variables and values must be upper limits for the standard deviation.
    std_detrend: bool
        whether or not to work with the fluctations of the data on the spikes and standard deviation test.
    std_detrend_kw:
        keywords to be passed to pymicra.detrend specifically to be used on the STD test.
    lower_limits: dict
        keys must be names of variables and values must be lower absolute limits for the values of each var.
    upper_limits: dict
        keys must be names of variables and values must be upper absolute limits for the values of each var.
    dif_limits: dict
        keys must be names of variables and values must be upper limits for the maximum difference
        of values that the linear trend of the run must have.
    maxdif_detrend: bool
        whether to detrend data before checking for differences.
    maxdif_detrend_kw: dict
        keywords to pass to pymicra.detrend when detrending for max difference test.
    maxdif_trend: bool
        whether to check for differences using the trend, instead of raw points (which can be the fluctuations
        or the original absolute values of data, depending if maxdif_detrend==True or False).
    maxdif_trend_kw: dict
        keywords to pass to pymicra.detrend when trending for max difference test.
        dictionary of keywords to pass to pymicra.data.trend. This is used in the max difference test, since
        the difference is taken between the max and min values of the trend, not of the series.
        Default = {'how':'linear'}.
    spikes_test: bool
        whether or not to check for spikes.
    spikes_detrend: dict
        keywords to pass to pymicra.detrend when detrending for spikes.
    visualize_spikes: bool
        whether or not to plot the spikes identification and interpolation (useful for calibration of spikes_func). Only
        one column is visualized at each time. This is set with the spikes_vis_col keyword.
    spikes_vis_col: str
        column to use to visualize spikes.
    spikes_func: function
        function used to look for spikes. Can be defined used numpy/pandas notation for methods with lambda functions.
        Default is: lambda x: (abs(x - x.mean()) > abs(x.std()*4.))
    replace_with: str
        method to use when replacing the spikes. Options are 'interpolation' and 'trend'.
    max_consec_spikes: int
        limit of consecutive spike points to be interpolated. After this spikes are left as they are in the output.
    accepted_percent: float
        limit percentage of spike points in the data. If spike points represent a higher percentage
        than the run fails the spikes check.
    chunk_size: str
        string representing time length of chunks used in the spikes and standard deviation check. Default is "2Min".
        Putting None will not separate in chunks. It's recommended to use rolling functions in this case (might be slow).
    RAT: bool
        whether or not to perform the reverse arrangement test on data.
    RAT_vars: list
        list containing the name of variables to go through the reverse arrangement test. If None, all variables are tested.
    RAT_points: int
        number of final points to apply the RAT. If 50, the run will be averaged to a 50-points run.
    RAT_significance:
        significance level to apply the RAT.
    RAT_detrend_kw:
        keywords to be passed to pymicra.detrend specifically to be used on the RA test. {"how":"linear"}
        is strongly recommended for this case.
    trueverbose: bool
        whether or not to show details on the successful runs.
    falseverbose: bool
        whether or not to show details on the failed runs.
    trueshow: bool
        whether of not to plot the successful runs on screen.
    trueshow_vars: list
        list of columns to plot if run is successfull.
    falseshow: bool
        whether of not to plot the failed runs on screen.
    outdir: str
        name of directory in which to write the successful runs. Directory must already exist.
    summary_file: str
        path of file to be created with the summary of the runs. Will be overwriten if already exists.

    Returns
    -------
    ext_summary: pandas.DataFrame
        dict with the extended summary, which has the path of the files that got "stuck" in each test along with the successful ones
    """
    from . import timeSeries
    from os.path import basename, join
    from dateutil.parser import parse
    import pandas as pd
    import numpy as np
    from . import algs
    from .io import write_as_fconfig

    if trueshow:
        import matplotlib.pyplot as plt

    if begin_date: begin_date=parse(begin_date)
    if end_date: end_date=parse(end_date)

    total_name='total'
    lines_name='failed lines test'
    nan_name='failed NaNs test'
    bound_name='failed boundaries test'
    spikes_name='failed spikes test'
    replacement_name = 'failed max. replacement test'
    STD_name='failed STD test'
    STD_stat_name='failed STD stationarity test'
    maxdif_name='failed maxdif test'
    RAT_name = 'failed RAT test'
    successful_name = 'passed all tests'
    replaced_nans_name = 'Runs with replaced nans'
    replaced_bound_name = 'Runs with replaced bound'
    replaced_spikes_name = 'Runs with replaced spikes'

    order = [total_name, lines_name, nan_name, bound_name, spikes_name, replacement_name, STD_name, STD_stat_name,
                maxdif_name, RAT_name, successful_name, replaced_nans_name, replaced_bound_name, replaced_spikes_name]

    #--------------
    # If the path to the dlc is provided, we read it as a dataloggerConfig object
    if isinstance(fileconfig, str):
        from . import fileConfig
        fileconfig = fileConfig(fileconfig)
    #--------------

    #-------------------------------------
    # Identifying columns that are not part of the datetime
    variables_list=fileconfig.variables
    if type(variables_list) == dict:
        usedvars=[ v for v in variables_list.values() if r'%' not in v ]
    else:
        raise TypeError('Check variables of the fileConfig object.')
    #-------------------------------------

    #--------------
    # We first create the dataframe to hols our limit values and the discarded dict, which
    # is what we use to produce our summary
    tables = pd.DataFrame(columns=usedvars)
    replaced = pd.DataFrame(columns=usedvars)
    control = pd.DataFrame({total_name: [], successful_name: []})
    #--------------

    #--------------
    # We update discarded and tables based on the tests we will perform.
    # If a test is not marked to be perform, it will not be on this list.
    if file_lines:
        control[ lines_name ] = []
    if nans_test:
        control[ nan_name ] = []
        control[ replaced_nans_name ] = []

    if upper_limits:
        tables = tables.append( pd.DataFrame(upper_limits, index=['upper_limits']) )
        control[ bound_name ] = []
        control[ replaced_bound_name ] = []
    if lower_limits:
        tables = tables.append( pd.DataFrame(lower_limits, index=['lower_limits']) )
        control[ bound_name ] = []
        control[ replaced_bound_name ] = []

    if spikes_test:
        control[ spikes_name ] = []
        control[ replaced_spikes_name ] = []

    if max_replacement_count:
        control[ replacement_name ] = []

    if std_limits:
        tables = tables.append( pd.DataFrame(std_limits, index=['std_limits']) )
        control[ STD_name ] = []

    if RAT:
        control[ RAT_name ] = []

    if dif_limits:
        tables = tables.append( pd.DataFrame(dif_limits, index=['dif_limits']) )
        control[ maxdif_name ] = []

    tables = tables.fillna(value=np.nan)
    #--------------

    #-------------------------------------
    # BEGINNING OF MAIN PROGRAM
    #-------------------------------------

    for idx, filepath in enumerate(files):
        filename=basename(filepath)
        print(filename)

        #---------------
        # DATE CHECK
        if begin_date or end_date:
            cdate = algs.name2date(filename, fileconfig)
            if begin_date:
                if cdate<begin_date:
                    print('Skipped because of begin_date.\n')
                    continue
            if end_date: 
                if cdate>end_date:
                    print('Skipped because of end_date.\n')
                    continue
        #----------------

        #----------------
        # If the test passes the date check then we include it in the total amount
        control.loc[idx, total_name ] = filename
        #----------------
    
        #-------------------------------
        # LINE NUMBERS TEST
        if file_lines:
            valid = tests.check_numlines(filepath, numlines=file_lines, falseverbose=falseverbose)

            result, failed = algs.testValid(valid, testname=lines_name, trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
            if result == False:
                #discarded[ lines_name ].append(filename)
                control.loc[idx, lines_name ] = filename
                continue
        #-------------------------------
    
        #-------------------------------
        # OPENNING OF THE FILE HAPPENS HERE
        # TRY-EXCEPT IS A SAFETY NET BECAUSE OF THE POOR DECODING (2015-06-21 00:00 appears as 2015-06-20 24:00)
        try:
            fin=timeSeries(filepath, fileconfig, **read_files_kw)
        except ValueError as e:
            if str(e)=='unconverted data remains: 0' and cdate.hour==23:
                continue
            else:
                raise ValueError(e)
        #-------------------------------

        #-------------------------------
        # We save the full input for writting it later and exclude unnused variables
        fullfin=fin.copy()
        fin=fin[usedvars].copy()
        replaced.loc[filename] = 0
        #-------------------------------

        #-----------------
        # CHECK NANS
        if nans_test:
            valid, nans_replaced = tests.check_nans(fin, max_percent=accepted_nans_percent, replace_with=replace_with)

            result, failed = algs.testValid(valid, testname=nan_name, trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
            control = algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=nan_name, filename=filename, falseshow=falseshow)

            #--------------
            # Add nans that were replaced to the full replaced list
            replaced.loc[ filename ] += nans_replaced
            control.loc[ idx, replaced_nans_name ] = nans_replaced.sum()
            #--------------

            if result==False: continue
        #-----------------

        #-------------------------------
        # BEGINNING OF LOWER AND UPPER VALUES CHECK (BOUNDARIES TEST)
        if lower_limits or upper_limits:
            fin, valid, limits_replaced = tests.check_limits(fin, tables, max_percent=accepted_bound_percent, replace_with=replace_with)

            result, failed = algs.testValid(valid, testname=bound_name, trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
            control = algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=bound_name, filename=filename, falseshow=falseshow)

            #--------------
            # Add high/low values that were replaced to the full replaced list
            replaced.loc[ filename ] += limits_replaced
            control.loc[ idx, replaced_bound_name ] = limits_replaced.sum()
            #--------------

            if result==False: continue
        #----------------------------------
 
        #-----------------
        # BEGINNING OF SPIKES CHECK
        if spikes_test:
            fin, valid, spikes_replaced = tests.check_spikes(fin, detrend=spikes_detrend,
                            visualize=visualize_spikes, vis_col=spikes_vis_col, chunk_size=chunk_size, replace_with=replace_with,
                            cut_func=spikes_func, max_consec_spikes=max_consec_spikes, max_percent=accepted_spikes_percent)

            result, failed = algs.testValid(valid, testname=spikes_name, trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
            control = algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=spikes_name, filename=filename, falseshow=falseshow)

            #--------------
            # Add spikes that were replaced to the full replaced list
            replaced.loc[ filename ] += spikes_replaced
            control.loc[ idx, replaced_spikes_name ] = spikes_replaced.sum()
            #--------------

            if result==False: continue
        #-----------------
   
        #-----------------
        # REPLACEMENT COUNT TEST
        if max_replacement_count:
            valid = tests.check_replaced(replaced.iloc[-1], max_count=max_replacement_count)

            result, failed = algs.testValid(valid, testname=replacement_name, trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
            control = algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=replacement_name, filename=filename, falseshow=falseshow)

            if result==False: continue
        #-----------------

        #----------------------------------
        # STANDARD DEVIATION TEST
        if std_limits:
            valid = tests.check_std(fin, tables, detrend=std_detrend, chunk_size=chunk_size, falseverbose=falseverbose)

            result, failed = algs.testValid(valid, testname=STD_name, trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
            control = algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=STD_name, filename=filename, falseshow=falseshow)

            if result==False: continue
        #----------------------------------
    
        #-------------------------------
        # STATIONARITY TEST
        if dif_limits:
            valid = tests.check_stationarity(fin, tables, detrend=maxdif_detrend, trend=maxdif_trend)

            result, failed = algs.testValid(valid, testname=maxdif_name, trueverbose=trueverbose, filepath=filepath)
            control=algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=maxdif_name, filename=filename, falseshow=falseshow)

            if result==False: continue
        #-------------------------------
    
        #------------------------------
        # REVERSE ARRANGEMENT TEST
        if RAT:
            valid = tests.check_RA(fin, detrend=RAT_detrend, detrend_kw=RAT_detrend_kw,
                                    RAT_vars=None, RAT_points=RAT_points, RAT_significance=RAT_significance)

            result, failed = algs.testValid(valid, testname=RAT_name, trueverbose=trueverbose, filepath=filepath)
            control = algs.applyResult(result, failed, fin, control=control, index_n=idx, testname=RAT_name , filename=filename, falseshow=falseshow)

            if result==False: continue
        #-------------------------------
 
        #-----------------
        # END OF TESTS
        print('Passed all tests')
        if trueshow:
            if trueshow_vars:
                fin.loc[:, trueshow_vars]
            else:
                fin.plot()
            plt.show()
        #discarded[ successful_name ].append(filename)
        control.loc[idx, successful_name ] = filename
        #-----------------

        #-----------------
        # FINALLY, we write the result in the output directory in the same format
        if outdir:
            print('Re-writing',filepath)
            fullfin[usedvars] = fin[usedvars]       # This is because some spikes were removed during the process
            write_as_fconfig(fullfin, join(outdir, basename(filepath)), fileconfig)
#            fullfin.to_csv(join(outdir, basename(filepath)),
#                       header=fileconfig.header, index=False, quoting=3, na_rep='NaN')
        #-----------------
        print()
    
    if replaced_report:
        replaced.to_csv(replaced_report, na_rep='NaN')
    if full_report:
        control.to_csv(full_report, na_rep='NaN')
    #-----------------
    # We create the summary dataframe which is the count of all non-NaN values
    summary = pd.DataFrame({'control': control.count()})
    #-----------------
    
    #-----------------
    # But this doesnt work for the count of runs that had replacement, because zero is also a non-NaN value
    if nans_test:
        summary.loc[ replaced_nans_name ] = control[ replaced_nans_name ].astype(bool).sum()
    if spikes_test:
        summary.loc[ replaced_spikes_name ] = control[ replaced_spikes_name ].astype(bool).sum()
    if lower_limits or upper_limits:
        summary.loc[ replaced_bound_name ] = control[ replaced_bound_name ].astype(bool).sum()
    #-------------

    #-------------
    # We calculate the percentage and print to the user
    summary['percent']=100.*summary['control']/summary.loc['total', 'control']
    summary = summary.loc[ order ].dropna()
    summary[ 'control' ] = summary[ 'control' ].astype(int)
    print(summary)
    #-------------

    summary.to_csv(summary_file, na_rep='NaN')
    return control
 





def _printUnit(string, mode='L', trim=True, greek=True):
    """
    Returns string formatted for LaTeX or other uses.

    Parameters
    ----------
    string: string
        string (unambiguous) that represents the unit
    mode: string
        {'L' : 'latex', 'P' : 'pretty', 'H' : 'html' }
    trim: bool
        if true it trims the output
    greek: bool
        yet to be implemented
    """
    from . import ureg as ur
    ur.default_format=mode
    u=string
    u='{:~L}'.format(u)
    if trim:
        u=u[3:].strip()
    return u


