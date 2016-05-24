#!/usr/bin/python
"""
Author: Tomas Chor
Date: 2015-08-07
-------------------------

TODO LIST
-INCLUDE DECODIFICAION OF DATA?

"""


def check_spikes(dfs, visualize=False, vis_col=1, interp_limit=3,
                 cut_func=lambda x: (abs(x - x.mean()) > 4.*abs(x.std())) ):
    '''
    Applies spikes-check according to Vickers and Mart (1997)

    Parameters:
    -----------
    dfs: list, tuple
        sequence of pandas.DataFrame objects
    visualize: bool
        whether of not to visualize the interpolation ocurring
    vis_col: str, int or list
        the column(s) to visualize when seeing the interpolation (only effective if visualize==True)
    interp_limit: int
        limit of consecutive spikes to interpolate
    cut_func: function
        function used to define spikes
    '''
    import pandas as pd
    import algs
    import matplotlib.pyplot as plt

    original = pd.concat(dfs)
    valid_cols=pd.Series(0, index=dfs[0].columns)

    for i in range(len(dfs)):
        #-------------------------------
        # We make a copy of the original df just in case
        chunk=dfs[i].copy()
        #-------------------------------

        #-------------------------------
        # This substitutes the spikes to NaNs so it can be interpolated later
        if len(chunk)>interp_limit:
            chunk=algs.limitedSubs(chunk, max_interp=interp_limit, func=cut_func)
        valid_cols = valid_cols + chunk.count()
        #-------------------------------

        #-------------------------------
        #Interpolation takes place here
        chunk=chunk.interpolate(method='time', axis=0)
        #-------------------------------

        #-------------------------------
        # We change the chunk in the original list of dfs to concatenate later
        dfs[i]=chunk.copy()
        #-------------------------------

    fou=pd.concat(dfs)
    fou=fou.interpolate(method='time', axis=0)
    fou_points=len(fou.index)
    valid_cols=valid_cols/fou_points

    #-------------------------------
    # Visualize what you're doing to see if it's correct
    if visualize:
        print 'Plotting de-spikes...'
        original[vis_col].plot(style='g-', label='original')
        #aux[ cut_func(aux) ].plot(style='ro-', label='spikes')
        fou[vis_col].plot(style='b-', label='final')
        plt.title('Column: {}'.format(vis_col))
        plt.legend()
        plt.show()
        plt.close()
    #-------------------------------

    return fou, valid_cols


def qcontrol(files, datalogger_config,
             detrend=False, accepted_percent=1.,
             file_lines=None, bdate=None, edate=None,
             std_limits={}, dif_limits={}, low_limits={}, upp_limits={},
             spikes_check=True, visualize_spikes = False, spikes_vis_col = 'u',
             spikes_func = lambda x: (abs(x - x.mean()) > 3.*abs(x.std())), 
             interp_limit=3, window_size=900, chunk_size='2Min',
             rev_arrang_test = False, RATvars = None,
             RAT_points = 50, RAT_significance = 0.01,
             trueverbose=False, falseverbose=True, falseshow=False, 
             trueshow=False, trueshow_vars=None,
             outdir='quality_controlled',
             summary_file='qcontrol_summary.csv'):

    """
    Function that applies various tests quality control to a set of datafiles and re-writes
    the successful files in another directory. A list of currently-applied tests is found 
    below in order of application. The only test available by default is the spikes test.
    All others depend on their respective keywords.


    Trivial tests:
    --------------
    date check:
        files outside a date_range are left out (edate and bdate keywords)
    lines test:
        files with a number of lines that is different from the correct number are out.

    Non-trivial tests:
    ------------------
    lowest value test:
        runs with values in any column lower then a pre-determined value are left out.
        Activate it by passing a low_limits keyword.
    highest value test:
        runs with values in any column higher then a pre-determined value are left out.
        Activate it by passing a upp_limits keyword.
    spikes test:
        runs with more than a certain percetage of spikes are left out. 
        Activate it by passing a spikes_check keyword. Adjust the test with the spikes_func
        visualize_spikes, spikes_vis_col, interp_limit, accepted_percent and chunk_size keywords.
    standard deviation check:
        runs with a standard deviation lower than a pre-determined value (generally close to the
        sensor precision) are left out.
        Activate it by passing a std_limits keyword.
    reverse arrangement test:
        runs that fail the reverse arrangement test for any variable are left out.
        Activate it by passing a rev_arrang_test keyword.
    maximum difference test:
        runs whose trend have a maximum different greater than a certain value are left out.
        This excludes non-stationary runs.
        Activate it by passing a dif_limits keyword.

    Parameters:
    -----------
    files: list
        list of filepaths
    datalogger_config: pymicra.datalogerConf object or str
        datalogger configuration object used for all files in the list of files or path to a dlc file.
    detrend: bool
        whether or not to work with the fluctations of the data in the tests where absolute values
        don't matter (spikes, standard deviation and reverse arrangement tests).
    file_lines: int
        number of line a "good" file must have. Fails if the run has any other number of lines.
    bdate: str
        dates before this automatically fail.
    edate: str
        dates after this automatically fail.
    std_limits: dict
        keys must be names of variables and values must be upper limits for the standard deviation.
    dif_limits: dict
        keys must be names of variables and values must be upper limits for the maximum difference
        of values that the linear trend of the run must have.
    low_limits: dict
        keys must be names of variables and values must be lower absolute limits for the values of each var.
    upp_limits: dict
        keys must be names of variables and values must be upper absolute limits for the values of each var.
    spikes_check: bool
        whether or not to check for spikes.
    visualize_spikes: bool
        whether or not to plot the spikes identification and interpolation (useful for calibration of spikes_func). Only
        one column is visualized at each time. This is set with the spikes_vis_col keyword.
    spikes_vis_col: str
        column to use to visualize spikes.
    spikes_func: function
        function used to look for spikes. Can be defined used numpy/pandas notation for methods with lambda functions.
        Default is: lambda x: (abs(x - x.mean()) > abs(x.std()*4.))
    interp_limit: int
        limit of consecutive spike points to be interpolated. After this spikes are left as they are in the output.
    accepted_percent: float
        limit percentage of spike points in the data. If spike points represent a higher percentage
        than the run fails the spikes check.
    chunk_size: str
        string representing time length of chunks used in the spikes and standard deviation check. Default is "2Min".
        Putting None will not separate in chunks. It's recommended to use rolling functions in this case.
    window_size: int
        window size for rolling mean used in the standard deviation test.
    rev_arrang_test: bool
        whether or not to perform the reverse arrangement test on data.
    RATvars: list
        list containing the name of variables to go through the reverse arrangement test. If None, all variables are tested.
    RAT_points: int
        number of final points to apply the RAT. If 50, the run will be averaged to a 50-points run.
    RAT_significance:
        significance level to apply the RAT.
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

    Returns:
    --------
    ext_summary: pandas.DataFrame
        dict with the extended summary, which has the path of the files that got "stuck" in each test along with the successful ones
    """
    from io import timeSeries
    from os.path import basename, join
    from dateutil.parser import parse
    import data
    import pandas as pd
    import numpy as np
    from algs import auxiliar as algs

    if trueshow:
        import matplotlib.pyplot as plt

    if bdate: bdate=parse(bdate)
    if edate: edate=parse(edate)

    #--------------
    # If the path to the dlc is provided, we read it as a dataloggerConfig object
    if isinstance(datalogger_config, str):
        from io import read_dlc
        datalogger_config = read_dlc(datalogger_config)
    #--------------

    #--------------
    # We first create the dataframe to hols our limit values and the numbers dict, which
    # is what we use to produce our summary
    tables=pd.DataFrame()
    numbers = {'total': [], 'successful': []}
    #--------------

    #--------------
    # We update numbers and tables based on the tests we will perform.
    # If a test is not marked to be perform, it will not be on this list.
    if bdate or edate:
        numbers['dates'] = []
    if spikes_check:
        numbers['spikes'] = []
    if file_lines:
        numbers['lines'] = []
    if std_limits:
        tables = tables.append( pd.DataFrame(std_limits, index=['std_limits']) )
        numbers['STD'] = []
    if low_limits:
        tables = tables.append( pd.DataFrame(low_limits, index=['low_limits']) )
        numbers['lowest value' ] = []
    if upp_limits:
        tables = tables.append( pd.DataFrame(upp_limits, index=['upp_limits']) )
        numbers['highest value' ] = []
    if dif_limits:
        tables = tables.append( pd.DataFrame(dif_limits, index=['dif_limits']) )
        numbers['difference'] = []
    if rev_arrang_test:
        numbers['RAT'] = []
    tables = tables.fillna(value=np.nan)
    #--------------

    #-------------------------------------
    # Identifying columns that are not part of the datetime
    variables_list=datalogger_config.varNames
    if type(variables_list) == dict:
        usedvars=[ v for v in variables_list.values() if r'%' not in v ]
    elif type(variables_list) == list:
        usedvars=[ v for v in variables_list[1:] if r'%' not in v ]
    else:
        raise TypeError('Check varNames of the dataloggerConfig object.')
    #-------------------------------------

    #-------------------------------------
    # BEGINNING OF MAIN PROGRAM
    #-------------------------------------

    for filepath in files:
        print
        filename=basename(filepath)
        print filename
        numbers['total'].append(filename)
        #---------------
        # BEGINNING OF DATE CHECK
        if bdate or edate:
            cdate = algs.name2date(filename, datalogger_config)
            if bdate:
                if cdate<bdate:
                    if falseverbose: print filename, 'failed dates check'
                    numbers['dates'].append(filename)
                    continue
            if edate: 
                if cdate>edate:
                    if falseverbose: print filename, 'failed dates check'
                    numbers['dates'].append(filename)
                    continue
            #-----------------
            # If we get here, then test is successful
            if trueverbose: print filename, 'passed dates check'
            #-----------------
        #----------------
    
        #-------------------------------
        # BEGINNING LINE NUMBERS TEST
        if file_lines:
            result=algs.check_numlines(filepath, numlines=file_lines)
            if result == False:
                if falseverbose: print filepath,'failed lines test'
                numbers['lines'].append(filename)
                continue
            else:
                if trueverbose: print filename, 'passed lines check'
        #-------------------------------
    
        #-------------------------------
        # OPENNING OF THE FILE HAPPENS HERE
        # TRY-EXCEPT IS A SAFETY NET BECAUSE OF THE POOR DECODING (2015-06-21 00:00 appears as 2015-06-20 24:00)
        try:
            fin=timeSeries(filepath, datalogger_config, parse_dates_kw={'correct_fracs':True, 'clean':False})
        except ValueError, e:
            if str(e)=='unconverted data remains: 0' and cdate.hour==23:
                continue
            else:
                raise ValueError, e
        #-------------------------------

        #-------------------------------
        # We save the full input for writting it later
        fullfin=fin.copy()
        fin=fin[usedvars]      # Exclude unnused variables
        #-------------------------------

        #-------------------------------
        # BEGINNING OF LOWEST VALUE CHECK
        if low_limits:
            valid= ~(fin < tables.loc['low_limits']).any(axis=0)

            result, failed=algs.testValid(valid, testname='lowest value', trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
            numbers=algs.applyResult(result, failed, fin, control=numbers, testname='lowest value', filename=filename, falseshow=falseshow)
            if result==False: continue
        #-------------------------------
    
        #-------------------------------
        # BEGINNING OF HIGHEST VALUE CHECK
        if upp_limits:
            valid= ~(fin > tables.loc['upp_limits']).any(axis=0)

            result, failed=algs.testValid(valid, testname='highest value', trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
            numbers=algs.applyResult(result, failed, fin, control=numbers, testname='highest value', filename=filename, falseshow=falseshow)
            if result==False: continue
        #-------------------------------

        #----------------------------------
        # BEGINNING OF STANDARD DEVIATION CHECK
        if std_limits:
            df=fin.copy()
            df=df-pd.rolling_mean(df,window=window_size, center=True)
            stds_list=df.resample(chunk_size, np.std).dropna()
            valid= ~(stds_list<tables.loc['std_limits']).any(axis=0)

            result, failed = algs.testValid(valid, testname='STD', trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
            #--------------
            # This specifies which part of the run failed on this test
            if falseverbose and (result==False):
                print (not result)*'The failed variables and times:\n{0}'.format(~(stds_list<tables.loc['std_limits']))
            #--------------
            numbers=algs.applyResult(result, failed, fin, control=numbers, testname='STD', filename=filename, falseshow=falseshow)
            if result==False: continue
        #----------------------------------
    
        #---------------------------------
        # BEGINNING OF REVERSE ARRANGEMENT TEST
        if rev_arrang_test:
            if RATvars:
                valid_chunks= fin[RATvars].apply(data.reverse_arrangement, axis=0, points_number=RAT_points, alpha=RAT_significance)
            elif RATvars==None:
                valid_chunks= fin.apply(data.reverse_arrangement, axis=0, points_number=RAT_points, alpha=RAT_significance)
            else:
                valid_chunks= fin.any(axis=0)

            result, failed = algs.testValid(valid_chunks, testname='reverse arrangement', trueverbose=trueverbose, filepath=filepath)
            numbers=algs.applyResult(result, failed, fin, control=numbers, testname='RAT', filename=filename, falseshow=falseshow)
            if result==False: continue
        #---------------------------------
    
        #--------------------------------
        # BEGINNING OF MAXIMUM DIFFERENCE METHOD
        if dif_limits:
            trend=data.trend(fin, mode='linear')
            maxdif=trend.iloc[0]-trend.iloc[-1]
            maxdif=maxdif.abs()
            chunks_valid= tables.loc['dif_limits'] - maxdif
            chunks_valid= ~(chunks_valid < 0)

            result, failed = algs.testValid(chunks_valid, testname='difference', trueverbose=trueverbose, filepath=filepath)
            numbers=algs.applyResult(result, failed, fin, control=numbers, testname='difference', filename=filename, falseshow=falseshow)
            if result==False: continue
        #--------------------------------
 
        #-------------------------------
        # BEGINNING OF SPIKES CHECK
        if spikes_check:
            chunks = algs.splitData(fin, chunk_size)
            fin,valid_cols=check_spikes(chunks, visualize=visualize_spikes, vis_col=spikes_vis_col, cut_func=spikes_func, interp_limit=interp_limit)

            valid= valid_cols >= (1.-(accepted_percent/100.))

            result, failed=algs.testValid(valid, testname='spikes', trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
            numbers=algs.applyResult(result, failed, fin, control=numbers, testname='spikes', filename=filename, falseshow=falseshow)
            if result==False: continue
        #-------------------------------
   
        #-----------------
        # END OF TESTS
        print 'Successful run!'
        if trueshow:
            if trueshow_vars:
                fin.loc[:, trueshow_vars]
            else:
                fin.plot()
            plt.show()
        numbers['successful'].append(filename)
        #-----------------

        #-----------------
        # FINALLY, we write the result in the output directory in the same format
        if outdir:
            print 'Re-writing',filepath
            fullfin[usedvars] = fin[usedvars]       # This is because some spikes were removed during the process
            fullfin.to_csv(join(outdir, basename(filepath)),
                       header=datalogger_config.header_lines, index=False, quoting=3, na_rep='NaN')
        #-----------------
    
    #-------------
    # We create the summary dataframe
    summary= {k: [len(v)] for k, v in numbers.items()}
    summary=pd.DataFrame({'numbers': map(len,numbers.values())}, index=numbers.keys())
    summary['percent']=100.*summary['numbers']/summary.loc['total','numbers']
    #-------------

    print summary
    summary.to_csv(summary_file, na_rep='NaN')
    return numbers
 

def printUnit(string, mode='L', trim=True, greek=True):
    """
    Returns string formatted for LaTeX or other uses.

    string: string
        string (unambiguous) that represents the unit
    mode: string
        {'L' : 'latex', 'P' : 'pretty', 'H' : 'html' }
    trim: bool
        if true it trims the output
    greek: bool
        yet to be implemented
    """
    try:
        from pint.unit import UnitRegistry
    except ImportError:
        raise ImportError('You must install python-pint in order to use this function')
    ur=UnitRegistry()
    ur.default_format=mode
    u=ur[string]
    u='{:~L}'.format(u)
    if trim:
        u=u[3:].strip()
    return u


def separateFiles(files, dlconfig, outformat='out_%Y-%m-%d_%H:%M.csv', outdir='',
                verbose=False, firstflag='.first', lastflag='.last', save_ram=False,
                frequency='30min', quoting=0, use_edges=False):
    """
    Separates files into (default) 30-minute smaller files. Useful for output files such
    as the ones by Campbell Sci, that can have days of data in one single file.

    Parameters:
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

    Returns:
    --------
    None
    """
    from os import path
    import io
    import pandas as pd
    from algs import auxiliar as algs

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
                if verbose: print 'Writting chunk to ',out
                chunk.to_csv(out, index=False, header=False)
        return
    
    #------------
    # If RAM consumption is an issue, then I won't use pandas. This saves RAM but takes a lot longer.
    else:
        header = dlconfig.header_lines
        if header==None: header=0
    #------------
        for fin in files:
            print 'Now opening',fin
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
                            print labeldates[0]
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
                    print 'Concatenating ', last, ' and ', first
                with open(root, 'wt') as fou:
                    fou.writelines(last_lines)
                    fou.writelines(first_lines)

        print 'List of bad files:', badfiles
        if verbose:
            print 'Done!'
        return


def correctDrift(drifted, correct_drifted_vars=None, correct=None,
                get_fit=True, write_fit=True, fit_file='correctDrift_linfit.params',
                apply_fit=True, show_plot=False, return_plot=False, units={}, return_index=False):
    """
    Parameters:
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

    Returns:
    --------
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
        for slw, fst in rwvars.iteritems():
            slow=correct[slw]
            fast=drifted[fst]
            #----------------
            # Check to see if the frequency in both datasets are the same. Otherwise we are comparing different things
            try:
                if pd.infer_freq(correct.index) == pd.infer_freq(drifted.index):
                    slow, fast = map(np.array, [slow, fast] )
                else:
                    print 'Frequencies must be the same, however, inferred frequencies appear to be different. Plese check.'
            except TypeError:
                print 'Cannot determine if frequencies are the same. We will continue but you should check'
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
        print cors
    #----------------

        #----------------
        # Writes the fit parameters in a file to be used later
        if write_fit:
            cors.index.name='correct_drifted'
            cors.to_csv(fit_file, index=True)
            if units:
                with open(fit_file+'.units', 'wt') as fou:
                    for key, item in units.iteritems():
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
        for slw, fst in rwvars.iteritems():
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
