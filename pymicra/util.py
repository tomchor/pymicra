#!/usr/bin/python
"""
Author: Tomas Chor
Date: 2015-08-07
-------------------------

Modifications:
CHECKLIST
-INCLUDE MAYBE DECODIFICAION OF DATA

"""
import algs
import data
import io
import pandas as pd
import numpy as np


def check_spikes(dfs, visualize=False, vis_col=1, interp_limit=3,
                 f=lambda x: (abs(x - x.mean()) > abs(x.std()*4.)) ):
    '''
    Applies spikes-check accorgin to Vickers and Mart (1997)

    Parameters:
    -----------

    dfs: list, tuple
        sequence of pandas.DataFrame objects
    visualize: bool
        whether of not to visualize the interpolation ocurring
    vis_col: str, int or list
        the column(s) to visualize when seeing the interpolation (only effective if visualize==True)
    '''
    valid_cols=pd.Series(0, index=dfs[0].columns)
    for i in range(len(dfs)):
        chunk=dfs[i].copy()
        if len(chunk)>interp_limit:
            chunk=algs.limitedSubs(chunk, max_interp=interp_limit, func=f)
        valid_cols= valid_cols + chunk.count()

        if visualize:
            aux=dfs[i][vis_col]
            aux2=(abs(aux - aux.mean()) < abs(aux.std()*cut_coef))
            aux.plot(marker='o', color='green', label='original')
            aux.mask(aux2).plot(marker='o', color='red')
        #Interpolation takes place here
        chunk=chunk.interpolate(method='time', axis=0)
        #
        if visualize:
            chunk[vis_col].plot(marker='', color='blue', label='controled')
            plt.show()
        dfs[i]=chunk.copy()

    fou=pd.concat(dfs)
    fou=fou.interpolate(method='time', axis=0)
    fou_points=len(fou.index)
    valid_cols=valid_cols/fou_points
    return fou, valid_cols


def qcontrol(files, datalogger_config,
             bdate='2003-01-01 00:00', edate='2023-12-31 20:00',
             spikes_func=None, interp_limit=3, accepted_percent=1,
             window_size=900, chunk_size='2Min', RATvars=None,
             trueverbose=False, falseverbose=False, falseshow=0, trueshow=0, 
             outdir='quality_controlled', date_format='"%Y-%m-%d %H:%M:%S.%f"',
             std_limits={}, dif_limits={}, low_limits={}, upp_limits={},
             summary='qcontrol_summary.csv'):

    """Program that applies quality control to a set of datafiles

    Parameters:
    -----------

    files: list
        list of filepaths to consider
    datalogger_config: pymicra.dataloger_configuration
        configuration used for all files in files
    bdate: str
        first date to be considered
    edate: str
        last date to be considered
    """
    from io import timeSeries
    from os.path import basename, join
    from itertools import izip_longest
    from datetime import datetime
    from dateutil.parser import parse

    bdate=parse(bdate)
    edate=parse(edate)

    tables=pd.DataFrame(dif_limits, index=['dif_limits'])
    tables.loc['std_limits']=std_limits
    tables.loc['low_limits']=low_limits
    tables.loc['upp_limits']=upp_limits

    numbers={'total': [],
    'spikes': [],
    'STD': [],
    'RAT': [],
    'difference': [],
    'lowest value': [],
    'highest value': [],
    'successful': []}

    variables_list=datalogger_config.varNames
    usedvars=[v for v in variables_list[1:] if r'%' not in v]
    
    
    #-------------------------------------
    # BEGINNING OF MAIN PROGRAM
    #-------------------------------------
    filename_format=datalogger_config.filename_format
    for filepath in files:
        print
        filename=basename(filepath)
        #--------------------------------
        # BEGINNING OF DATE CHECK
        #-------------------------------
        f=''.join([ s for s,v in izip_longest(filename, filename_format) if v!='?' ])
        fmt=filename_format.replace('?','')
        cdate=datetime.strptime(f, fmt)
        print cdate
        if cdate<bdate or cdate>edate:
            print cdate,':', filename, 'discarded because date is not valid'
            continue
        else:
            if trueverbose: print cdate,':', filename, 'included because date is valid'
            pass
    
        #-------------------------------
        # BEGINNING LINE NUMBERS TEST
        #-------------------------------
        result=algs.check_numlines(filepath, numlines=18000)
        if result == False:
            print filepath,'was skipped for not having the correct number of lines!'
            continue
    
        # TRY-EXCEPT IS A SAFETY NET BECAUSE OF THE POOR DECODING (2015-06-21 00:00 appears as 2015-06-20 24:00)
        try:
            fin=timeSeries(filepath, datalogger_config, correct_fracs=True, complete_zeroes='%H%M')
        except ValueError, e:
            if str(e)=='unconverted data remains: 0' and cdate.hour==23:
                continue
            else:
                raise ValueError, e
        fin=fin[usedvars]      # exclude unnused variables
        numbers['total'].append(filename)

        #-------------------------------
        # BEGINNING OF LOWEST VALUE CHECK
        #-------------------------------
        valid= ~(fin < tables.loc['low_limits']).any(axis=0)
        result, failed=algs.testValid(valid, testname='lowest value', trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
        numbers=algs.applyResult(result, failed, fin, control=numbers, testname='lowest value', filename=filename, falseshow=falseshow)
        if result==False: continue
    
        #-------------------------------
        # BEGINNING OF HIGHEST VALUE CHECK
        #-------------------------------
        valid= ~(fin > tables.loc['upp_limits']).any(axis=0)

        result, failed=algs.testValid(valid, testname='highest value', trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
        numbers=algs.applyResult(result, failed, fin, control=numbers, testname='highest value', filename=filename, falseshow=falseshow)
        if result==False: continue

        #-------------------------------
        # BEGINNING OF SPIKES CHECK
        #-------------------------------
        chunks=algs.splitData(fin, chunk_size)
        fin,valid_cols=check_spikes(chunks, visualize=False, vis_col='u', f=spikes_func, interp_limit=interp_limit)
        valid= valid_cols >= (1.-(accepted_percent/100.))

        result, failed=algs.testValid(valid, testname='spikes', trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
        numbers=algs.applyResult(result, failed, fin, control=numbers, testname='spikes', filename=filename, falseshow=falseshow)
        if result==False: continue

        #----------------------------------
        # BEGINNING OF STANDARD DEVIATION CHECK
        #----------------------------------
        df=fin.copy()
        df=df-pd.rolling_mean(df,window=window_size, center=True)
        stds_list=df.resample(chunk_size, np.std).dropna()
        valid= ~(stds_list<tables.loc['std_limits']).any(axis=0)

        result,failed=algs.testValid(valid, testname='STD', trueverbose=trueverbose, filepath=filepath, falseverbose=falseverbose)
        if falseverbose:
            print (not result)*'The failed variables and times:\n{0}'.format(~(stds_list<tables.loc['std_limits']))
        numbers=algs.applyResult(result, failed, fin, control=numbers, testname='STD', filename=filename, falseshow=falseshow)
        if result==False:
            continue
    
        #---------------------------------
        # BEGINNING OF REVERSE ARRANGEMENT TEST
        #---------------------------------
        if RATvars:
            valid_chunks= fin[RATvars].apply(data.reverse_arrangement, axis=0, points_number=50, alpha=.01)
        elif RATvars==None:
            valid_chunks= fin.apply(data.reverse_arrangement, axis=0, points_number=50, alpha=.01)
        else:
            valid_chunks= fin.any(axis=0)

        result,failed=algs.testValid(valid_chunks, testname='reverse arrangement', trueverbose=trueverbose, filepath=filepath)
        numbers=algs.applyResult(result, failed, fin, control=numbers, testname='RAT', filename=filename, falseshow=falseshow)
        if result==False:
            continue
    
        #--------------------------------
        # BEGINNING OF MAXIMUM DIFFERENCE METHOD
        #--------------------------------
        trend=data.trend(fin, mode='linear')
        maxdif=trend.iloc[0]-trend.iloc[-1]
        maxdif=maxdif.abs()
        chunks_valid= tables.loc['dif_limits'] - maxdif
        chunks_valid= ~(chunks_valid < 0)

        result,failed=algs.testValid(chunks_valid, testname='maximum difference', trueverbose=trueverbose, filepath=filepath)
        numbers=algs.applyResult(result, failed, fin, control=numbers, testname='maximum difference', filename=filename, falseshow=falseshow)
        if result==False: continue
    
        #--------------------------------
        # END OF TESTS
        #--------------------------------
        print 'Successful run!'
        if trueshow:
            fin.plot()
            plt.show()
        numbers['successful'].append(filename)
        #--------------------------------
        # FINALLY
        #--------------------------------
        print 'Re-writng',filepath
        fin.to_csv(join( outdir, basename(filepath) ),
                   header=datalogger_config.header_lines, date_format=date_format, quoting=3, na_rep='NaN')
    
    summary= {k: [len(v)] for k, v in numbers.items()}
    summary=pd.DataFrame({'numbers': map(len,numbers.values())}, index=numbers.keys())
    summary['percent']=summary['numbers']/summary.loc['total','numbers']
    print summary
    summary.to_csv('qcontrol_summary.csv', na_rep='NaN')
    return summary
 

def printUnit(string, mode='L', trim=True, greek=True):
    """
    string: string
        string (unambiguous) that represents the unit
    mode: string
        {'L' : 'latex', 'P' : 'pretty', 'H' : 'html' }
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
                frequency='30min', quoting=0):
    """
    Separates files into (default) 30 minute smaller files

    ADD FIRST AND LAST LABELS TO PANDAS METHOD
    NEEDS TO ADD CHUNK-GLUEING SO AS NOT TO LOSE ANY DATES
    """
    from os import path
    #-----------------------
    # First considerations
    outpath=path.join(outdir, outformat)
    firstpath= path.join(outdir, outformat+firstflag)
    lastpath = path.join(outdir, outformat+ lastflag)
    parser = lambda x: algs.line2date(x, dlconfig)
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
    # If RAM consumption is an issue, then I won't use pandas, which saves RAM but takes longer
    else:
        header = dlconfig.header_lines
        if header==None: header=0
    #------------
        for fin in files:
            
            #------------
            # Creates a sequence of equaly-spaced dates based on the frequency and nicely rounded-up
            ft, lt = algs.first_last(fin)
            ft, lt = map(parser, [ft, lt])
            labeldates = pd.Series(index=pd.date_range(start=ft, end=lt, freq='min')).resample(frequency).index
            nfiles=len(labeldates)
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
                        print cdate
                        print labeldates[0]
                        print
                        #------------
                        # labels first and last runs separately, so they can be identified and put together later
                        if len(labeldates) == nfiles:
                            fou = open((labeldates[0]).strftime(firstpath), 'wt')
                        elif len(labeldates) == 1:
                            fou = open((labeldates[0]).strftime(lastpath), 'wt')
                        else:
                            fou = open((labeldates[0]).strftime(outpath), 'wt')
                        #------------
                        fou.writelines(lines)
                        fou.close()
                        labeldates=labeldates.drop(labeldates[0])
                        #---------
                        # this is in case the data has a big jump of more than the value of the frequency
                        while True:
                            if cdate>=labeldates[1]:
                                labeldates=labeldates.drop(labeldates[0])
                            elif (cdate>=labeldates[0]) and (cdate<labeldates[1]):
                                break
                            else:
                                raise ValueError('SOMETHING WRONG')
                        #---------

                        #---------
                        # Starts over the lines
                        lines = [line]
                        #---------
                fou=open((labeldates[0]).strftime(lastpath), 'wt')
                fou.writelines(lines)
                fou.close()
        print 'Done!'
            
            
        #for filename in files:
        #    df=pd.read_csv(filename, header=dlConfig.header, index_col=None, columns=dlConfig.varNames)
        return


