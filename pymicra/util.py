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






def reverse_arrangement(array, points_number=None, alpha=0.05):
    '''
    Performs the reverse arrangement test
    according to Bendat and Piersol - Random Data - 4th edition, page 96

    Parameters
    ----------

    Array: np.array, list, tuple, generator
    array which to test for the reverse arrangement test

    points_number: integer
    number of chunks to consider to the test. Maximum is the length of the array.
    If it is less, then the number of points will be reduced by application of a mean

    alpha: float
    Significance level for which to apply the test

    WARNING! This fuction approximates table A.6 from Bendat&Piersol as a normal distribution.
    This may no be true, since they do not express which distribution they use to construct
    their table. However, in the range 9<N<101, this approximation is as good as 5% at N=10
    and 0.1% at N=100.
    '''
    if points_number==None:
        points_number=len(array)
        mean_Matrix=array
    elif points_number==len(array):
        mean_Matrix=array
    else:
        pts = len(array)//points_number
        mean_Matrix = []
        for j in range(0,points_number):
            mean_Matrix.append( mean(array[(j*pts):((j+1)*pts)]) ) # Calculo a media de cada um dos 50 intervalos
    A = []
    for i in range(len(mean_Matrix)):
        h = []
        for j in range(i,len(mean_Matrix)):
            if(mean_Matrix[i] > mean_Matrix[j]):
                h.append(1)
        A.append(sum(h))
    Atot = sum(A)
    N=len(mean_Matrix)
    mu=N*(N-1)/4
    variance=N*(2*N+5)*(N-1)/72
    #
    # USE BELOW FUNCTION
    def mu_var(N):
        mu=N*(N-1.)/4.
        variance=N*(2.*N+5.)*(N-1.)/72.
        return mu, variance

    f=inverse_normal_cdf(mu, math.sqrt(variance))
    phi1=1-alpha/2.
    phi2=alpha/2.
    A1=f(phi1)
    A2=f(phi2)
    if A1 < Atot < A2:
        return True
    else:
        return False




def qcontrol(files, datalogger_config,
             bdate='2003-01-01 00:00', edate='2023-12-31 20:00',
             cut_coef=4.0,
             interp_limit=3, accepted_percent=1,
             window_size=900, chunk_size='2Min',
             trueverbose=0, falseshow=0, trueshow=0, 
             outdir='quality_controlled', date_format='"%Y-%m-%d %H:%M:%S.%f"',
             std_limits={}, dif_limits={}, lower_boundaries=None,
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
    from os.path import basename
    from itertools import izip_longest
    from datetime import datetime
    from dateutil.parser import parse

    bdate=parse(bdate)
    edate=parse(edate)

    tables=pd.DataFrame(dif_limits, index=['dif_limits'])
    tables.loc['std_limits']=std_limits

    numbers={'total': [],
    'spikes': [],
    'STD': [],
    'RAT': [],
    'difference': [],
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
            fin=timeSeries(filepath, datalogger_config, correct_fracs=True)
        except ValueError, e:
            if str(e)=='unconverted data remains: 0' and cdate.hour==23:
                continue
            else:
                raise ValueError, e
        fin=fin[usedvars]      # exclude unnused variables
        numbers['total'].append(filename)
    
        #-------------------------------
        # BEGINNING OF SPIKES CHECK
        #-------------------------------
        chunks=algs.splitData(fin, chunk_size)
        fin,valid_cols=check_spikes(chunks, visualize=False, vis_col='u')
        chunks_valid= valid_cols >= (1.-(accepted_percent/100.))
        result, failed=algs.testValid(chunks_valid, testname='spikes', trueverbose=trueverbose)
        if result==False:
            if falseshow:
                fin[failed].plot()
                plt.show()
            numbers['spikes'].append(filename)
            continue
    
        #----------------------------------
        # BEGINNING OF STANDARD DEVIATION CHECK
        #----------------------------------
        df=fin.copy()
        df=df-pd.rolling_mean(df,window=window_size, center=True)
        stds_list=df.resample(chunk_size, np.std).dropna()
        for l,stds in enumerate(stds_list.values):
            std=stds_list.iloc[l]
            lim=tables.loc['std_limits']
            chunks_valid=std-lim
            chunks_valid=chunks_valid>=0
            result,failed=algs.testValid(chunks_valid, testname='STD', trueverbose=0)
            if result==False:
                break
            else:
                pass    
        if result==False:
            if falseshow:
                print 'Failed chunk was:'
                print std
                fin[failed].plot()
                plt.show()
            numbers['STD'].append(filename)
            continue
        else:
            if trueverbose: print filepath,'Passed STD test'
    
        #---------------------------------
        # BEGINNING OF REVERSE ARRANGEMENT TEST
        #---------------------------------
        valid_chunks=fin.apply(reverse_arrangement, points_number=360, axis=0, broadcast=False)
        result,failed=algs.testValid(chunks_valid, testname='reverse arrangement', trueverbose=trueverbose)
        if result==False:
            if falseshow:
                fin[failed].plot()
                plt.show()
            numbers['RAT'].append(filename)
            continue
    
        #--------------------------------
        # BEGINNING OF MAXIMUM DIFFERENCE METHOD
        #--------------------------------
        trend=pm.data.trend(fin, mode='linear')
        maxdif=trend.iloc[0]-trend.iloc[-1]
        maxdif=maxdif.abs()
        chunks_valid= tables.loc['dif_limits'] - maxdif
        chunks_valid=chunks_valid >= 0
        result,failed=algs.testValid(chunks_valid, testname='maximum difference', trueverbose=trueverbose)
        if result==False:
            if falseshow:
                fin[failed].plot()
                plt.show()
            numbers['difference'].append(filename)
            continue
    
        print 'Successful run!'
        if trueshow:
            fin.plot()
            plt.show()
        numbers['successful'].append(filename)
        #--------------------------------
        # FINALLY
        #--------------------------------
        print 'Re-writng',filepath
        fin.to_csv(os.path.join( outdir, os.path.basename(filepath) ),
                   header=datalogger_config.header_lines, date_format=date_format, quoting=3, na_rep='NaN')
    
    summary= {k: [len(v)] for k, v in numbers.items()}
    summary=pd.DataFrame({'numbers': map(len,numbers.values())}, index=numbers.keys())
    summary['percent']=summary['numbers']/summary.loc['total','numbers']
    print summary
    summary.to_csv('qcontrol_summary.csv', na_rep='NaN')
    return
 

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

