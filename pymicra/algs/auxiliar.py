from __future__ import print_function
"""
"""



def stripDown(str, final='', args=['_', '-']):
    """
    Auxiliar function to strip down keywords from symbols
    """
    for arg in args:
        str=str.replace(arg,final)
    return str
        

def lenYear(year):
    """
    Calculates the length of a year in days
    Useful to figure out if a certain year is a leap year
    """
    import calendar
    feblen=calendar.monthrange(year,2)[1]
    otherlens=365-28
    return feblen + otherlens


def completeHM(string):
    """
    Completes %H%M strings for cases when 2hours 0 minutes appear
    as 020. Should be dropped eventually because this is pretty much a hack that
    corrects for file configuration
    """
    if string.isdigit():
        pass
    else:
        raise TypeError('String passed must contain only digits. Check the argument')
    if len(string)==3:
        string='0'+string
    return string



def testValid(df_valid, testname='', falseverbose=True, trueverbose=True, filepath=None):
    """
    Tests a boolean DataFrane obtained from the test and prints standard output

    Parameters
    -----------
    df_valid: pandas.Series
        series contaning only True or False values for each of the variables, which should be the indexes
    testname: string
        the name of the test that generated the True/False values
    falseverbose: bool
        whether to return which variables caused a false result
    trueverbose: bool
        whether to print something successful cases

    Returns
    --------
    result: bool
        True if the run passed the passed
    failed: list
        list of failed variables if result==False. None otherwise.
    """
    from os.path import basename
    
    if False in df_valid.values:
        failed = df_valid[ df_valid==False ].index
        print(basename(filepath), ': !FAILED',testname,'test!\n')
        if falseverbose:
            print('Failed variable(s):', ', '.join(failed),'\n')
            print
        return False, failed
    else:
        if trueverbose: print(basename(filepath),'passed',testname,'test')
        return True, None



def applyResult(result, failed, df, control=None, testname=None, filename=None, falseshow=False, index_n=None):
    """
    Auxiliar function to be used with util.qcontrol

    Parameters
    -----------
    result: bool
        whether the test failed and succeeded
    failed: list
        list of failed variables. None object if the test was successful
    control: dictionary
        dictionary whose keys are the names of the tests and items are lists
    testname: string
        name of the test (has to match control dict)
    filename: string
        name or path or identifier of the file tested
    falseshow: bool
        whether to show the failed variables or not
    """
    import matplotlib.pyplot as plt
    if result==False:
        if falseshow:
            df[failed].plot()
            plt.show()
        if type(control)==dict:
            control[testname].append(filename)
        else:
            control.loc[ index_n, testname ] = filename
    return control


def first_last(fname):
    """
    Returns first and last lines of a file
    """
    with open(fname, 'rb') as fin:
        first=fin.readline()
        for line in fin:
            pass
        last=line
    return first, last

