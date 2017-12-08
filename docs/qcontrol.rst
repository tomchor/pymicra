

Quality control procedures
--------------------------

Here we give a brief example of how to perform some simple quality control
procedures using Pymicra, which is generally the first thing to be done with a
new dataset on our hands. We will use some 1-hour files as examples for this
procedure. The following lines simply load the necessary packages.

.. ipython:: python
      :suppress:

   import pymicra as pm
   pm.notation = pm.Notation()
   import pandas as pd
   pd.options.display.max_rows = 30 # This line just alters the visualization



The first thing is to read the files so as to have a list of file paths. We
will need a configuration file since we must know what is in each column to be
able to perform the quality control correctly.

.. ipython:: python

   import pymicra as pm
   from glob import glob
   fconfig = pm.fileConfig('../examples/lake.config')
   fnames = glob('../examples/ex_data/***.csv')
   print(fnames)


There are two functions that currently do the quality control in Pymicra:
`pm.util.qc_replace()` and `pm.util.qc_discard()`. They are separated only
for simplicity, but are meant to be used in sequence (first `qc_replace()` and
then `qc_discard()`).

The first one, `pm.util.qc_replace()`, is intended for quality control
procedures that can "fix" certain data by removing or interpolating a few "bad
points". At this time it checks for file dates and file lines and it filters
the data for spikes, NaNs and out-of-bounds values. Every time one of those
latter one is found they are interpolated with either the linear trend of the
dataset or a simple linear interpolation.  If at the end the number of
interpolations is higher than a user-defined threshold than the run fails the
quality control and is discarded. Thus we can do the first part of the quality
control with


.. ipython:: python

    pm.util.qc_replace(fnames, fconfig,
        file_lines=72000,
        nans_test=True,
        lower_limits=dict(theta_v=10, mrho_h2o=0, mrho_co2=0),
        upper_limits=dict(theta_v=45),
        spikes_test=True,
        spikes_func=lambda x: (abs(x - x.median()) > (5./0.6745)*x.mad()),
        visualize_spikes=True,
        chunk_size=2400,
        max_replacement_count=1440,
        outdir='../examples/passed_1st',
        replaced_report='../examples/rrep.txt')

This command effectively does the first part of the quality control. There are other options
that aren't used here but let's go through the options we used

- `file_lines` tells the function the amount of lines you expect the files to have. Anything different from that
 number means that there was a problem while collecting the data and that file is discarded
- `nans_test=True` means you want to check for "Not a Number"s, or missing numbers. If one is found it is
 interpolated
- `lower_limits` and `upper_limits` provide the bounds for the data. Any data
 that falls outside these bounds is interpolated.  Usually it is suggested that
 at least the concentrations have `0` as lower limits. In the previous example we considered that
 virtual temperatures lower than 10 C or higher than 45 C, given the site atmospheric conditions.
- `spikes_test` tells the function whether to test for spikes or not. Spikes are interpolated over.
- `spikes_func` every point is tested for spikes using this function. If the result if `True`, it
 is considered a spike. Any Pandas function works here, as `x` in this case is a pandas `DataFrame`.
- `visualize_spikes` decides if you want to see the points you are considering as spikes or not. A matplotlib
 plot appears on screen. This is good for the first iterations of the quality control, so you can calibrate your
 spike parameters and see if your spikes test is really doing what you think it's doing.
- `max_replacement_count` sets the maximum number of replaces values a run can have before it is discarded.
 This included replacements from the NaNs test, bounds test and spikes test.
- `outdir` is the path to the directory where the quality-controlled files will be copied.
- `replaced_report` is the path to a file that will be created with a detailed report on the replacements that were made.

.. note::
 More options are available for the function. Please read the documentation for
 `pymicra.util.qc_replace()` for details.

With this function, all the runs that passed our quality control were fixed for
the spikes and NaNs that were found and were copied to `outdir`. The second
step is to get these files and apply the second part of the quality control
procedure, `pymicra.util.qc_discard()`. This second part applies tests that, if
turn out to fail, there is no way to recover the data and the dataset is simply
discarded. The tests are the Standard Deviation (STD) test and the Maximum
difference (or stationarity) test.

The STD test separates the run into chunks (of length given by `chunk_size`)
and takes the STD of each chunk. If it is smaller than a certain threshold,
than the run is discarded.  The Max. diff. test can test for different things,
depending on the parameters. The most common use perhaps is to use the
`maxdif_trend` keyword to get the low-frequency variations (low-frequency
trend) by doing a moving median and taking out the absolute difference between
the higher and the lower values of this trend. If it is higher than a certain
value the run is considered non-stationary and discarded. It is kept otherwise.
There are more options and more things to do, so be sure to read the docs for
`pymicra.util.qc_discard()`.

Below is an example of a simple quality-control procedure done using the data
that passed the previous procedure.

.. ipython:: python

    fnames = sorted(glob('../examples/passed_1st/***.csv'))
    print(fnames)

    pm.util.qc_discard(fnames, fconfig,
    std_limits = dict(u=0.03, v=0.03, w=0.01, theta_v=0.02),
    std_detrend=dict(how='linear'), 
    chunk_size=1200, 
    dif_limits = dict(u=5.0, v=5.0, theta_v=2.0),
    maxdif_trend=dict(how='movingmedian', window= 600), 
    failverbose=True, 
    failshow=True, 
    outdir='../examples/passed_2nd', 
    summary_file='2nd_filter_summary.csv', 
    full_report='2nd_full_report.csv')


After the previous example you should be left with reasonably
quality-controlled dataset. We used the following options in this case:

- `std_limits` gives the minimum value for the STD before it is considered too small and the run is discarded. So,
 values of STD larger than `std_limits` are fine.
- `std_detrend` gives the detrending options to use before taking the STD in the STD test.
- `chunk_size` size of chunks in which to separate the data for the STD test.
- `dif_limits` are the maximum difference between the largest and lowest value in the Max dif (stationarity) test. If the
 difference is larger the run is discarded.
- `maxdif_trend` is the trend to consider when taking the difference between the lowest and highest value.
- `failverbose` shows extra info on runs that fail any test (like which variable(s) failed).
- `failshow` plots runs that failed some test on the screen for checking purposes.

As usual, there are more options to this, but this is a basic introduction.
With the dataset already corrected and filtered we move on to the processing of
data.

