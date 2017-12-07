

Quality control procedures
--------------------------

Here we give a brief example of how to perform some simple quality control
procedures using Pymicra, which is generally the first thing to be done with a
new dataset on our hands. We will use some 1-hour files as examples for this
procedure.

.. ipython:: python
      :suppress:

   import pymicra as pm
   pm.notation = pm.Notation()
   import pandas as pd
   pd.options.display.max_rows = 30



The first thing is to read the files so as to have a list of file paths. We
will need a configuration file since we must know what is in each column to be
able to perform the quality control correctly.

.. ipython:: python

   import pymicra as pm
   from glob import glob
   fconfig = pm.fileConfig('../examples/lake.config')
   fnames = glob('../examples/ex_data/***.csv')
   print(fnames)


There are two functions that currently do the quality control in Pymicra. The
first one, `pm.util.qc_replace()`, is intended for quality control procedures
that can "fix" certain data by removing or interpolating a few "bad points". At
this time it checks for file dates and file lines and it filters the data for
spikes, NaNs and out-of-bounds values. Every time one of those latter one is
found they are interpolated with either the linear trend of the dataset or a
simple linear interpolation.  If at the end the number of interpolations is
higher than a user-defined threshold than the run fails the quality control and
is discarded. Thus we can do the first part of the quality control with


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
 that falls outside theese bounds is interpolated.  Usually it is suggested that
 at least the concentrations have `0` as lower limits. Or that very high values
 should be used that would surely indicate sensor malfunction (like temperatures
 of over 70 degrees Celsius).
- `spikes_test` tells the function whether to test for spikes or not. Spikes are interpolated over.
- `spikes_func` every point is tested for spikes using this function. If the result if `True`, it
 is considered a spike. Any Pandas function works here, as `x` in this case is a pandas `DataFrame`.
- `visualize_spikes` decides if you want to see the points you are considering as spikes or not. A matplotlib
 plot appears on screen. This is good for the first iterations of the quality control, so you can calibrate your
 spike parameters and see if your spikes test is really doing what you think it's doing.
- `max_replacement_count` sets the maximum number of replaces values a run can have before it is discarded.
 This included replacements from the NaNs test, bounds test and spikes test.
- `outdir` is the path to the directory where the quality-controled files will be copied.
- `replaced_report` is the path to a file that will be created with a detailed report
 on the replacements that were made.


More options are available


