# Supplemental Information - physiological model #

All figures based on simulations with the EMD model can be generated
with this code.  It will run in the virtual machine defined in our
[Vagrantfile](https://github.com/strawlab/asymmetric-motion/blob/master/Vagrantfile).
Everything is described in detail there.

We recommend running this code natively on a multiprocessor machine
with all required packages installed.

On a 24 core (2x Intel(R) Xeon(R) CPU E5-2630 v2 @ 2.60GHz)
machine generating all data takes about 3 hours.

If you intend to run it on your own machine, you need to
fulfill the following requirements:


## Requirements ##


### Generate figures ###

Make sure you have the following installed:

 * Python version >=2.7.3, <3.x
 * numpy >=1.7.1
 * matplotlib >=1.3.0
 * h5py >=2.0.1
 * scipy >=0.9.0
 * progressbar >=2.3

Generate figure data with:
```
python calculate_figure_data.py
```

When the figure data is calculated, display the figures interactively with:
```
python plot_figure3B_bottom.py
python plot_figure3B_top.py
python plot_figureS3.py
```

