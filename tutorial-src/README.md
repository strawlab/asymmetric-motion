# Supplemental Information - Phenomenological Model #

The tutorial for the phenomenological model can be viewed and
generated in various different formats. It's recommended to run
everything in the virtual machine defined in our
[Vagrantfile](https://github.com/strawlab/asymmetric-motion/blob/master/Vagrantfile).
Everything is described in detail there.

If you intend to run it on your own machine, you need to
fulfill the following requirements:

## Requirements ##


### Generate figures only ###

Make sure you have the following installed:

 * Python version >=2.7.3, <3.x
 * numpy >=1.7.1
 * matplotlib >=1.3.0
 * sympy >=0.7.1.rc1
 * scipy >=0.9.0

Generate the figures with:
```
python phenomenological_model_tutorial_code.py
```

Display the figure interactively with:
```
python phenomenological_model_tutorial_code.py --show
```


### Generate the PDF from the markdown file ###

You need to have everything installed to create the figures.
And you need additional packages. If you are on _Ubuntu 14.04_,
install the following packages:

 * texlive-fonts-recommended
 * texlive-latex-extra
 * texlive-latex-base
 * inkscape
 * pandoc
 * jq

For other distributions check with your package-manager.
Make sure that your _pandoc_ version is >=1.11.

Then just run:
```
make pdf-standalone
```


### Generate and run the ipython notebook ###

You need to have everything installed to generate the PDF.
And you need additional packages. If you are on _Ubuntu 14.04_,
install the following packages:

 * ipython
 * ipython-notebook

For other distributions check with your package-manager.
Make sure that your _ipython_ version is >=1.0.0.

To generate the ipython notebook and launch the ipython notebook
server, run:
```
make ipynb
ipython notebook
```

Your default browser will open up and display the notebook.

