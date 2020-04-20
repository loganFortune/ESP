# Epilepsy Seizures Prediction (ESP)
Epilepsy Seizures Prediction (ESP) Using Non-Linear Dynamics And Machine Learning.

The code provided is a wrapper for some Tisean tools. It has been optimized for our purpose: __the prediction of epileptic seizures__.

# ESP functions

## EEG Visualization 

![Data Visualization](/images/all_data_skull.png)

## EEG Data Analysis

![Power Spectrum](/images/power_spectrum.png)
![Spectogram](/images/spectogram.png)

## Embedding Delay Analysis (Autocorrelation and Mutual Information)

![Power Spectrum](/images/mutual_info.png)

## Embedding dimension (Taken's theorem) - False nearest neighbors

## Correlation Dimension Analysis

## Maximum Lyapunov exponent

![Maximum Lyapunov Exponent](/images/lyap.png)

## Synchrony Index

# Installation Procedure

## Python Packages

In order to run the ESP software, you must intall some Python Packages :

- _numpy, scipy, matplotlib_
- _pyedflib_:
```
pip install pyEDFlib
```
- _PyGnuplot_:
```
pip install PyGnuplot
```
- _nolds_: 'NOnLinear measures for Dynamical Systems' software is a Python Package for Times Series Analysis (https://cschoel.github.io/nolds/)
```
pip install nolds
```
## Tisean Package

- You must follow the Tisean [installation process](https://www.pks.mpg.de/~tisean/archive_3.0.0.html) according to your operating system. 
- It is easier to download binary executables. 
- Follow [this procedure](https://www.pks.mpg.de/~tisean/Tisean_3.0.1/index.html) for **Linux users**.
- Don't forget to add /bin file in your 'Path' (environment variables) for **Windows users** (_only thing to do_).

You should be able to see that kind of result after the installation:

![Henon -h shell command](/images/install_Tisean.png)

## Gnuplot

- Gnuplot is already installed in Linux environment. 

- For Windows users, [install Gnuplot](https://sourceforge.net/projects/gnuplot/files/). 
**Don't forget to check the** _Add to 'Path' environment variables_ **option !**
