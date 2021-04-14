# python_visualisation
Meteorological data visualisation scripts <br>
language: Python, bash
<br>

## 1. Hough_wind_field.py
<table>
  <tr><td>Description: </td><td> plotting wind fields of Hough harmonics for MODES website (outcome figures: https://modes.cen.uni-hamburg.de/hough)</td><tr>
  <tr><td> Input: </td><td> reanalysis wind data (u,v) and geopotential height data (Z) for different waves, decomposed by normal-mode function software MODES (https://modes.cen.uni-hamburg.de/) </td><tr>
  <tr><td>Output: </td><td> wind vectors and geopotential height contour on a global map (longitude range shifted automatically according to the wave center in each data file) </td><tr>
</table>

## 2. mean_spectra.py
<table>
  <tr><td> Description: </td><td> calculate the mean energy spectra from single spectra over a week, and plot them (1-D);
<br>calculate the irrot./total ratio and plot it </td><tr>
  <tr><td> Input: </td><td> binary data of En_k(n,l), En_n(k,l) </td><tr>
  <tr><td> Output: </td><td> 1-D mean energy spectra figures (x: log(k), log(n), y: log(En), irrot./total ratio </td><tr>
</table>


## 3. plot_2D.py
<table>
  <tr><td> Description: </td><td> calculate the mean spectra from single spectra over a week, and plot them (2-D) </td><tr>
  <tr><td> Input: </td><td> binary data of En(k,n,l) </td><tr>
  <tr><td> Output: </td><td> 2-D mean energy spectra figures </td><tr>
</table>

## 4. loop_plot.sh
<table>
  <tr><td> Description: </td><td> Running a python script for 6-hourly data files </td><tr>
</table>

## 5. Rossby_contourf_sph.py, convert_gif.py
<table>
  <tr><td> Description: </td><td> plot contourf of Rossby wave; convert a series of Rossby wave images into a .gif </td><tr>
  <tr><td> Input: </td><td> decomposed Rossby wave data (6-hourly) </td><tr>
  <tr><td> Output: </td><td> see "Gallery-Rossby.gif" </td><tr>
</table>

## 6. Gallery
examples of outcome figures
