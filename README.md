# Description

This is the code used to generate results presented in "173rd Meeting of the Acoustical Society of America and the 8th Forum Acusticum", as "Numerical study of nonlinear distortion of resonant states in acoustic tubes". Please note this is a part of the project concerning only with the excitation of the system.

Authors : Roberto Velasco-Segura and Pablo L. Rendón

Affiliations : Grupo de Acústica y Vibraciones, Centro de Ciencias
Aplicadas y Desarrollo Tecnológico (CCADET), Universidad Nacional
Autónoma de México (UNAM).

Registration numbers: 03-2014-022811451900-01 INDAUTOR México

Source repository : https://github.com/rvelseg/ASA17

# Some technical info

Expected usage is like
```Octave
> study
```
which takes about 10 hours to complete. It is also possible, and handy when making changes to the code, to execute like
```Octave
> study(24)
```
which will start in run 24.

The whole presented study has 6 calls to the script `study` see details in
the script.

The code has been tested most extensivelly in Octave 4.2 over
GNU/Linux/Debian/Jessie, but also in matlab 2016b, over Jessie, and
over windows 7.

If using Octave the "signal" package is needed.

# Known issues

To generate videos matlab is needed, by now. I expect to fix this soon
for Octave in Linux.

If you are working in Octave, and make changes to `propagator.m`, before each call to `study` you will need to
```Octave
> clear all;
```
or something similar, in order to load again the new deffinition of the class.
