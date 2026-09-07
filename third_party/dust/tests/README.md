This folder contains a rudimentary set of regression tests to check the
consistency of the results over some cases after modifications.

The set is not yet complete, and the execution not completely automated.

At the moment different python scripts are employed for the different tests and
must be envoked with:
~~~python
python run_x.py exe_location
~~~
providing the location of the executables to test (both preprocessor and solver are tested at the moment, but not the postprocessor).

Should work both with python 2 and python 3 although python 3 is preferred. 
In addition a few python packages are required to execute the scripts. 

The error with respect to the reference results is printed to screen. Note that due to the parallel execution and compilers optimization, a non-zero error is expected, however the error should remain tiny (more or less around 1e-10).

| |DUST regression tests| | | |.in| | |
|:----|:----|:----|:----|:----|:----|:----|:----|
| | | | |status|case|h5ref|python|
|1|mirror|static|pan| | | | |
|2| | |ll| | | | |
|3| | |vl| | | | |
|4| | |nl_vl| | | | |
|6| |dynamic|pan| | | | |
|7| | |ll| | | | |
|8| | |vl| | | | |
|9| | |nl_vl| | | | |
|11| |coupled|pan| | | | |
|12| | |ll| | | | |
|13| | |vl| | | | |
|14| | |nl_vl| | | | |
|16|symmetry|static|pan| | | | |
|17| | |ll| | | | |
|18| | |vl| | | | |
|19| | |nl_vl| | | | |
|21| |dynamic|pan| | | | |
|22| | |ll| | | | |
|23| | |vl| | | | |
|24| | |nl_vl| | | | |
|26| |coupled|pan| | | | |
|27| | |ll| | | | |
|28| | |vl| | | | |
|29| | |nl_vl| | | | |
|31|plunge|dynamic|pan| | | | |
|32| | |ll| | | | |
|33| | |vl| | | | |
|34| | |nl_vl| | | | |
|35| |coupled|pan| | | | |
|36| | |ll| | | | |
|37| | |vl| | | | |
|38| | |nl_vl| | | | |
|39|1blade|dynamic|pan| | | | |
|40| | |ll| | | | |
|41| | |vl| | | | |
|42| | |nl_vl| | | | |
|43| |coupled|pan| | | | |
|44| | |ll| | | | |
|45| | |vl| | | | |
|46| | |nl_vl| | | | |
|47|hinge|dust|pan| | | | |
|48| | |vl| | | | |
|49| |coupled|pan| | | | |
|50| | |vl| | | | |
