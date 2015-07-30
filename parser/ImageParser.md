
Image Parser
================

We rely on a custom C++ programm to convert cell tracking data into a format that is suitable to build a database.

To run the image parser on a cluster you need to emulate X (see [here](http://stackoverflow.com/questions/13215120/how-do-i-make-python-qt-and-webkit-work-on-a-headless-server)

    xvfb-run ./imageParser

You can use the provided binary on linux64. For other platforms you need to compile it manually. Two required libraries need to be installed separately

    sudo apt-get install netcdf-bin
    sudo apt-get install libnetcdf-dev