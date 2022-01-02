#!/usr/bin/bash

doxygen ../doc/documentation.dox

cd ../doc/latex/

make

cp refman.pdf ../../fluid_simulation.pdf
