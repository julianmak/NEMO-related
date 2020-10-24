#!/bin/bash
#! Script to clean up the NEMO outputs

echo " "
echo "Christopher is going to run the model for you ......"
echo "                  ,-.____,-.          "
echo "                  /   ..   \          "
echo "                 /_        _\         "
echo "                |'o'      'o'|        "
echo "               / ____________ \       "
echo "             , ,'    '--'    '. .     "
echo "            _| |              | |_    "
echo "          /  ' '              ' '  \  "
echo "         (    ',',__________.','    ) "
echo "          \_    ' ._______, '     _/  "
echo "             |                  |     "
echo "             |    ,-.    ,-.    |     "
echo "              \      ).,(      /      "
echo "         gpyy   \___/    \___/        "


# remove the existing transport files otherwise diadct complains
rm -rv *transport

mpiexec -n 1 ./opa &
