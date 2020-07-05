#!/bin/tcsh

foreach file (`cat list`)
    git rm ${file}.rst
end
