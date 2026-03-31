#!/bin/bash

# remove unicode characters in make.log for abirules
#

sed -i -e "s/['"‘’"']/'/g" make.log
