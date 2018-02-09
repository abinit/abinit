#!/bin/bash

mv $1 tmp
exec > $1
sed '1,3!d' tmp
sed '4,/^\*/!d' tmp | head -n-1 | sort
tail -3 tmp
rm tmp
