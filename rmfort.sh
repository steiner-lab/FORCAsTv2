#!/bin/bash

string="fort.3"
echo $string
for i in {01..99} 
do
S="$string$i"
echo $S
rm $S
done


