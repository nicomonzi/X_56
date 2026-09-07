#!/bin/sh

awk '$1==99999 {printf "%+15.8e %+15.8e %+15.8e\n", .8e-2*i++,$2,$8}' ./x.mov > xxx
awk '$1==0 {printf "%+15.8e %+15.8e %+15.8e %+15.8e\n", J, 0, 0, 0; J=-35.835} $1<99999 {J += 1*$4*9.81 + .5*($8*$8+$10*$10) + .5/12*$12*$12}' ./x.mov > yyy
paste xxx yyy > 6barmech.txt
rm -f xxx yyy
