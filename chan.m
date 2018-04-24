## Copyright (C) 2008-2012 Ben Abbott
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {} chan ()
##   Returns the value of the data in a CSV Tektronik scope file named 
## "channame"
##
## @code{[t,vol,points] = chan (@var{channame})} returns two vectors, t and vol, 
## with the time and voltage data and the number of points in the file (points)
## @end deftypefn

## Author: Gonzalo Rodríguez Prieto (gonzalo#rprietoATuclm#es)
##       (Mail: change "#" by "." and "AT" by "@")
## Created: May 2013

function [t,vol,points] = chan(channame)

if (ischar(channame)!=1)
 error("chan: The name must be a string!");
endif;

#Reading the file and placing the data in a "cell" structure. 
[fileid, msg] = fopen(channame,'r'); 

if (fileid == -1) 
   error ("chan: Unable to open file name: %s, %s",channame, msg); 
endif; 
r_rows =  "%s %f %f %f %f";
if feof(fileid)==0 #Read the file until it finds the EndOfFile character.
   data = textscan(fileid, r_rows, "delimiter", ",", "endofline", "");
endif;
fclose(fileid);

#Transforming the data in numbers:
t = cell2mat(data(:,4)); #Time vector (in seconds)

vol = cell2mat(data(:,5)); #Voltage channel value (in volts)
vol = supsmu(t,vol,"span",0.01); #Smoothing the voltage to work with it.

points = length(vol);


endfunction;

#That's...that's all folks!!!
