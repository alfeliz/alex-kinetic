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
## @deftypefn {Function File} {} deri ()
##   Derivates numerically a vector, given the space to be used.
##
## @code{der = deri (@var{x},@var{h})} returns the derivative of the vector x, 
## with h the space to be used on the derivative.
## @end deftypefn

## Author: Gonzalo Rodríguez Prieto (gonzalo.rprietoATuclm.es)
## Created: May 2013
##############################################################################################################
##
##                                   FUNCTION TO FIND THE NUMERICAL DERIVATIVE 
##                              Made from algorithm in:
## http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
##
##
##            It shlould be cited as:
##@misc{snrd, author = {Pavel Holoborodko},title = {Smooth Noise Robust Differentiators},howpublished = {http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/}, year = {2008} }
##############################################################################################################


function der = deri(x,h)

   if (nargin!=2) #If not enough parameters are given.
     error("deri: 2 parameters are needed.");
   endif

   if (isscalar(x) ==0) #Works when a vector is passed, not a escalar.
     for i = 3 : rows(x)-2
        der(i) = (  2*(x(i+1) - x(i-1)) + (x(i+2) - x(i-2)) ) ./ ( 8 * h);
     endfor;
   else
     error("deri: the first parameter must be a vector!");
   endif

endfunction
