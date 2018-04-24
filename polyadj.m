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
## @deftypefn {Function File} {} polyadj ()
##   Adjust to a polynomium a vector, with a minimum correlation factor of 0.1.
##
## @code{padj = polyadj (@var{x},@var{fx},@var{rmin})} returns polynomium that adjusts to the vector
## @var{fx} with @var{x} values with a correlation factor larger than @var{rmin}.
## @code{[padj, rcor, padjstruc] = polyadj (@var{x},@var{fx},@var{rmin})} returns a polynomium that adjusts to the vector
## @var{fx} with @var{x} values with a correlation factor larger than @var{rmin} in @var{padj}, the correlation factor in @var{rcor} and the
## polynomium structure in @var{padjstruc}.
## @end deftypefn

## Author: Gonzalo Rodríguez Prieto (gonzalo.rprietoATuclm.es)
## Created: Nov 2013



function [padj, rcor, padjstruc] = polyadj(x,fx,rmin)


#########
# Control error part
#########
if (nargin!=3) #If not enough parameters are given.
 error("poladj: 3 parameters are needed.");
endif;

if (isscalar(x) ==1) #Works when a vector is passed, not a escalar or other type.
 error("poladj: Variable x must be a vector.");
endif;

if (isscalar(fx) ==1) #Works when a vector is passed, not a escalar or other type.
 error("poladj: Variable fx must be a vector.");
endif;

if ( length(x) != length(fx) ) #Length of both vectors are not the same.
 error("poladj: Lengths of vectors x and fx must coincide.");
endif;

if (rmin<0.01)
 warning('poladj: Value of minimum correlation factor switched to 0.1');
 rmin = 0.1;
endif;


#########
# Adjusting loop part
#########
rcor = 0.01;
n = 1;

#Fitting loop :
# (Use the polyfit function until the correlation factor is larger than rmin
while ( (rcor<rmin) && (n<10) ) 
	#Fitting:
	[padj, padjstruc] = polyfit(x,fx, n);

	#Calculation of the goodnes of the fit:
	denom = (length(fx) -  1) * var(fx);
	rcor = abs(1 - padjstruc.normr.**2 / denom);
        if (rcor>1) %To avoid problems in the correlation factor.
	  rcorv(n) = 0;        
	else
	  rcorv(n) = rcor;
        endif;
        n = n + 1;
endwhile;


if (n>=10) #Only when there is no fitting with corr. factor less than rmin.
  #Polynomium order of the best fit:
  npol = find(rcorv == max(rcorv)); %Find the position of the minumum value of the correlation.
					%(The best polynomial fit of 10)
  [padj, padjstruc] = polyfit(x, fx, npol);
  #Calculation of the goodnes of the fit:
  denom = (length(fx) -  1) * var(fx);
  rcor = 1 - padjstruc.normr.**2/denom;
endif; 


endfunction;
