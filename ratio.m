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
## @deftypefn {Function File} {} ratio (@var{x},@var{v}, @var{opt})
##   Calculates the ratio between kinetic energy and plasma expansion work in a plasma and other useful functions based on the dimensionless radial expansion @var{x} and the velocity @var{v} following different options given by @var{opt}.
##
## @code{rela = ratio (@var{x},0,1)} returns the ratio between kinetic energy 
## and plasma expansion work in a plasma without the dimensions.
## @code{rela = ratio (@var{x},0,2)} returns the plasma expansion working energy 
## calculated from the radius expansion without the proper dimensions.
## @code{rela = ratio (@var{x},0,3)} returns the plasma kinetic energy 
## calculated from the radius expansion without the proper dimensions.
## @code{rela = ratio (@var{x},@var{v},4)} returns the plasma internal energy from the Mie-Gruneisen state equation.
## @code{rela = ratio (@var{x},0,5)} returns the plasma internal energy from dimensionless simplified state eq. of paper.
## @end deftypefn

function rela = ratio(r,v, opt)

if ( (v==0) && (opt==4) )
 error("ratio: Velocity vector is cero, so no possible to use in State Equation.");
endif;

%Kinetic energy:
kin = r.**2 -1;
%Plasma expansion:
volexp = ( (1 - r).**2./120 ) .* ( 24.*(r.**3) + 18.*(r.**2) + 12.*r + 6 );
%Mie-Gruneisen state equation energy:
nu = 1./r; %Densities ratio
gamma0 = 2.0; %A constant.
rho0 = 8960; %Kg/m³. Copper density.
s =1.5; %The linear Hugoniot slope coefficient.
c0 = 3933; %m/s Sound velocity on the metal.
a = nu - ( 0.5 .* gamma0 .* (nu -1 ) );
b = nu - ( s .* (nu-1) );
part1 = (1-s) .* 2 .* (gamma0 - 1) .* b;
mie = (1./part1) .* ( (v.**2 .* b) - (rho0 .* c0.**2 .* a) - (rho0 .* c0.**2 .* (1 - 0.5.*gamma0) .* (nu -1) )  + ( 2 .* (1-s) .* rho0 .* c0.**2 .* (nu-1) .*(a./b) ) ); 
%Dimensionless simplified EOS:
Bi = 140e6; %Pascals (Copper bulk modulus)
B = 2; %Constant between 2 and 4.
rini = 0.05 .* 1e-3; %For ALEX095 (m)
l = 2e-2; %m. Wire length
e1 = ( pi .* rini.**2 .* l .* Bi) ./  (B-1); 
simplyene = e1 .*  r.**2 .* ( (1./B).*(1./r).**(2.*B) -(1./B) -(1./r).**2 +1  );



switch opt
  case {1}
    rela = kin ./ volexp;
  case {2}
    rela = volexp;
  case {3}
    rela = kin;
  case {4}
    rela = mie;
  case {5}
    rela = simplyene;
endswitch;


endfunction

