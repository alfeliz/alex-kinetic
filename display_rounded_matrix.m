## Copyright (C) 2008-2012 Ben Abbott
##
## This file is part of Octave.
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
## @deftypefn {Function File} {} display_rounded_matrix ()
##   Displays a matrix with a engineering format, 
## with different precision for different columns.  
##
## @code{display_rounded_matrix (@var{matrix},@var{significative},@var{output})} 
##  displays the matrix @var{matrix} with @var{significative} cyphers in each column.
##  If not @var{output} is given, the stdout is used.
##  @var{significative} could be a column vector, with the same number of columns as the matrix.
## @end deftypefn

## Modified by: Gonzalo Rodríguez Prieto (gonzalo.rprietoATuclm.es)
## In: July 2013


function display_rounded_matrix(matrix, significative, outputFile)


###
# Checking if parameters are Ok:
###
  
#Adequate number of them
  if (nargin < 2)
	error("display_runded_matrix: Wrong number of parameters, at least matrix and significative must be provided.");
  endif;
#Significative values are a vector:
  if (rows(significative)!=1)
	disp("significative:");
	disp(significative);
	error("display_runded_matrix: Significative vector not a vector");
  endif;


###
# Initializing some variables:
###

#Console exit:
  if (exist("outputFile","var")!=1) #Assuming console.
     outputFile = stdout;
  endif;

#Intercolumns space:
  space_between_columns = " "; #One space character.

#Matrix rows and columns:
  row = rows(matrix);
  cols = columns(matrix);

#Significative ciphers transfer toa cell:
  precision_format = cell(columns(significative), 1);
  for i = 1:columns(significative) #transfering the precision into a cell of numbers.
    precision_format{i,1} = strcat("%.", num2str(significative(1,i)), "e");
  endfor;


###
#Writing the matrix:
###

if (columns(significative) == 1) #One significative value for all the matrix columns.
  formato = strcat("%.", num2str(significative(1,i)), "e"); #Format line. Change it of you want a different output.
  for i = 1:row #For every row in the matrix:
	for j = 1:cols #Write every data column in the row i.
		fprintf(outputFile, sprintf(formato, matrix(i,j)));
		if (j ~= cols) #Whenever we are not at the end of the column, put the space between columns.
		  fprintf(outputFile, space_between_columns);
		endif;
	endfor;
  	fprintf(outputFile, "\n"); #Row finished.
  endfor;
elseif (columns(significative) == columns(matrix)) #Every column has its own number of cifers.
  for i = 1:row #For every row in the matrix:
	for j = 1:cols #Print the column, each with its format (from the cell precision_format)
		  fprintf(outputFile, sprintf(precision_format{j,1}, matrix(i,j)));
		  if (j ~= cols) #Whenever we are not at the end of the column, put the space between columns.
		    fprintf(outputFile, space_between_columns);
		  endif;
	endfor;
  	fprintf(outputFile, "\n"); #Row finished.
  endfor;
  else #In case the function does not understand the significative variable.
	disp("significative:");
	disp(significative);
	error("display_rounded_matrix: A problem with the significative variable.");
endif;



endfunction;

#That's...that's all folks!!!
