###########################################################################################
#
#
#  OCTAVE SCRIPT TO:
#       - READ ALEX DATA FROM TEKTRONIK SCOPES,
#	- READ DATA FROM STREAK IMAGES,
#	- OBTAIN FROM THEM RADIUS, AND VELOCITY DATA AND
#	- TRANSFORM THESE DATA INTO INTERNAL ENERGY AND ELECTRICAL ENERGY DATA TO COMPARE THEM.  
#    Made by Gonzalo Rodríguez Prieto
#       (gonzalo#rprietoATuclm#es)
#       (Mail: change "#" by "." and "AT" by "@")
#              Version 2.0
#
#
#########################################################
#
#  It uses the functions: 
#                        chan
#                        supsmu
#                        baseline
#			 ratio
#			 deri
#			 poladj
#			 display_rounded_matrix
#  They must be in the same directory.
#
###########################################################################################

more off; %To make the lines display when they appear in the script, not at the end of it.

clear; %Just in case there is some data in memory.

tic; %Total time of the script.


#########################################################################################################
#
#	MAIN PARAMETERS FOR ALL THE SCRIPT:
#
#########################################################################################################

#Streak camera sweep range (For later calibration)
sweep = 097.66; %nanosec/px (For 50 microseconds in 512 rows in the picture.) [50 000 / 512]

#Metal and wire parameters:
rho = 8960; %Kg/m³ Copper density.
lengthwire = 1e-2; %m Wire length
r0 = 0.05e-3; %m (Wire radius)
B = 2; %Constant between 2 and 4 for the EOS
Bi = 140e6; %Pascals (Copper bulk modulus)

#########################################################################################################
#
#	SCOPE DATA INTERPRETATION PART.
#
#########################################################################################################



###
# Reading information over the shot features:
###

#Choosing the shot:
disp("On the ALEX shoot name, use alway capital letters.")
shot = input("Which ALEX shoot to transform?\n","s");
#shot = "ALEX095"; %For testing purposes

#Reading the TXT file with info over the shot:
#Shot name transformation:
infoshot = horzcat(shot,".txt");
#Open it:
[fileinfo, msg] = fopen(infoshot, "r");
if (fileinfo == -1) 
   error ("Script alex-kinetic: Unable to open file name: %s, %s",infoshot, msg); 
endif;
#Putting the information into a cell structure:
txtinfo = textscan(fileinfo,"%s", "treatasempty","/");
fclose(fileinfo); #Closing the file
txtinfo = cell2mat(txtinfo); #To convert the cell structure in a matrix of strings.

#Putting the channels info in matrices of strings:
#(This is valid if the format of file SHOT.txt follows the procedure of ALEX data)
channels = ""; #Channel name
chantype = ""; #Channel probe type
for i=1 : length(txtinfo)
 if (cell2mat(strfind(txtinfo(i),"CH")) ==1) %We have the channels here, with the type following.
   % Standard format of the *.txt file for ALEX.
   channels = strvcat(channels,txtinfo(i));
   chantype = strvcat(chantype,txtinfo(i+1));
 elseif(cell2mat(strfind(txtinfo(i),"Del")) ==1) %We should stop to avoid channels confusion.
   #After channels info cames other shot info, not stored or used in this version of ALEX.m
   break;
 endif;
endfor;


###
#Taking all channel data in one matrix (VOL): 
###

for i=1 : rows(channels)
 channame = strcat(shot,".",channels(i,:),".CSV");
 fprintf(stdout,"Taking data from %s\n", channame);
 [t(:,i),vol(:,i),points(i)] = chan(channame); #chan is the function to take the channel file CHANNAME info.
endfor;


###
#Transforming scope data from (t,volts) into (t,dimensional units) taking into account the type:
###

#Selecting the variable type from the file information:
selection = []; #Defining the variable
for i=1 : rows(channels)
#Choosing which adaptation is necessary from the channel type characters:
if strcmp(chantype(i,1:4),"Div.")==1
 selection = [selection,1]; %Number 1: Resistive divider
endif;
if strcmp(chantype(i,1:5),"Sonda")==1
 selection = [selection,2]; %Number 2: Tektronik voltage probe
endif;
if strcmp(chantype(i,1:3),"Rog")==1
 selection = [selection,3]; %Number 3: Rogowsky coil probe
endif;
if ( length(chantype(i,:))>8 )
    if strcmp(chantype(i,1:10),"Photodiode")==1
     selection = [selection,4]; %Number 4: Photodiode signal
    endif;
endif;
endfor;

# Data transformation of the channels from previous selection:
for i=1 : columns(selection)
 switch(selection(i))
    case 1 %My resistive divider data
      volt(:,2) = vol(:,i) ./ 0.00032; %Wire volts
      volt(:,1) = t(:,i).*1e6; %time in us
      volt(:,2) = supsmu(volt(:,1),volt(:,2),"span",0.01); %Smoothing the voltage data.
	#It allows for a better treatment of the data, without losing too much information.
      disp("Transforming voltage data from resistive divider");
    case 2 %The tektronik probe
      volt(:,2) = vol(:,i); %Wire volts. No transformation, but not clear if correct.
      volt(:,1) = t(:,i).*1e6; %time in us
      volt(:,2) = supsmu(volt(:,1),volt(:,2),"span",0.01); %Smoothing
      disp("Transforming voltage data from Tektronik probe");
    case 3 %Rogowsky coil
      dint(:,2) = (vol(:,i).*63.095)./5.85e-9; 
        %Following calibration page 33 of "Logbook of diverse projects 01"
      dint(:,1) = t(:,i).*1e6; %time in us
      dint(:,2) = supsmu(dint(:,1),dint(:,2),"span",0.01); %Smoothing
      disp("Transforming Rogowsky signal");
    case 4 %Photodiode signal (Normalized to 1)
      phot(:,2) = vol(:,i);
      phot(:,2) = phot(:,2)./max(phot(:,2)); 
      phot(:,1) = t(:,i).*1e6; %time in us
      phot(:,2) = supsmu(phot(:,1),phot(:,2),"span",0.02); %Smoothing
      disp("Transforming photodiode signal");
 endswitch
endfor;

#Checking and warning for missed data on the shot:
if (exist("volt","var")==0)
  warning("Script alex-kinetic: There is no voltage signal in shot %s",shot);
elseif (exist("dint","var")==0)
  warning("Script alex-kinetic: There is no intensity signal in shot %s",shot);
elseif (exist("phot","var")==0)
  warning("Script alex-kinetic: There is no photodiode signal in shot %s",shot);
endif;


###
#Removing values on the baseline for integration or further operations of the data:
###

disp("Removing the baselines form the scope data.");
#Finding the baseline values (if the variables exist):
#Putting the voltage baseline to zero:
if (exist("volt","var")==1)
  [basev,mediav,puntosv] = baseline(volt(1:250,2)); %Voltage signal
  if (puntosv<=10)
    warning("Script alex-kinetic: only %u points in voltage baseline",puntosv);
  endif;
  volt(:,2) = volt(:,2) - mediav;
endif;
#The Rogowsky signal baseline to zero:
if (exist("dint","var")==1)
  [basei,mediai,puntosi] = baseline(dint(1:250,2)); %Intensity signal
  if (puntosi<=10)
    warning("Script alex-kinetic: only %u points in intensity derivative baseline",puntosi);
  endif;
  dint(:,2) = dint(:,2) - mediai;
endif;
#Photodiode signal to zero:
if (exist("phot","var")==1)
  [basep,mediap,puntosp] = baseline(phot(1:200,2)); %Photodiode signal
  if (puntosp<=10)
    warning("Script alex-kinetic: only %u points in photodiode baseline",puntosp);
  endif;
  phot(:,2) = phot(:,2) - mediap;
endif;


###
# Integrating the intensity derivative:
###

if (exist("dint","var")==1)
  int = - cumtrapz(dint(:,1)./1e6,dint(:,2)); %The minus because of the signal in the voltage and the derivative.
  #The division of 1e6 is because we need to use ALL unit in I.S!!!!!
endif;


####
## Finding the electrical power and energy delivered by time and the total: 
####

if ( (exist("int","var")==1) && (exist("volt","var")==1) ) 
  elpowr = volt(:,2).*int; #Electrical power (Watts)
  #Electrical energy delivered:
  elenergia(:,1) = dint(:,1); %Time (microseconds)
  elenergia(:,2) = cumtrapz(elpowr).*(abs(dint(1,1)-dint(2,1)).*1e-6); %Energy (Joules)
endif;







#########################################################################################################
#
#	STREAK IMAGE INTERPRETATION PART.
#
#########################################################################################################


####
## Opening the picture and choosing some parameters.
####

filename = horzcat(shot,".dat"); %Compatibility with previous code and data treatment.

disp("Opening streak image.");

[file, msg] = fopen(filename, "r");
if (file == -1) 
   error ("Unable to open file name: %s, %s",filename, msg); 
endif; 

#Standard deviation parameter (Value to make the binarization)
std_dev_par = 35; %This seems to work, change if results are very bad.

#Charge the image as a matrix: (IT MUST BE ALREADY FORMATTED AS SUCH)
T = dlmread(file);
fclose(file);

#Rounding values for the text final file:
redond = [0 0 0 0 0 0 0 0];

#The binarized matrix. empty with zeroes now:
Out = zeros(rows(T),columns(T));

#The output vector. As before, it starts with zeroes inside:
move = zeros(rows(T),8);

#Maximum of the matrix/2. For showing the results in images and finding edges.
valout = max(max(T))/2;


###
# Making the binarized matrix.
###
#   Loop to fill the OUT matrix with the positions and binarized image: 
#  It places "valout" on the output matrix. So this matrix has only two values: 0
# and "valout". Then it is easy to look for the external edges.
#  It is done on a double loop because with vector logics and so on,
# the program gives not right results.

disp("Making the binarized matrix.");

for i = 1 : rows(T) #For every row of the picture...
    row = T(i,:); # put the row on a variable.
    if ( std(row) > std_dev_par ) #When there is data on the row, no just noise ->
    % Signaled by a standard deviation higher than "std_dev_par".
      halfrow = max(row)/2; # Take half of the maximum value for this row and use it as mark for margins
	for j = 1 : columns(T)
           %  In every column of the OUT row change the value to "valout" if 
           % the column value is higher than halfrow.
	   if ( row(j) > halfrow )
	     Out(i,j) = valout;
	   endif;
	endfor;
    endif;
endfor;


###
# Find the positions of the edges in the binary matrix
###
#Loop for finding the real positions pixels table:

disp("Finding the edges in the matrix.");

j = 0; %Initialize the "j" variable.
for i = 1 : rows(Out) %For every row of the picture...
   pos = find(Out(i,:)==valout); %Find the positions of "valout" in the output matrix and
   # store them in an index vector(POS)
   if (columns(pos) > 1) %When there are really values
     j = j+1;
     move(j,1) = i; %First column is time (In pixels units)
     move(j,2) = pos(1); %Second column is space on "first" side
     move(j,3) = pos(end); %Third column is space on the "last" side
     if (j==1) %Finding the center position for the first column:
       zenter = abs(move(1,2)-move(1,3))/0.5 + min(move(1,2),move(1,3));
     endif;
   endif;
endfor;


###
# Streak data conversion in radius and velocity over time.
###

disp("Streak data treatment and transformation into radius and velocity over time.");

#Remove the zeros on the matrix results:
 move = move( (move(:,2)!=0),:); %Take out the positions were the second column has zero value.
#(Made using logical indexing)

#Calibration of pixels in time and space and centering:
#CAREFUL: THIS CALIBRATION IS VALID FROM ALEX086 TO ALEX099 ONLY !!!!!
move(:,4) = move(:,1) .* sweep; #Passing pixels to nanoseconds.
move(:,5) = (move(:,2)-zenter) .* 0313; #Passing pixels to micrometers and centering.
move(:,6) = (move(:,3)-zenter) .* 0313; #Passing pixels to micrometers and centering.

#smoothing radius data with the function "supsmu", check function help to see how it works:
move(:,5) = supsmu(move(:,4),move(:,5),"span",0.01);
move(:,6) = supsmu(move(:,4),move(:,6),"span",0.01);


###
# Data transforming into radius expansion and its velocity
###

#Putting in zero the radius displacement (different for each side)
#This configuration depends on the streak image orientation. 
#Is now made for "up" left side and "down" right side.
move(:,5) = -(move(:,5) - max(move(:,5))); 
move(:,6) = move(:,6) - min(move(:,6)); 

#Name the displacement vector in a more convinient way:
rup = ( (move(:,5) .* 1e-6) + r0) ./ r0; %transforming micrometers into dimensionless space.
rdown = ( (move(:,6) .* 1e-6) + r0) ./ r0; %transforming micrometers into dimensionless space.
t = move(:,4) .* 1e-9; %transforming nanoseconds in seconds.

#Deriving to obtain velocity (in micrometers/nanosecond):
  h = abs(move(10,4)-move(11,4));
dev = deri(move(:,5),h);
dev2 = deri(move(:,6),h);

move(1:columns(dev),7) = 1000 .* dev; %Transforming it in m/s from um/ns and placing it in the "move" matrix.
move(1:columns(dev),8) = 1000 .* dev2; %Transforming it in m/s from um/ns and placing it in the "move" matrix.
#Renaming velocities:
vup = move(:,7); 
vdown = move(:,8);


###
# Finding the kinetic energy:
###

kinenergiaup(:,1) = move(:,4).*1e-3; %Time in microseconds
kinenergiadown(:,1) = kinenergiaup(:,1);
interenergiaup(:,1) = kinenergiaup(:,1);
internergiadown(:,1) = kinenergiaup(:,1); 

#Polynomia from velocities upper side of the expansion:
[vup_pol, rcorup, polestr] = polyadj(rup.*r0,vup,0.99); %Adjusting to polynomium the velocity data with the correct dimensions.
vsqup = conv(vup_pol,vup_pol); %Making the product of these polynomia.
kint_up = polyint(vsqup); %Integrating v^2.
kinup(:,2) = polyval(kint_up,rup.*r0); %Kinetic energy...
#Adjusting the dimensions:
kinenergiaup(:,2) = rho .* r0 .* lengthwire .* pi .* kinup(:,2);

#Polynomia from velocities down side of the expansion:
[vdown_pol, rcordown, polestr] = polyadj(rdown.*r0,vdown,0.99); %Adjusting to polynomium the velocity data with the correct dimensions.
vsqdown = conv(vdown_pol,vdown_pol); %Making the product of these polynomia.
kint_down = polyint(vsqdown); %Integrating v^2.
kindown(:,2) = polyval(kint_down,rdown.*r0); %Kinetic energy...
#Adjusting the dimensions:
kinenergiadown(:,2) = rho .* r0 .* lengthwire .* pi .* kindown(:,2);


###
# Finding the internal plasma energy:
###

#Simplified EOS (Eq. 30  paper: 
#"Adiabatic plasma heating and fusion-energy production by a compressible fast liner", Nucl. Fusion Vol. 19, pp. 155-177 (1979) )
e1 = ( pi .* r0.**2 .* lengthwire .* Bi) ./  (B-1); 
intenergiaup(:,2) = e1 .*  rup.**2 .* ( (1./B).*(1./rup).**(2.*B) -(1./B) -(1./rup).**2 +1  );
intenergiadown(:,2) = e1 .*  rdown.**2 .* ( (1./B).*(1./rdown).**(2.*B) -(1./B) -(1./rdown).**2 +1  );



#########################################################################################################
#
#	PLOTING AND SAVING INTERNAL, KINETIC AND ELECTRICAL ENERGY PARTS AND PARAMETERS.
#
#########################################################################################################

####
# CAREFUL!!!!! THE ADJUSTMENT FOR THE TIME ORIGIN IS VERY "CHAPUZAS": Just based on pure experience
####


#####
## General figure features.
#####

h = figure(1); %Header to the figure
set (h,'paperunits','centimeters'); %Units to be used.
set (h,'paperorientation','landscape'); %Orientation of the file/paper obtained. NOT VISIBLE ON SCREEN!!!!
set (h,'papersize',[10, 7]); %Paper size in the units defined.
set (h,'position', [0.15, 0.15, 1, 1] .* [10, 7, 10, 7]); %Set bottom-left position, width and height of figure on the paper/file output.
set (h,'defaultaxesposition', [0.15, 0.18, 0.8, 0.7]); %Set the axis position lower left corner and their relative size on the plot.
set (h,'defaultaxesfontsize', 11);
set (h,'visible','off'); %Do not show the plot. Plot to be printed.

h2 = figure(2); %Header to the figure
set (h2,'paperunits','centimeters'); %Units to be used.
set (h2,'paperorientation','landscape'); %Orientation of the file/paper obtained. NOT VISIBLE ON SCREEN!!!!
set (h2,'papersize',[10, 7]); %Paper size in the units defined.
set (h2,'position', [0.15, 0.15, 1, 1] .* [10, 7, 10, 7]); %Set bottom-left position, width and height of figure on the paper/file output.
set (h2,'defaultaxesposition', [0.15, 0.18, 0.8, 0.7]); %Set the axis position lower left corner and their relative size on the plot.
set (h2,'defaultaxesfontsize', 11);
set (h2,'visible','off'); %Do not show the plot. Plot to be printed.

#####
## Plotting the data and placing info on the plot.
#####

disp("Creating and saving data files.");

#Upper radius:
figure(1); %Setting the window plot to the plot window 1.
plot(elenergia(:,1),elenergia(:,2),"*r","markersize",3, %Electrical energy data
kinenergiaup(:,1)-6.4,kinenergiaup(:,2),".-b","markersize",7, %kinetic energy data (time adjustment)
kinenergiaup(:,1)-6.77,intenergiaup(:,2),".-g"); %Internal energy data (time adjustment)

text(5,50,'.- kinetic','color','blue'); %Labels on the plot. First options, coordinates on the plot dimensions
text(5,30,'* electrical','color','red'); 
text(5,10,'.- EOS internal','color','green');
xlabel("time (microseconds)");
title(horzcat(shot," energy data. Upper radius."));
axis([-2 10 0 200]); %Axis limits. First x axis, then y axis (In the units of the vector files).
%close(); %To close the figure and avoid interference

%sleep(1);

#down radius:
figure(2); %Setting the window plot to the plot window 1.
plot(elenergia(:,1),elenergia(:,2),"*r","markersize",3, %Electrical energy data
kinenergiadown(:,1)-6.4,kinenergiadown(:,2),".-b","markersize",7, %kinetic energy data
kinenergiadown(:,1)-6.77,intenergiadown(:,2),".-g"); %Internal energy data

text(5,50,'.- kinetic','color','blue'); %Labels on the plot. First options, coordinates on the plot dimensions
text(5,30,'* electrical','color','red'); 
text(5,10,'.- EOS internal','color','green');
xlabel("time (microseconds)");
tit = horzcat(shot," energy data. Down radius.");
title(tit);
axis([-2 10 0 200]); %Axis limits. First x axis, then y axis (In the units of the vector files).

#Printing the plots:
print(h,horzcat(shot,"-up.jpg"), "-r300"); %Print a jpg file with 300ppi resolution.
print(h2,horzcat(shot,"-down.jpg"), "-r300"); %Print a jpg file with 300ppi resolution.

###
#Saving the data in files:
###
#Voltage (if there is)
if (exist("volt","var"))
  #Output file name:
  name = horzcat(shot,"_voltage.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"time(micros)  voltage(V)");
  redond = [2 0]; %Saved precision 
  display_rounded_matrix(volt, redond, output); %This function is not made by my.
  fclose(output); %Closing the file.
  disp("Voltage saved.");
endif;
#Intensity (if there is)
if (exist("dint","var"))
  #Output file name:
  name = horzcat(shot,"_intensity.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"time(micros)  intensity(A)");
  redond = [2 0]; %Saved precision 
  intensity = [dint(:,1),int];
  display_rounded_matrix(intensity, redond, output); 
  fclose(output); %Closign the file.
  disp("Intensity saved.");
endif;
#Photodiode signal (if there is)
if (exist("phot","var"))
  #Output file name:
  name = horzcat(shot,"_photodiode.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"time(micros)  luminosity(A.U.)");
  redond = [2 0]; %Saved precision 
  display_rounded_matrix(phot, redond, output); %This function is not made by my.
  fclose(output); %Closing the file.
  disp("Relative luminosity saved.");
endif;
#Electrical energy (if there is)
if (exist("elenergia","var"))
  #Output file name:
  name = horzcat(shot,"_electrical_energy.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"time(micros)  energy(J)");
  redond = [2 3]; %Saved precision 
  display_rounded_matrix(elenergia, redond, output); %This function is not made by my.
  fclose(output);%Closing file
  disp("Electrical energy saved.");
endif;
#Kinetic energy of the upper side (if there is)
if (exist("kinenergiaup","var"))
  #Output file name:
  name = horzcat(shot,"_KinEner_up_energy.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"time(micros)  energy(J)");
  redond = [2 3]; %Saved precision 
  display_rounded_matrix(kinenergiaup, redond, output); %This function is not made by my.
  fclose(output); %Closing file
  disp("Kinetic upper energy saved.");
endif;
#Kinetic energy of the down side (if there is)
if (exist("kinenergiadown","var"))
  #Output file name:
  name = horzcat(shot,"_KinEner_down_energy.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"time(micros)  energy(J)");
  redond = [2 3]; %Saved precision 
  display_rounded_matrix(kinenergiadown, redond, output); %This function is not made by my.
  fclose(output);
  disp("Kinetic down energy saved.");
endif;
#Internal energy of the upper side (if there is)
if (exist("intenergiaup","var"))
  #Output file name:
  name = horzcat(shot,"_IntEner_up_energy.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"time(micros)  energy(J)");
  redond = [2 3]; %Saved precision 
  display_rounded_matrix(intenergiaup, redond, output); %This function is not made by my.
  fclose(output);
  disp("Internal up energy saved.");
endif;
#Internal energy of the down side (if there is)
if (exist("intenergiadown","var"))
  #Output file name:
  name = horzcat(shot,"_IntEner_down_energy.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"time(micros)  energy(J)");
  redond = [2 3]; %Saved precision 
  display_rounded_matrix(intenergiadown, redond, output); %This function is not made by my.
  fclose(output);
  disp("Internal up energy saved.");
endif;

more on; #Revert more normal comments behaviour.

###
# Total processing time
###
timing = toc;
disp("Script alex execution time:")
disp(timing)
disp(" seconds")

#That's...that's all, folks! 
