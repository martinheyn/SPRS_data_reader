%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function calculates the ice drift between two points             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Input data:
%  
% Output data:
%  estination
%
%    Copyright:     NTNU
%    Project:	    Arctic Ocean 2016
%    Author:        Hans-Martin Heyn
%    Date created:  2016-08-17  Hans-Martin Heyn (NTNU)
%    Change Log: 


% This function calculates the ice drift
function drift = icedriftcalc(startcoordinateLAT,startcoordinateLONG,stopcoordinateLAT,stopcoordinateLONG)
    
startdegreesLAT = startcoordinateLAT;
stopdegreesLAT = stopcoordinateLAT;

startdegreesLONG = startcoordinateLONG;
stopdegreesLONG = stopcoordinateLONG;

[startx,starty,~] = llh2ecef(startdegreesLONG*(pi/180),startdegreesLAT*(pi/180),0); 
[stopx,stopy,~] = llh2ecef(stopdegreesLONG*(pi/180),stopdegreesLAT*(pi/180),0);

deltax = stopx-startx;
deltay = stopy-starty;

drift = sqrt(deltax^2+deltay^2)/1854;

end