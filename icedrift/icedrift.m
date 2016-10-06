%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script asks questions for icedrift calculation                   %%
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



% This is a GUI for asking questions for icedrift calculation
prompt = {'LatDegrees','LongDegrees'};
userinput = inputdlg(prompt,'Start Coordinates',1,{'LatDegrees','LongDegrees'});
startcoordinate(1) = str2double(cell2mat(userinput(1)));
startcoordinate(2) = str2double(cell2mat(userinput(2)));

prompt = {'LatDegrees','LongDegrees'};
userinput = inputdlg(prompt,'Final Coordinates',1,{'LatDegrees','LongDegrees'});
stopcoordinate(1)= str2double(cell2mat(userinput(1)));
stopcoordinate(2) = str2double(cell2mat(userinput(2)));

icedriftdistance = icedriftcalc(startcoordinate(1),startcoordinate(2),stopcoordinate(1),stopcoordinate(2));
msgbox(strcat('The ice drifted =>',num2str(icedriftdistance),'<= Nautical Miles'));

clear prompt userinput