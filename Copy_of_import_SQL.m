%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This scripts runs the import of SQL data from Frej or Oden           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script allows the import of SQL data from Frej or Oden into a
% structure for later processing
% 
%
% Input data:
%   A GUI requests the .csv files from the SQL database export
%   
% Output data:
%   Excelfile
%  
%
%    Copyright:     NTNU
%    Project:	    SAmCoT, AMOS
%    Author:        Hans-Martin Heyn
%    Date created:  2016-08-21  Hans-Martin Heyn (NTNU)
%    

%---------------------------------------------------------------------%
function import_SQL()

fprintf('o---------------------------------------------o\n')
fprintf('|\t The SQL import tool for Michelle V1.0   \t|\n')
fprintf('o---------------------------------------------o\n\n');

%% User input 
% Request file for shipdata
% [CSVPressureName,CSVPressureDir] = uigetfile('.csv','Select the pressure file','*.csv');
% [CSVRotoName,CSVRotoDir] = uigetfile('.csv','Select the rototronic file',strcat(CSVPressureDir,'*.csv'));
% [CSVName,CSVDir] = uigetfile('.csv','Select the shipdata file',strcat(CSVPressureDir,'*.csv'));
% [CSVWaterName,CSVWaterDir] = uigetfile('.csv','Select the water file',strcat(CSVPressureDir,'*.csv'));

% Request folder with shipdata
CSVDir = uigetdir(pwd,'Please select the folder with shipdata');

% Create stringvector with dates
datestringvec = [{'2016-08-09.matlab.csv'},{'2016-08-10.matlab.csv'},{'2016-08-11.matlab.csv'},{'2016-08-12.matlab.csv'},{'2016-08-13.matlab.csv'},{'2016-08-14.matlab.csv'},{'2016-08-15.matlab.csv'},...
    {'2016-08-16.matlab.csv'},{'2016-08-17.matlab.csv'},{'2016-08-18.matlab.csv'},{'2016-08-19.matlab.csv'},{'2016-08-20.matlab.csv'},{'2016-08-21.matlab.csv'},{'2016-08-22.matlab.csv'},...
    {'2016-08-23.matlab.csv'},{'2016-08-24.matlab.csv'},{'2016-08-25.matlab.csv'},{'2016-08-26.matlab.csv'},{'2016-08-27.matlab.csv'},...
    {'2016-08-28.matlab.csv'},{'2016-08-29.matlab.csv'},{'2016-08-30.matlab.csv'},{'2016-08-31.matlab.csv'},{'2016-09-01.matlab.csv'},{'2016-09-02.matlab.csv'},{'2016-09-03.matlab.csv'},{'2016-09-04.matlab.csv'},{'2016-09-05.matlab.csv'},{'2016-09-06.matlab.csv'},{'2016-09-07.matlab.csv'},{'2016-09-08.matlab.csv'},...
    {'2016-09-09.matlab.csv'},{'2016-09-10.matlab.csv'},{'2016-09-11.matlab.csv'},{'2016-09-12.matlab.csv'},{'2016-09-13.matlab.csv'},]; % To be continued

% Where should the excel file go
[ExcelName,ExcelDir] = uiputfile('.xlsx','Please specify where you want to save the excel file for on the spot values',strcat(pwd,'\awesomeshipdata.xlsx'));
[ExcelNameAV,ExcelDirAV] = uiputfile('.xlsx','Please specify where you want to save the excel file for averaged values',strcat(ExcelDir,'\awesomeshipdata_average.xlsx'));

% Request desired time interval
prompt = {'Which time interval in seconds do you want (in 5 second steps please)?'};
userinput = inputdlg(prompt,'Time interval [s]',1,{'300'});
interval = str2double(cell2mat(userinput(1)));

t1 = 0; % Running number for all
t2 = 0;
t3 = 0;
t4 = 0;
t5 = 0;
%% Reading of raw data and converting
for u = 1:1:length(datestringvec)
    display(strcat('Processing ',datestringvec{u},':) :) :) :) :)'))
    
    display('Reading Shipdata')
    try
        [shipdata.timestamp,shipdata.heading,shipdata.COG,shipdata.SOG,shipdata.GPS_lon,shipdata.GPS_lat,shipdata.windDirTrue,shipdata.windSpeedTrue,shipdata.windDirRel,shipdata.windSpeedRel]=Tool_Import_SQL_Oden(strcat(CSVDir,'\ship_',datestringvec{u}),2,inf);
    catch
        display('Ups cannot read ship data and that is serious because I have to abort')
        error(strcat('Sorry without shipdata I am screwed. Remove the ',datestringvec{m},'from the vector of dates'))
    end 
    display('Reading Rotronicdata')
    try
        [roto.timestamp,roto.tempport,roto.tempstbd,roto.humport,roto.humstbd,roto.temperature,roto.humidity] = Tool_Import_SQL_Rototronic(strcat(CSVDir,'\rotronic_',datestringvec{u}),2,inf);
    catch
        roto.timestamp = shipdata.timestamp;
        roto.tempport = zeros(size(shipdata.timestamp));
        roto.tempstbd = zeros(size(shipdata.timestamp));
        roto.humport = zeros(size(shipdata.timestamp));
        roto.humstbd = zeros(size(shipdata.timestamp));
        roto.temperature = zeros(size(shipdata.timestamp));
        roto.humidity = zeros(size(shipdata.timestamp));
        display('Aaaaaaaaah something went wrong with Roto')
    end
    
    display('Reading Waterdata')
    try
        [water.timestamp,water.tempsea,water.templab,water.conductivitylab,water.salinitylab,water.soundvelocitylab] = Tool_Import_SQL_water(strcat(CSVDir,'\water_',datestringvec{u}),2,inf);
    catch
        water.timestamp = shipdata.timestamp;
        water.tempsea = zeros(size(shipdata.timestamp));
        water.templab = zeros(size(shipdata.timestamp));
        water.conductivitylab = zeros(size(shipdata.timestamp));
        water.salinitylab = zeros(size(shipdata.timestamp));
        water.soundvelocitylab = zeros(size(shipdata.timestamp));
        display('Big problems with the water')
    end
    
    display('Reading Pressuredata')
    try
        [pressure.timestamp,pressure.pressure] = Tool_Import_SQL_pressure(strcat(CSVDir,'\pressure_',datestringvec{u}), 2, inf);
    catch
        pressure.timestamp = shipdata.timestamp;
        pressure.pressure = zeros(size(shipdata.timestamp));
        display('Pressure seems to be wrong')
    end
    
    display('Reading Insolationdata')
    try
        [insolation.timestamp,insolation.insolation,insolation.boardTemp,insolation.intVoltage] = Tool_Import_SQL_insolation(strcat(CSVDir,'\insolation_',datestringvec{u}), 2, inf);
    catch
        insolation.timestamp = shipdata.timestamp;
        insolation.insolation = zeros(size(shipdata.timestamp));
        insolation.boardTemp = zeros(size(shipdata.timestamp));
        insolation.intVoltage = zeros(size(shipdata.timestamp));
        display('And here we have a problem with the insolation data')
    end
    
    % Convert date to matlab format
    shipdate = datenum(shipdata.timestamp);
    rotodate = datenum(roto.timestamp);
    waterdate = datenum(water.timestamp);
    pressuredate = datenum(pressure.timestamp);
    insolationdate = datenum(insolation.timestamp);
    datenuminterval = datenum(0,0,0,0,0,interval); % Convert interval to Matlab date format

    %% Initialisation of merging process
    % Init buffer vectors
    water_temp_sea_temp = 0;
    salinity_temp = 0;
    windspeedtrue_temp = 0;
    air_pressure_temp = 0;
    air_temp_temp = 0;
    air_humi_temp = 0;
    insolation_temp = 0;
    water_temp_lab_temp = 0;

    display('Starting of data merging')
    % Resample with desired interval for shipdata

    %% Merging of shipdata
    display('Merging ship data')
    stop = shipdate(1)+datenuminterval-datenum(0,0,0,0,0,10); % First stop
    k = 1; % Counter for averaging
    m = 1; % Counter for output
    for i = 1:1:length(shipdata.timestamp)

        windspeedtrue_temp = windspeedtrue_temp + shipdata.windSpeedTrue(i);
    %     GPS_lon_temp = GPS_lon_temp + shipdata.GPS_lon(i); % Really bad idea!!!
    %     GPS_lat_temp = GPS_lat_temp + shipdata.GPS_lat(i); 
    %     
        k = k + 1;

        if (shipdate(i) >= stop) % Stop

            %display(strcat('Processing timestep ',datestr(stop)));
            % No averaging
            t1 = t1 + 1; % GPS Counter one up for collecting all datapoints
            excelzeit(m,1) = exceltime(shipdata.timestamp(i));
            GPS_lon(m,1) = shipdata.GPS_lon(i);
            GPS_lat(m,1) = shipdata.GPS_lat(i);
            winddirtrue(m,1) = shipdata.windDirTrue(i);
            windspeedtrue(m,1) = shipdata.windSpeedTrue(i);
            speedoverground(m,1) = shipdata.SOG(i);
            
            excelzeit_all(t1,1) = exceltime(shipdata.timestamp(i));
            GPS_lon_all(t1,1) = shipdata.GPS_lon(i);
            GPS_lat_all(t1,1) = shipdata.GPS_lat(i);
            winddirtrue_all(t1,1) = shipdata.windDirTrue(i);
            windspeedtrue_all(t1,1) = shipdata.windSpeedTrue(i);
            speedoverground_all(t1,1) = shipdata.SOG(i);

            % Averaging
            windspeedtrue_av(m,1) = windspeedtrue_temp / k;
            
            windspeedtrue_av_all(t1,1) = windspeedtrue_temp / k;

            
            k = 1; % Reset counter and temps
            windspeedtrue_temp = 0;

            m = m + 1;
            stop = stop + datenuminterval; % Set new interval
        end
    end

    %% Merging of Rototronic data
    display('Merging rototronic data')
    stop =  rotodate(1)+datenuminterval-datenum(0,0,0,0,0,10);% First stop
    k = 1; % Counter for averaging
    m = 1; % Counter for output
    for i = 1:1:length(roto.timestamp)

        air_temp_temp = air_temp_temp + roto.temperature(i);
        air_humi_temp = air_humi_temp + roto.humidity(i);

    %     GPS_lon_temp = GPS_lon_temp + shipdata.GPS_lon(i); % Really bad idea!!!
    %     GPS_lat_temp = GPS_lat_temp + shipdata.GPS_lat(i); 
    %     
        k = k + 1;

        if (rotodate(i) >= stop) % Stop
            
            t2 = t2+1;
            % No averaging
            air_temp(m,1) = roto.temperature(i);
            air_humi(m,1) = roto.humidity(i);
            air_temp_all(t2,1) = roto.temperature(i);
            air_humi_all(t2,1) = roto.humidity(i);
            
            %display(strcat('Processing timestep ',datestr(stop)));
            % Averaging
            air_temp_av(m,1) = air_temp_temp /k;
            air_humi_av(m,1) = air_humi_temp /k;
            air_temp_av_all(t2,1) = air_temp_temp /k;
            air_humi_av_all(t2,1) = air_humi_temp /k;
            
            k = 1; % Reset counter and temps
            air_temp_temp = 0;
            air_humi_temp = 0;

            m = m + 1;
            stop = stop + datenuminterval; % Set new interval
        end
    end

    %% Merging of Water data
    display('Merging water data')
    stop = waterdate(1)+datenuminterval-datenum(0,0,0,0,0,10); % First stop
    k = 1; % Counter for averaging
    m = 1; % Counter for output
    for i = 1:1:length(water.timestamp)

        water_temp_sea_temp = water_temp_sea_temp + water.tempsea(i);
        water_temp_lab_temp = water_temp_lab_temp + water.templab(i);
        salinity_temp = salinity_temp + water.salinitylab(i);

        k = k + 1;

        if (waterdate(i) >= stop) % Stop

            %display(strcat('Processing timestep ',datestr(stop)));
            t3 = t3+1;
            % No averaging
            water_temp_sea(m,1) = water.tempsea(i);
            water_temp_lab(m,1) = water.templab(i);
            salinity(m,1) = water.salinitylab(i);
            water_temp_sea_all(t3,1) = water.tempsea(i);
            water_temp_lab_all(t3,1) = water.templab(i);
            salinity_all(t3,1) = water.salinitylab(i);
            
            % Averaging
            water_temp_sea_av(m,1) = water_temp_sea_temp / k;
            water_temp_lab_av(m,1) = water_temp_lab_temp / k;
            salinity_av(m,1) = salinity_temp / k;
            water_temp_sea_av_all(t3,1) = water_temp_sea_temp / k;
            water_temp_lab_av_all(t3,1) = water_temp_lab_temp / k;
            salinity_av_all(t3,1) = salinity_temp / k;

            k = 1; % Reset counter and temps
            water_temp_sea_temp = 0;
            salinity_temp = 0;

            m = m + 1;
            stop = stop + datenuminterval; % Set new interval
        end
    end

    %% Merging of Pressure data
    display('Merging pressure data')
    stop = pressuredate(1) + datenuminterval-datenum(0,0,0,0,0,10); % First stop
    k = 1; % Counter for averaging
    m = 1; % Counter for output
    for i = 1:1:length(pressure.timestamp)

        air_pressure_temp = air_pressure_temp + pressure.pressure(i);

        k = k + 1;

        if (pressuredate(i) >= stop) % Stop
            
            %display(strcat('Processing timestep ',datestr(stop)));
            t4 = t4 + 1;
            % No averaging
            air_pressure(m,1) = pressure.pressure(i);
            air_pressure_all(t4,1) = pressure.pressure(i);
            
            % Averaging
            air_pressure_av(m,1) = air_pressure_temp / k;
            air_pressure_av_all(t4,1) = air_pressure_temp / k;

            k = 1; % Reset counter and temps
            air_pressure_temp = 0;

            m = m + 1;
            stop = stop + datenuminterval; % Set new interval
        end
    end
        %% Merging of Insolation data
    display('Merging insolation data')
    stop = insolationdate(1) + datenuminterval-datenum(0,0,0,0,0,10); % First stop
    k = 1; % Counter for averaging
    m = 1; % Counter for output
    for i = 1:2:length(insolation.timestamp) % Double steps because of double values in files

        insolation_temp = insolation_temp + insolation.insolation(i);

        k = k + 1;

        if (insolationdate(i) >= stop) % Stop
    
            t5 = t5 + 1;
            % No averaging
            insolationwr(m,1) = insolation.insolation(i);
            insolationwr_calc(m,1) = (insolation.insolation(i) - 0.0097)/(6.418*10^(-4));
            insolationwr_all(t5,1) = insolation.insolation(i);
            insolationwr_calc_all(t5,1) = (insolation.insolation(i) - 0.0097)/(6.418*10^(-4));
            
            % Averaging
            insolationwr_av(m,1) = insolation_temp / k;
            insolationwr_av_calc(m,1) = (insolation.insolation(i) - 0.0097)/(6.418*10^(-4));
            insolationwr_av_all(t5,1) = insolation_temp / k;
            insolationwr_av_calc_all(t5,1) = (insolation.insolation(i) - 0.0097)/(6.418*10^(-4));

            k = 1; % Reset counter and temps
            insolation_temp = 0;

            m = m + 1;
            stop = stop + datenuminterval; % Set new interval
        end
    end

    %% Output to Excel file
    % Write to excel file
    display('Writing to Excel file')
    xlswrite(strcat(ExcelDirAV,ExcelNameAV),[excelzeit,GPS_lon,GPS_lat,speedoverground,water_temp_sea_av,water_temp_lab_av,salinity_av,windspeedtrue_av,winddirtrue,air_temp_av,air_pressure_av,air_humi_av,insolationwr_av,insolationwr_av_calc],datestringvec{u},'A2');
    xlswrite(strcat(ExcelDirAV,ExcelNameAV),{'Time','GPS lon','GPS lat','Speed over Ground not-Averaged','Sea-Water Temperature Averaged','Lab-Water Temperature Averaged','Salinity Lab Averaged','Windspeed true Averaged','Winddirection true Not Averaged','Air Temperature Averaged','Air Pressure Averaged','Air Humidity Averaged','Insolation Averaged (Volt)','Insolation Averaged (\muE/m^2s)'},datestringvec{u},'A1');
    xlswrite(strcat(ExcelDir,ExcelName),[excelzeit,GPS_lon,GPS_lat,speedoverground,water_temp_sea,water_temp_lab,salinity,windspeedtrue,winddirtrue,air_temp,air_pressure,air_humi,insolationwr,insolationwr_calc],datestringvec{u},'A2');
    xlswrite(strcat(ExcelDir,ExcelName),{'Time','GPS lon','GPS lat','Speed over Ground','Sea-Water Temperature','Lab-Water Temperature','Salinity','Windspeed true','Winddirection true','Air Temperature','Air Pressure','Air Humidity','Insolation','Insolation Averaged (\muE/m^2s)'},datestringvec{u},'A1');
    
    xlswrite(strcat(ExcelDirAV,ExcelNameAV),[excelzeit_all,GPS_lon_all,GPS_lat_all,speedoverground_all,water_temp_sea_av_all,water_temp_lab_av_all,salinity_av_all,windspeedtrue_av_all,winddirtrue_all,air_temp_av_all,air_pressure_av_all,air_humi_av_all,insolationwr_av_all,insolationwr_av_calc_all],'all days','A2');
    xlswrite(strcat(ExcelDirAV,ExcelNameAV),{'Time','GPS lon','GPS lat','Speed over Ground not-Averaged','Sea-Water Temperature Averaged','Lab-Water Temperature Averaged','Salinity Lab Averaged','Windspeed true Averaged','Winddirection true Not Averaged','Air Temperature Averaged','Air Pressure Averaged','Air Humidity Averaged','Insolation Averaged (Volt)','Insolation Averaged (\muE/m^2s)'},'all days','A1');
    xlswrite(strcat(ExcelDir,ExcelName),[excelzeit_all,GPS_lon_all,GPS_lat_all,speedoverground_all,water_temp_sea_all,water_temp_lab_all,salinity_all,windspeedtrue_all,winddirtrue_all,air_temp_all,air_pressure_all,air_humi_all,insolationwr_all,insolationwr_calc_all],'all days','A2');
    xlswrite(strcat(ExcelDir,ExcelName),{'Time','GPS lon','GPS lat','Speed over Ground','Sea-Water Temperature','Lab-Water Temperature','Salinity','Windspeed true','Winddirection true','Air Temperature','Air Pressure','Air Humidity','Insolation','Insolation Averaged (\muE/m^2s)'},'all days','A1');
end

display('Done. Have a nice day. :)')
end