function u=physical_units_struct(System_of_units,base_time_units)
% Create a struct with physical unit conversion factors.
%
% function u=physical_units_struct(System_of_units,base_time_units)
%     
% 
% Inputs:
% --System_of_units
%
% if System_of_units=='US' or System_of_units==3
%    basic assumed units are American Engineering:
%    LENGTH=FT, TIME=SEC, MASS=SLUG TEMPERATURE=RAN FORCE=LB
%
% elseif System_of_units=='CGS' or System_of_units==2, 
%    basic assumed units are Centimeter,Gram,Second:
%    LENGTH=CM, TIME=SEC, MASS=GM TEMPERATURE=K FORCE=DYNE
%
% elseif System_of_units=='IMPERIAL' or System_of_units==4, 
%    basic assumed units are Imperial:
%    LENGTH=FT, TIME=SEC, MASS=SLUG TEMPERATURE=RAN FORCE=LB
%
% otherwise, 
%    basic assumed units are SI (default):
%    LENGTH=M , TIME=SEC, MASS=KG   TEMPERATURE=K   FORCE=N
%
% --base_time_units defaults to SEC
%
% Outputs:
%    u = struct with fields      
%         SEC,MIN,HR,DAY,WK,MONTH,YR,FT,SLUG,RAN,...
%         M,KG,K,CM,GM,DEG,REV,...
%         IN,NMI,YD,MM,MILE,G,LB,NT,OZ,LBM,KT,...
%         MPH,PSI,PA,BAR,ATM,TORR,mmHG,J,...
%         CAL,MEV,ERG,BTU,W,MW,HP,A,V,...
%         HZ,RPS,RPM,L,GAL,GPM
%
% SEC = second
%
% % time conversions
% MIN = Minute
% HR = hour
% DAY = Day
% WK = Week
% MONTH = Month
% YR = Year
%
%
% % prefixes
% NANO=10^(-9)
% MICRO=10^(-6);
% MILLI=10^(-3);
% KILO=10^3;
% MEGA=10^6;
% GIGA=10^9;
% TERA=10^12;
%
% %  angles
% RAD=Radian (angle)
% DEG=Degree (angle)
% REV=Measure of angles in terms of revolutions
%
% % Length
% IN=Inch
% NMI=Nautical mile
% CM=Centimeter
% YD=Yard
% MM=Millimeter
% MILE=Mile
% G=Gravity acceleration
%
% % Temperature
% K= Degrees Kelvin
% RAN= Degrees Rankine
%
% % Force
% LB=Pound
% NT=newton
% OZ=ounce
% LBM= Pound of mass
%
% % Speed
% KT=knot
% MPH=mile per hour
%
% % Pressure
% PSI=Pound per square inch
% PA= Pascal
% BAR=Bar
% ATM=Atmosphere
% TORR=torr
% mmHG=mm Hg
% BA=CGS unit of pressure, dyne/cm^2
%
% % Work and power
% J= Joule
% CAL = Calorie
% MEV = Mega electron volt
% ERG = erg
% BTU = BTU
% W=watt
% MW = Megawatt
% HP=Horsepower
%
% % Electric quantities
% COUL=Coulomb (charge);
% A = ampere
% V = volt
%
% HZ = Frequency (Hertz)
% RPS = Revolutions per second
% RPM =  Revolutions per minute
%
% L = Liter
% GAL = Gallon
% GPM = Gallons per minute 
% 
%
% Examples: 
%     pu= physical_units_struct;
%     E= 210000*pu.MEGA*pu.PA;% Young's modulus in MPa
%
%
% See also: physical_units

% Defaults:
if (~exist('System_of_units'))
    System_of_units='SI';%
end
if (~exist('base_time_units'))
    base_time_units='SEC';%
end


[uSEC,uMIN,uHR,uDAY,uWK,uMONTH,uYR,uFT,uSLUG,uRAN,...
    uM,uKG,uK,uCM,uGM,uDEG,uRAD,uREV,...
    uIN,uNMI,uYD,uMM,uMILE,uG,uLBF,uNT,uOZ,uLBM,uKT,...
    uMPH,uPSI,uPA,uBAR,uATM,uTORR,ummHG,uJ,...
    uCAL,uMEV,uERG,uBTU,uW,uMW,uHP,uA,uV,...
    uHZ,uRPS,uRPM,uL,uGAL,uGPM,...
    uNANO,uMICRO,uKILO,uMILLI,uMEGA,uGIGA,uTERA,uBBL,uBA,uF]=physical_units(System_of_units,base_time_units);
 
% Check the number of output arguments and report error if different from
% that expected.
if (nargout >1)
    error('Mismatched number of output arguments');
end
for j={'uSEC','uMIN','uHR','uDAY','uWK','uMONTH','uYR','uFT','uSLUG','uRAN',...
    'uM','uKG','uK','uCM','uGM','uDEG','uRAD','uREV',...
    'uIN','uNMI','uYD','uMM','uMILE','uG','uLBF','uNT','uOZ','uLBM','uKT',...
    'uMPH','uPSI','uPA','uBAR','uATM','uTORR','ummHG','uJ',...
    'uCAL','uMEV','uERG','uBTU','uW','uMW','uHP','uA','uV',...
    'uHZ','uRPS','uRPM','uL','uGAL','uGPM',...
    'uNANO','uMICRO','uKILO','uMILLI','uMEGA','uGIGA','uTERA', 'uBBL','uBA','uF'}
    % .(j{1})=eval (j{1});
    n=j{1};
    u.(n(2:end))=eval (j{1});
end

end



function [uSEC,uMIN,uHR,uDAY,uWK,uMONTH,uYR,uFT,uSLUG,uRAN,...
    uM,uKG,uK,uCM,uGM,uDEG,uRAD,uREV,...
    uIN,uNMI,uYD,uMM,uMILE,uG,uLBF,uNT,uOZ,uLBM,uKT,...
    uMPH,uPSI,uPA,uBAR,uATM,uTORR,ummHG,uJ,...
    uCAL,uMEV,uERG,uBTU,uW,uMW,uHP,uA,uV,...
    uHZ,uRPS,uRPM,uL,uGAL,uGPM,...
    uNANO,uMICRO,uKILO,uMILLI,uMEGA,uGIGA,uTERA,uBBL,uBA,uF]=physical_units(System_of_units,base_time_units)
% Create physical unit conversion factors.
%
% [uSEC,uMIN,uHR,uDAY,uWK,uMONTH,uYR,uFT,uSLUG,uRAN,...
%     uM,uKG,uK,uCM,uGM,uDEG,uRAD,uREV,...
%     uIN,uNMI,uYD,uMM,uMILE,uG,uLBF,uNT,uOZ,uLBM,uKT,...
%     uMPH,uPSI,uPA,uBAR,uATM,uTORR,ummHG,uJ,...
%     uCAL,uMEV,uERG,uBTU,uW,uMW,uHP,uA,uV,...
%     uHZ,uRPS,uRPM,uL,uGAL,uGPM,...
%     uNANO,uMICRO,uKILO,uMILLI,uMEGA,uGIGA,uTERA, uBBL,uBA,uF]=physical_units(System_of_units,base_time_units)
%     
%
% Inputs:
% --System_of_units defaults to SI
% if System_of_units=='US', basic assumed units are American Engineering:
%    LENGTH=uFT, TIME=uSEC, MASS=uSLUG TEMP=uRAN FORCE=uLBF
% if System_of_units=='CGS', basic assumed units are Centimeter,Gram,Second:
%    LENGTH=uCM, TIME=uSEC, MASS=uGM TEMP=uK FORCE=uDYNE
% if System_of_units=='IMPERIAL', basic assumed units are Imperial:
%    LENGTH=uFT, TIME=uSEC, MASS=uSLUG TEMP=uRAN FORCE=uLBF
% otherwise, basic assumed units are SI:
%    LENGTH=uM , TIME=uSEC, MASS=uKG   TEMP=uK   FORCE=uNT
% --base_time_units defaults to SEC
%
% Outputs:
% 
% uSEC = second
% % time conversions
% uMIN = Minute
% uHR = hour
% uDAY = Day
% uWK = Week
% uMONTH = Month
% uYR = Year
% uCOUL=Coulomb (charge);
% uRAD=Radian (angle)
% % prefixes
% uNANO=10^(-9)
% uMICRO=10^(-6);
% uMILLI=10^(-3);
% uKILO=10^3;
% uMEGA=10^6;
% uGIGA=10^9;
% uTERA=10^12;
% uDEG=Degree (angle)
% uREV=Measure of angles in terms of revolutions
% uIN=Inch
% uNMI=Nautical mile
% uCM=Centimeter
% uYD=Yard
% uMM=Millimeter
% uMILE=Mile
% uG=Gravity acceleration
% uLBF=Pound (Force)
% uNT=newton
% uOZ=ounce
% uLBM= Pound of mass
% uKT=knot
% uMPH=mile per hour
% uPSI=Pound per square inch
% uPA= Pascal
% uBAR=Bar
% uATM=Atmosphere
% uTORR=torr
% ummHG=mm Hg
% uJ= Joule
% uCAL = Calorie
% uMEV = Mega electron volt
% uERG = erg
% uBTU = BTU
% uW=watt
% uMW = Megawatt
% uHP=Horsepower
% uA = ampere
% uV = volt
% uHZ = Frequency (Hertz)
% uRPS = Revolutions per second
% uRPM =  Revolutions per minute
% uL = Liter
% uGAL = Gallon
% uGPM = Gallons per minute 
% 
% Based on UNITS/GL_UNITS by Glenn Gomes, Version 2011.
%
% Examples: 
%    [uSEC,uMIN,uHR,uDAY,uWK,uMONTH,uYR,uFT,uSLUG,uRAN,...
%     uM,uKG,uK,uCM,uGM,uDEG,uRAD,uREV,...
%     uIN,uNMI,uYD,uMM,uMILE,uG,uLBF,uNT,uOZ,uLBM,uKT,...
%     uMPH,uPSI,uPA,uBAR,uATM,uTORR,ummHG,uJ,...
%     uCAL,uMEV,uERG,uBTU,uW,uMW,uHP,uA,uV,...
%     uHZ,uRPS,uRPM,uL,uGAL,uGPM,...
%     uMICRO,uKILO,uMILLI,uMEGA,uGIGA,uTERA,uBBL,uBA]=physical_units;
%     E= 210000*uMEGA*uPA;% Young's modulus in MPa

if (~exist('System_of_units')) %#ok<*EXIST>
    System_of_units='SI';
end
if (~exist('base_time_units'))
    base_time_units='SEC';
end

if (strcmp(upper(base_time_units),'SEC')) %#ok<*STCI>
    uSEC = 1.0;
elseif (strcmp(upper(base_time_units),'MIN'))
    uSEC = 1/60;
elseif (strcmp(upper(base_time_units),'HR'))
    uSEC = 1/(60*60);
elseif (strcmp(upper(base_time_units),'DY'))
    uSEC = 1/(60*60*24);
elseif (strcmp(upper(base_time_units),'YR'))
    uSEC = 1/(60*60*24*365);
elseif (strcmp(upper(base_time_units),'WK'))
    uSEC = 1/(60*60*24*7);
else
    disp('The only valid entries for the time base units are: SEC|MIN|HR|DY|YR')
    return
end

% time conversions
uMIN = uSEC*60;
uHR = uSEC*60*60;
uDAY = uSEC*60*60*24;
uWK = uSEC*60*60*24*7;
uMONTH = uSEC*60*60*24*365/12;
uYR = uSEC*60*60*24*365;

% base units for charge, time, angle (all systems)
uCOUL=1.0;uRAD=1.0;

% base units for length, mass, temperature (different for each system)

if (strcmp(upper(System_of_units),'US'))
    strflag='US';
    uFT=1.0; uSLUG=1.0; uRAN=1.0;                     % US base units
    uM=uFT*3*(1/0.9144);uKG= uSLUG/14.593903591998501; uK=1.8*uRAN; % SI base units
    uCM=uM/100.0;uGM= uKG/1000.0;                     % CGS base units
elseif (strcmp(upper(System_of_units),'IMPERIAL'))
    strflag='IMPERIAL';
    uFT=1.0; uSLUG=1.0; uRAN=1.0;                     % US base units
    uM=uFT*3*(1/0.9144);uKG= uSLUG/14.593903591998501; uK=1.8*uRAN; % SI base units
    uCM=uM/100.0;uGM= uKG/1000.0;                     % CGS base units
elseif (strcmp(upper(System_of_units),'CGS'))
    strflag='CGS'; %#ok<*NASGU>
    uCM=1.0; uGM=1.0; uK=1.0;						  % CGS base units
    uM=100.0*uCM;uKG=1000.0*uGM; 					  % SI base units
    uFT=uM/(3*(1/0.9144));uSLUG=uKG*14.593903591998501;uRAN=uK/1.8; % US base units
else %SI is the default
    strflag='SI';
    uM=1.0; uKG=1.0; uK=1.0;						  % SI base units
    uCM=uM/100.0; uGM=uKG/1000.0;               % CGS base units
    uFT=uM/(3*(1/0.9144));uSLUG=uKG*14.593903591998501;uRAN=uK/1.8; % US base units
end

% prefixes
uNANO=10^(-9);
uMICRO=10^(-6);
uMILLI=10^(-3);
uKILO=10^3;
uMEGA=10^6;
uGIGA=10^9;
uTERA=10^12;

% angle conversions (degrees of angle)
uDEG=pi*uRAD/180.0;
uREV=2.0*pi;% Measure of angles in terms of revolutions

% length conversions
uIN=uFT/12;
uNMI=1852*uM;
uCM=uM/100.0;
uYD=3.0*uFT;
uMM=uM/1000.0;
uMILE=5280*uFT;

% acceleration conversions
uG=32.17405*uFT/uSEC^2;

% force conversions
uLBF=uSLUG*uFT/uSEC^2;
uNT=uKG*uM/(uSEC^2); %newton
uOZ=uLBF/16.0;

% mass conversions
uLBM=uLBF/(uG);

% speed conversions
uKT=uNMI/uHR;
uMPH=uMILE/uHR;

% pressures
uPSI=uLBF/(uIN^2);
uPA= uNT/(uM^2);
uBAR=10^5*uPA;
uATM=14.6959*uPSI;
uTORR=133.322*uPA;
ummHG=uTORR;
uBA=0.1*uPA;%cgs unit of pressure, Barye    

% work
uJ= uNT*uM;
uCAL = uJ * 4.1868;
uMEV = uJ / 6.242e12;
uERG = uJ * 1E-7;
uBTU = 1054*uJ;

% power
uW=uJ/uSEC;
uMW = 1e6*uW;
uHP=745.7*uW;

% electricity
uA = uCOUL/uSEC;
uV = uJ/uCOUL;

% frequency
uHZ = 1.0/uSEC;
uRPS = uHZ;% Revolutions per second
uRPM = 1/uMIN;% Revolutions per minute

% volume
uL = 1000.0*uCM^3;

if (strcmp(upper(System_of_units),'IMPERIAL')) %#ok<*STCUL>
    uGAL = 4.54609*uL;% Imperial gallon
else
    uGAL = 3.785411784*uL;% US gallon
end
uBBL = 159*uL;% standard barrel for measuring volume of oil http://en.wikipedia.org/wiki/Barrel_(unit)

% flow rate
uGPM = uGAL/uMIN;% gallons per minute

% Temperature, degrees Fahrenheit
uF=uRAN;

% Check the number of output arguments and report error if different from that expected.
if (nargout ~=62)
    error('Mismatched number of output arguments');
end
end