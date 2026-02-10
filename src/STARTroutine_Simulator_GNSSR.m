%% ************************************************************************
% START routine for HydroGNSS SAVERS SIMULATOR
%% ************************************************************************
% MODULE NAME:      STARTroutine_Simulator_GNSSR.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Dente, D. Comite, N. Pierdicca, L. Guemrriero
% @copyright:   Tor Vergata University of Rome and Sapienza University of Rome
% Original version:      September 2019 by Laura Dente
% Main updates:          Sep 2019 – Dec 2022 by L. Dente, D. Comite, Aneesha Musunuri
%                        Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel, Tx pattern)
%                        August 2023 by Aneesha Musunuri(Output files)
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module is the main routine of the simulator.
% It opens the main GUI, load all inputs and configuration files,
% calls the routines to read the orbit file, to read the land cover map and to define the
% input bio geo parameters, to run the simulations and to write the outputs in nc files.
% A number of tag can be applied, see lines 28-55.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

 function []=STARTroutine_Simulator_GNSSR(filenamedefault)
%% starting counter

tic
% MP: clear all variables except the file to read
clearvars -except filenamedefault ; 
clc
close all
beep off

%% declaration of global variables
declaration_global_variables

%% Istruction to make exe
rng(10);
%remind to set tag_exe = 1 when not compiling
tag.exe            =1 ; %not zero only to generate exe
%tag.exe=0 for normal simulations and tag.exe=1 when creating the
%executable
%%% mcc -m STARTroutine_Simulator_GNSSR.m
%%to create an executable - ansha 
%% mcc -m ../src/STARTroutine_Simulator_GNSSR.m -a ../src -a ../conf -a ../log -a ../bin/matrici_eq48 -a ../bin/
%%
% if the file (argument) is present NO GUI appears or all windows will
% popup. - added by ansha
if nargin < 1
   tag.GUI            = 0; % if 0 all windows will be opened (added by ansha)
else 
   tag.GUI            = 1; % if 1 the code sets default values from the file presented as the argument
end
tag.debug          = 0; %if 1 the code sets default values, if 0 all windows will be opened

tag.spchoice       = 1; %if 1 you want to select it

tag.DDMexN         = 1; %if 1 it exports the entire time series

tag.saveCoverMap   = 0; %1 to save the .mat file of the cover map resized over the site

tag.saveReflectivityMat = 0;

tag.SNR            = 1; %1 to compute the SNR and save it in the inOutReferenceFile

tag.Cov            = 0; %if 1 the covariance module and ddm through convolution will be evaluated

tag.PattInterp     = 1; %1 to interpolate all antenna patterns

tag.plotPatterns   = 0;

tag.NoNoise        = 0; %if 1 the no noise configuration file will be used.


if tag.GUI == 1          %% added by ansha for not opening GUI and loading the default previous values from conf folder.
    closeIniGUI = 1;
%     filename='../conf/UserInputs.mat'
%      load('../conf/UserInputs.mat');
     load(filenamedefault);
      if mode1 == "Single Point" || mode1 == "Trackwise"
        disp('.......................................')
        disp(' ### SIMULATION CANNOT PROCEED WITH SINGLE POINT OR TRACKWISE MODE ### ')
        disp('........................................')  
        goToSimulation=0;
        close all
        return
      elseif  mode1 == "Monte Carlo" || mode1 == "Global"
        goToSimulation=1;
      end
else
     disp('.......................................')
     disp(' ### No Input file present opening the main GUI### ')
     disp('........................................') 
     closeIniGUI=UI_HSAVERS_FirstGUI(tag);
     if closeIniGUI==1
            disp('.......................................')
            disp(' ### SIMULATION STOPPED BY THE USER ### ')
            disp('........................................') 
            return
     end
 end

 %% AMIR
 tag.directsignalonly = 0;
 %% AMIR

%     tag.debug          = 1;  %if 1 the code sets default values, if 0 all windows will be opened
%% ------ OPEN THE FIRST (MAIN) INTERACTIVE WINDOW -----------------------
%closeIniGUI=UI_HSAVERS_FirstGUI(tag); % commented by ansha
% if closeIniGUI==1
%     disp('.......................................')
%     disp(' ### SIMULATION STOPPED BY THE USER ### ')
%     disp('........................................') 
%     return
% end

if isempty(data_directory) | data_directory == 0
    errordlg('Data directory not found','Error directory');
    disp('.......................................')
    disp(' ### SIMULATION STOPPED BECAUSE OF AN ERROR ### ')
    disp('........................................')
    error('Fatal error')
    %return
end

if isempty(outfile_directory) | outfile_directory == 0
    errordlg('Output directory not found','Error directory');
    disp('.......................................')
    disp(' ### SIMULATION STOPPED BECAUSE OF AN ERROR ### ')
    disp('........................................')
    error('Fatal error')
    %return
end

%% ------- Orbit(s) selection ----------
 switch  satellite  %  ADDED BY ANSHA 
     case 'HydroGNSS-1'
            filenameData = [data_directory, '\OrbitData\', orbitdataHG1];
            indexPlatform = 1;
            if orbitdataHG1 == 0
                errordlg('Input file not found','File error');
                disp('.......................................')
                disp(' ### SIMULATION STOPPED BECAUSE OF AN ERROR ### ')
                disp('........................................')
                error('Fatal error')
                %return
            end
     case  'HydroGNSS-2'
            filenameData = [data_directory, '\OrbitData\', orbitdataHG2]; 
            indexPlatform = 2;
            if orbitdataHG2 == 0      %% ADDED BY ANSHA               
                errordlg('Input file not found','File error');
                disp('..................TRAp.Rx_pitch.....................')
                disp(' ### SIMULATION STOPPED BECAUSE OF AN ERROR ### ')
                disp('........................................')
                error('Fatal error')
                %return
            end
 end

%% -----check if satellite 1 or satellite 2
% if contains(filenameData,'Sat2')
%     indexPlatform = 2;
% else
%     indexPlatform = 1;
% end

%% --------- load the configuration file ---------

HSAVERSConfigurationFilename = '../conf/HSAVERSConfiguration.xlsx';
if tag.NoNoise == 0
    if indexPlatform == 1
        InstrumentConfigurationFilename = '../conf/InstrumentConfiguration_L1E1_L5E5_Sat1.xlsx';
    else
        InstrumentConfigurationFilename = '../conf/InstrumentConfiguration_L1E1_L5E5_Sat2.xlsx';
    end
else
    InstrumentConfigurationFilename = '../conf/InstrumentConfiguration_NoNoise_L1E1_L5E5.xlsx';
end

[configParam,instrumentConfParam] = readConfigurationFile(HSAVERSConfigurationFilename,InstrumentConfigurationFilename);

%% --------- operational mode --------------------
OperationModeString_temp = string(OperationMode);
OperationModeString      = OperationModeString_temp(OperationModeString_temp ~= "0");

tag.varibleBioGeoInputs         = 0;
DEMinfo.flatsurface             = 0;
tag.MC                          = 0;
tag.CoarseMode                  = 0;
geoSYSp.Mode_extra              = 'Extended simulation area';

% sensitivity analysis mode, variable bio geo inputs
indexVerboseMode = find(OperationModeString=='Variable bio/geo params');
if ~isempty(indexVerboseMode)
    tag.varibleBioGeoInputs       = 1;
    geoSYSp.Mode_extra              = 'Variable bio/geo parameters';
end
%smaller simulation area to shorten the simulation time
indexCoarseMode = find(OperationModeString=='Coarse mode');
if ~isempty(indexCoarseMode)
    tag.CoarseMode               = 1; %some simplification
    %Total Size of the side of the simulation area under consideration
    geoSYSp.AsideMC = configParam.sizeIntegrationAreaCoarseMode;
    geoSYSp.Mode_extra              = 'Coarse mode';
end
%
indexFlatTopo = find(OperationModeString=='Flat Topography');
if ~isempty(indexFlatTopo)
    DEMinfo.flatsurface  = 1; %1 flat surface; DEM elevation and surface curvature set to zero
    geoSYSp.Mode_extra              = 'Flat Topography';
end

geoSYSp.Mode = mode1;

if strcmp(geoSYSp.Mode,'Monte Carlo') == 1
    tag.MC  = 1;
end

if strcmp(Sampling_all_Text{2},'All SPs')
    tag.AllSP=1;
else
    tag.AllSP=0;
end


%% --------- tags for plots -------------------
OutputPlotTextString = string(OutputPlotText);
%OutputPlotText(1) is Map of selected SPs
if OutputPlotTextString(1) ~= "0"
    tag.plotSPsMap = 1;
else
    tag.plotSPsMap = 0;
end
%OutputPlotText(2) is DEM
if OutputPlotTextString(2) ~= "0"
    tag.plotDEM = 1;   
else
    tag.plotDEM = 0;
end
%OutputPlotText(3) is Cover map
if OutputPlotTextString(3) ~= "0"
    tag.plotCoverMap = 1;
else
    tag.plotCoverMap = 0;
end
%OutputPlotText(4) is Angles & Isolines
if OutputPlotTextString(4) ~= "0"
    tag.plotAngles_Isolines = 1;
else
    tag.plotAngles_Isolines = 0;
end
%OutputPlotText(5) is Sigma0
if OutputPlotTextString(5) ~= "0"
    tag.plotSigma0 = 1;
else
    tag.plotSigma0 = 0;
end
%OutputPlotText(6) is DDM
if OutputPlotTextString(6) ~= "0"
    tag.plotDDM = 1;
else
    tag.plotDDM = 0;
end


%% Selection of Tx features (the signal can be either GPS L1 or GALILEO E5a/GPS L5 or Galileo E1C)

geoSYSp.TxSignalButton = signal;
TxPowerError_dBW       = TxAnswers(1);
Tx_phi_roll            = TxAnswers(2);
Tx_theta_pitch         = TxAnswers(3);
Tx_psi_yaw             = TxAnswers(4);
Tx_Euler_phi_Ant       = TxAnswers(5);
Tx_Euler_theta_Ant     = TxAnswers(6);
Tx_Euler_psi_Ant       = TxAnswers(7);
buttonTxMismatch       = Txswitch;

TRAp.TxPowerError_dBW   = TxPowerError_dBW;
%TRAp.TxSignalPower_GUI  = TRAp.TxPower_dBW + TxPowerError_dBW;

%Set transmitter working frequency (GHz), power (dBW), antenna pattern
%according to the selected Tx system.
switch geoSYSp.TxSignalButton
    case 'GPS L1'
       geoSYSp.indexTxSignal   = 0;
       TRAp.TxFrequency        = configParam.frequencyGPSL1_GalileoE1c_GHz;
       TRAp.TxPower_dBW        = configParam.TxPower_L1_dBW;
       TxGainPattern           = configParam.TxAntennaPattern_L1L5;

    case 'GPS L1&L5'
       geoSYSp.indexTxSignal   = 1;
       TRAp.TxFrequency        = [configParam.frequencyGPSL1_GalileoE1c_GHz configParam.frequencyGPSL5_GalileoE5a_GHz];
       TRAp.TxPower_dBW        = [configParam.TxPower_L1_dBW configParam.TxPower_L5_dBW_temp];
       TxGainPattern           = configParam.TxAntennaPattern_L1L5;
%     case 'GPS L5'
%        TRAp.TxFrequency      = configParam.frequencyGPSL5_GalileoE5a_GHz;
%        geoSYSp.indexTxSignal = 1;
%        TxPower_dBW           = configParam.TxPower_L5_dBW_temp;
%        TxGainPattern         = configParam.TxAntennaPattern_L5_temp;

    case 'Galileo E1c'
       geoSYSp.indexTxSignal = 2;
       TRAp.TxFrequency        = configParam.frequencyGPSL1_GalileoE1c_GHz;
       TRAp.TxPower_dBW        = configParam.TxPower_E1bc_dBW;
       TxGainPattern           = configParam.TxAntennaPattern_E1E5;

    case 'Galileo E1c&E5a'
       geoSYSp.indexTxSignal = 3;
       TRAp.TxFrequency        = [configParam.frequencyGPSL1_GalileoE1c_GHz configParam.frequencyGPSL5_GalileoE5a_GHz];
       TRAp.TxPower_dBW        = [configParam.TxPower_E1bc_dBW configParam.TxPower_E5a_dBW];
       TxGainPattern           = configParam.TxAntennaPattern_E1E5;
%     case 'Galileo E1c'
%        TRAp.TxFrequency      = configParam.frequencyGPSL1_GalileoE1c_GHz;
%        geoSYSp.indexTxSignal = 3;
%        TxPower_dBW           = configParam.TxPower_E1bc_dBW;
%        TxGainPattern         = configParam.TxAntennaPattern_E1;
end

%Defintion of the wavelength
lightSpeed = 0.299792458; %m/ns
TRAp.waveL = lightSpeed./TRAp.TxFrequency;

%Set Tx Pol Mismatch
switch buttonTxMismatch
    case 'Ideal'
        TRAp.flagTxMismatch = 0;
    case 'Mismatch'
        TRAp.flagTxMismatch = 1;
end

%DC: added to consider attitude of Tx - not controlled by GUI but needed
%for the structure of the code
TRAp.Tx_nomAtt_angles      = configParam.Tx_nomAtt_angles;
TRAp.TxAttitude_angles_GUI = [Tx_phi_roll, Tx_theta_pitch, Tx_psi_yaw];
TRAp.TxEuler_angles        = [Tx_Euler_phi_Ant, Tx_Euler_theta_Ant, Tx_Euler_psi_Ant];

%% ------- Rx settings --------

Rx_phi_roll          = RxAnswers(1);
Rx_theta_pitch       = RxAnswers(2);
Rx_psi_yaw           = RxAnswers(3);
Rx_nomAtt_roll       = RxAnswers(4); 
Rx_nomAtt_pitch      = RxAnswers(5);
Rx_nomAtt_yaw        = RxAnswers(6);
buttonRxDOWNmismatch = Rxswitch_Nadir;
buttonRxUPmismatch   = Rxswitch_Zenith;

%DC: added to consider attitude of Rx
TRAp.Rx_nomAtt_angles      = [Rx_nomAtt_roll, Rx_nomAtt_pitch, Rx_nomAtt_yaw];
TRAp.RxAttitude_angles_GUI = [Rx_phi_roll, Rx_theta_pitch, Rx_psi_yaw];
TRAp.RxEuler_angles        = configParam.RxEuler_angles;

%Set Rx Down Antenna Polarization Mismatch 
switch buttonRxDOWNmismatch
    case 'Ideal'
        TRAp.flagRxDOWNmismatch = 0;
    case 'NadirMismatch'
        TRAp.flagRxDOWNmismatch = 1;
end

%Set Rx UP Antenna Polarization Mismatch 
switch buttonRxUPmismatch
    case 'Ideal'
        TRAp.flagRxUPmismatch = 0;
    case 'ZenithMismatch'
        TRAp.flagRxUPmismatch = 1;
end


%% Read Rx DOWN antenna pattern
switch  satellite  %  ADDED BY ANSHA 
     case 'HydroGNSS-1'
            filenameRxDOWNGainPattern = [data_directory, '\AntennaPattern\', nadirantennadataHG1];
            if nadirantennadataHG1 == 0
                errordlg('Input file not found','File error');
                disp('.......................................')
                disp(' ### SIMULATION STOPPED BECAUSE OF AN ERROR ### ')
                disp('........................................')
                error('Fatal error')
                %return
            end
     case  'HydroGNSS-2'
             filenameRxDOWNGainPattern = [data_directory, '\AntennaPattern\', nadirantennadataHG2];
            if nadirantennadataHG2 == 0      %% ADDED BY ANSHA               
                errordlg('Input file not found','File error');
                disp('.......................................')
                disp(' ### SIMULATION STOPPED BECAUSE OF AN ERROR ### ')
                disp('........................................')
                error('Fatal error')
                %return
            end
 end
% filenameRxDOWNGainPattern = [data_directory, '\AntennaPattern\', nadirantennadata];

if indexPlatform == 1 && contains(filenameRxDOWNGainPattern, 'Sat2')
    CreateStruct.Interpreter = 'tex';
    CreateStruct.WindowStyle = 'modal';
    uiwait(msgbox({'WARNING: The platform of the orbit file and of the nadir antenna pattern do not match (SAT1/SAT2).';...
                    'The antenna pattern matching the platform of the orbit file will be used.'},CreateStruct))
    filenameRxDOWNGainPattern = replace(filenameRxDOWNGainPattern,'Sat2','Sat1');
end

if indexPlatform == 2 && contains(filenameRxDOWNGainPattern, 'Sat1')
    CreateStruct.Interpreter = 'tex';
    CreateStruct.WindowStyle = 'modal';
    uiwait(msgbox({'WARNING: The platform of the orbit file and of the nadir antenna pattern do not match (SAT1/SAT2).';...
                    'The antenna pattern matching the platform of the orbit file will be used.'},CreateStruct))
    filenameRxDOWNGainPattern = replace(filenameRxDOWNGainPattern,'Sat1','Sat2');
end

%NOTE: all patterns files should be excel files with the following structure
%table where cell (1,1) is text/empty; first raw, starting from cell (1,2)
%reports theta values; the first column, starting from cell (2,1) reports phi values.
%The gain is reported in a table starting from cell (2,2).
%The gain should be in dB. Theta and phi should be in degree.
%Phi should be wrapped to 180.
%Please see HydroGnss_Nadir_L1E1_L5E5_LHCPandRHCP.xlsx for
%reference.

%The nadir antenna pattern file is an .xls file with 4 spreadsheets, two
%for LHCP, L1/E1 and L5/E5 frequency, and two for RHCP, L1/E1 and L5/E5 frequency.
%Each spreadsheet is a table of 362x92, where the first raw reports theta values and the first column
%reports phi values. The cell 1,1 is text.
%The data is stored as 1 degree steps in theta and phi.
%Theta ranges from 0° to 90° where 0 degree is at boresight.
%Phi ranges from -180° to 180°, where the gain at -180° is equal to the gain at 180°
%The gain is in dBiC.

%read the nadir antenna pattern
optsReadTableNadir=detectImportOptions(filenameRxDOWNGainPattern);
% Read mismatch parameters
optsReadTableNadir.Sheet='MismatchParams' ; 
Mismatch_All_tab= readtable(filenameRxDOWNGainPattern, optsReadTableNadir, "UseExcel", false); 
Mismatch_All=table2array(Mismatch_All_tab);
% LHCP, L1/E1
optsReadTableNadir.Sheet = 'LHCP_L1E1';
nadirAntennaPatternArray_LHCP_tab = readtable(filenameRxDOWNGainPattern, optsReadTableNadir, "UseExcel", false);
nadirAntennaPatternArray_LHCP     = table2array(nadirAntennaPatternArray_LHCP_tab);

%% amir for rongowai
nadirAntennaPatternArray_LHCP(2:end,2:end) = nadirAntennaPatternArray_LHCP(2:end,2:end) - 33;

nadirPattern_azimuth_0       = nadirAntennaPatternArray_LHCP(:,1);
nadirPattern_azimuth         = nadirPattern_azimuth_0(2:length(nadirPattern_azimuth_0));
nadirPattern_offPointing_0   = nadirAntennaPatternArray_LHCP(1,:);
nadirPattern_offPointing     = nadirPattern_offPointing_0(2:length(nadirPattern_offPointing_0));
nadirPattern_gain_linear_LHCP= db2pow(nadirAntennaPatternArray_LHCP(2:length(nadirPattern_azimuth_0),2:length(nadirPattern_offPointing_0)));
Mismatch_LHCP=Mismatch_All(1:2,2) ; 
% RHCP, L1/E1
optsReadTableNadir.Sheet = 'RHCP_L1E1';
nadirAntennaPatternArray_RHCP_tab = readtable(filenameRxDOWNGainPattern, optsReadTableNadir, "UseExcel", false);
nadirAntennaPatternArray_RHCP     = table2array(nadirAntennaPatternArray_RHCP_tab);

%% amir for rongowai
nadirAntennaPatternArray_RHCP(2:end,2:end) = nadirAntennaPatternArray_RHCP(2:end,2:end) - 33;


nadirPattern_gain_linear_RHCP= db2pow(nadirAntennaPatternArray_RHCP(2:length(nadirPattern_azimuth_0),2:length(nadirPattern_offPointing_0)));
Mismatch_RHCP=Mismatch_All(3:4,2) ; 


%% Cross-pol Antenna pattern - Amirreza
% G_LR
optsReadTableNadir.Sheet = 'LR_L1E1';
nadirAntennaPatternArray_GLR_tab = readtable(filenameRxDOWNGainPattern, optsReadTableNadir, "UseExcel", false);
nadirAntennaPatternArray_GLR     = table2array(nadirAntennaPatternArray_GLR_tab);

%% amir for rongowai
nadirAntennaPatternArray_GLR(2:end,2:end) = nadirAntennaPatternArray_GLR(2:end,2:end) - 33;


GLR_nadirPattern_azimuth_0       = nadirAntennaPatternArray_GLR(:,1);
GLR_nadirPattern_azimuth         = GLR_nadirPattern_azimuth_0(2:length(GLR_nadirPattern_azimuth_0));
GLR_nadirPattern_offPointing_0   = nadirAntennaPatternArray_GLR(1,:);
GLR_nadirPattern_offPointing     = GLR_nadirPattern_offPointing_0(2:length(GLR_nadirPattern_offPointing_0));
nadirPattern_gain_linear_GLR= db2pow(nadirAntennaPatternArray_GLR(2:length(GLR_nadirPattern_azimuth_0),2:length(GLR_nadirPattern_offPointing_0)));
% G_RL
optsReadTableNadir.Sheet = 'RL_L1E1';
nadirAntennaPatternArray_GRL_tab = readtable(filenameRxDOWNGainPattern, optsReadTableNadir, "UseExcel", false);
nadirAntennaPatternArray_GRL     = table2array(nadirAntennaPatternArray_GRL_tab);

%% amir for rongowai
nadirAntennaPatternArray_GRL(2:end,2:end) = nadirAntennaPatternArray_GRL(2:end,2:end) - 33;



GRL_nadirPattern_azimuth_0       = nadirAntennaPatternArray_GRL(:,1);
GRL_nadirPattern_azimuth         = GRL_nadirPattern_azimuth_0(2:length(GRL_nadirPattern_azimuth_0));
GRL_nadirPattern_offPointing_0   = nadirAntennaPatternArray_GRL(1,:);
GRL_nadirPattern_offPointing     = GRL_nadirPattern_offPointing_0(2:length(GRL_nadirPattern_offPointing_0));
nadirPattern_gain_linear_GRL= db2pow(nadirAntennaPatternArray_GRL(2:length(GRL_nadirPattern_azimuth_0),2:length(GRL_nadirPattern_offPointing_0)));

%%

%if L1 and L5 (or E1 and E5) are simulated in parallel, read also the L5/E5
%sheet and save it in the (:,:,2) of the pattern array
if geoSYSp.indexTxSignal==1 || geoSYSp.indexTxSignal==3
    %LHCP, L5/E5, to linear unit
    optsReadTableNadir.Sheet = 'LHCP_L5E5';
    nadirAntennaPatternArray_LHCP_5_tab = readtable(filenameRxDOWNGainPattern, optsReadTableNadir, "UseExcel", false);
    nadirAntennaPatternArray_LHCP_5     = table2array(nadirAntennaPatternArray_LHCP_5_tab);
    nadirPattern_gain_linear_LHCP_5= db2pow(nadirAntennaPatternArray_LHCP_5(2:length(nadirPattern_azimuth_0),2:length(nadirPattern_offPointing_0)));
    nadirPattern_gain_linear_LHCP(:,:,2) = nadirPattern_gain_linear_LHCP_5;
    %RHCP, L5/E5, to linear unit
    optsReadTableNadir.Sheet = 'RHCP_L5E5';
    nadirAntennaPatternArray_RHCP_5_tab = readtable(filenameRxDOWNGainPattern, optsReadTableNadir, "UseExcel", false);
    nadirAntennaPatternArray_RHCP_5     = table2array(nadirAntennaPatternArray_RHCP_5_tab);
    nadirPattern_gain_linear_RHCP_5 = db2pow(nadirAntennaPatternArray_RHCP_5(2:length(nadirPattern_azimuth_0),2:length(nadirPattern_offPointing_0)));
    nadirPattern_gain_linear_RHCP(:,:,2) = nadirPattern_gain_linear_RHCP_5;
end

TRAp.phiRxDown          = nadirPattern_azimuth;
TRAp.thetaRxDown        = nadirPattern_offPointing;
TRAp.GainRxDown_LHCP    = nadirPattern_gain_linear_LHCP;
TRAp.GainRxDown_RHCP    = nadirPattern_gain_linear_RHCP;
TRAp.GainRxDown_LR    = nadirPattern_gain_linear_GLR;
TRAp.GainRxDown_RL    = nadirPattern_gain_linear_GRL;

TRAp.MM_LHCP=Mismatch_LHCP ;
TRAp.MM_RHCP=Mismatch_RHCP ; 

[NadirPattern_offPointing, NadirPattern_azimuth] = meshgrid(nadirPattern_offPointing, nadirPattern_azimuth);

%pattern interpolation
if tag.PattInterp == 1
    nadirPattern_offPointing_new  =  0:configParam.dang:max(max(NadirPattern_offPointing));
    nadirPattern_azimuth_new      = -180:configParam.dang:180;
    [NadirPattern_offPointing_new, NadirPattern_azimuth_new_new] = meshgrid(nadirPattern_offPointing_new, nadirPattern_azimuth_new);

    nadirPattern_gain_linear_LHCP_new = interp2(NadirPattern_offPointing, NadirPattern_azimuth,...
                                                nadirPattern_gain_linear_LHCP(:,:,1), NadirPattern_offPointing_new, NadirPattern_azimuth_new_new);
    nadirPattern_gain_linear_RHCP_new = interp2(NadirPattern_offPointing, NadirPattern_azimuth,...
                                                nadirPattern_gain_linear_RHCP(:,:,1), NadirPattern_offPointing_new, NadirPattern_azimuth_new_new);
  
    nadirPattern_gain_linear_GLR_new = interp2(NadirPattern_offPointing, NadirPattern_azimuth,...
                                                nadirPattern_gain_linear_GLR(:,:,1), NadirPattern_offPointing_new, NadirPattern_azimuth_new_new);

    nadirPattern_gain_linear_GRL_new = interp2(NadirPattern_offPointing, NadirPattern_azimuth,...
                                                nadirPattern_gain_linear_GRL(:,:,1), NadirPattern_offPointing_new, NadirPattern_azimuth_new_new);


    %if L1 and L5 (or E1 and E5) are simulated in parallel, interpolate also the L5/E5
    %pattern and save it in the (:,:,2) of the interp pattern array
    if geoSYSp.indexTxSignal==1 || geoSYSp.indexTxSignal==3
        nadirPattern_gain_linear_LHCP_new(:,:,2) = interp2(NadirPattern_offPointing, NadirPattern_azimuth,...
                                                    nadirPattern_gain_linear_LHCP(:,:,2), NadirPattern_offPointing_new, NadirPattern_azimuth_new_new);

        nadirPattern_gain_linear_RHCP_new(:,:,2) = interp2(NadirPattern_offPointing, NadirPattern_azimuth,...
                                                    nadirPattern_gain_linear_RHCP(:,:,2), NadirPattern_offPointing_new, NadirPattern_azimuth_new_new);
    end

    %overwrite with interpolated variables
    TRAp.phiRxDown          = nadirPattern_azimuth_new;
    TRAp.thetaRxDown        = nadirPattern_offPointing_new;
    TRAp.GainRxDown_LHCP    = nadirPattern_gain_linear_LHCP_new;
    TRAp.GainRxDown_RHCP    = nadirPattern_gain_linear_RHCP_new;

    TRAp.GainRxDown_LR    = nadirPattern_gain_linear_GLR_new;
    TRAp.GainRxDown_RL    = nadirPattern_gain_linear_GRL_new;
end

if tag.plotPatterns == 1
   figure
   imagesc(NadirPattern_offPointing,NadirPattern_azimuth,nadirPattern_gain_linear_LHCP(:,:,1))
   ylabel('\phi'), xlabel('\theta'), colorbar, title('Rx DOWN pattern for LHCP L1/E1') 

   figure
   imagesc(nadirPattern_offPointing_new, nadirPattern_azimuth_new, nadirPattern_gain_linear_LHCP_new(:,:,1))
   ylabel('\phi'), xlabel('\theta'), colorbar, title('Rx DOWN interp. pattern for LHCP L1/E1') 
end


%% Read Rx UP antenna pattern
switch  satellite  %  ADDED BY ANSHA 
     case 'HydroGNSS-1'
            filenameRxUPGainPattern   = [data_directory, '\AntennaPattern\', zenithantennadataHG1];
            if zenithantennadataHG1 == 0
                errordlg('Input file not found','File error');
                disp('.......................................')
                disp(' ### SIMULATION STOPPED BECAUSE OF AN ERROR ### ')
                disp('........................................')
                error('Fatal error')
                %return
            end
     case  'HydroGNSS-2'
            filenameRxUPGainPattern   = [data_directory, '\AntennaPattern\', zenithantennadataHG2];
            if zenithantennadataHG2 == 0      %% ADDED BY ANSHA               
                errordlg('Input file not found','File error');
                disp('.......................................')
                disp(' ### SIMULATION STOPPED BECAUSE OF AN ERROR ### ')
                disp('........................................')
                error('Fatal error')
                %return
            end
 end
% filenameRxUPGainPattern   = [data_directory, '\AntennaPattern\', zenithantennadata];

if indexPlatform == 1 && contains(filenameRxUPGainPattern, 'Sat2')
    CreateStruct.Interpreter = 'tex';
    CreateStruct.WindowStyle = 'modal';
    uiwait(msgbox({'WARNING: The platform of the orbit file and of the zenith antenna pattern do not match (SAT1/SAT2).';...
                    'The antenna pattern matching the platform of the orbit file will be used.'},CreateStruct))
    filenameRxUPGainPattern = replace(filenameRxUPGainPattern,'Sat2','Sat1');
end

if indexPlatform == 2 && contains(filenameRxUPGainPattern, 'Sat1')
    CreateStruct.Interpreter = 'tex';
    CreateStruct.WindowStyle = 'modal';
    uiwait(msgbox({'WARNING: The platform of the orbit file and of the zenith antenna pattern do not match (SAT1/SAT2).';...
                    'The antenna pattern matching the platform of the orbit file will be used.'},CreateStruct))
    filenameRxUPGainPattern = replace(filenameRxUPGainPattern,'Sat1','Sat2');
end

%NOTE: all patterns files should be excel files with a predefined structure

%The zenith antenna pattern file is an .xls file with two spreadsheets,
% one for L1 and E1 frequency and one for L5 and E5 frequency. 
%It is a table of 362x92, where the first raw reports theta values and the first column
%reports phi values. The cell 1,1 is text.
%The data is stored as 1 degree steps in theta and phi.
%Theta ranges from 0° to 90° where 0 degree is at boresight.
%Phi ranges from -180° to 180°, where the gain at -180° is equal to the gain at 180°
%The gain is in dBiC.

%read the zenith antenna pattern
optsReadTableZenith=detectImportOptions(filenameRxUPGainPattern);

%read L1/E1 pattern, to linear unit
optsReadTableZenith.Sheet = 'L1E1';
zenithAntennaPatternArray_tab = readtable(filenameRxUPGainPattern, optsReadTableZenith, "UseExcel", false);
zenithAntennaPatternArray     = table2array(zenithAntennaPatternArray_tab);
zenithPattern_azimuth_0       = zenithAntennaPatternArray(:,1);
zenithPattern_azimuth         = zenithPattern_azimuth_0(2:length(zenithPattern_azimuth_0));
zenithPattern_offPointing_0   = zenithAntennaPatternArray(1,:);
zenithPattern_offPointing     = zenithPattern_offPointing_0(2:length(zenithPattern_offPointing_0));
zenithPattern_gain_linear     = db2pow(zenithAntennaPatternArray(2:length(zenithPattern_azimuth_0),2:length(zenithPattern_offPointing_0)));

%if L1 and L5 (or E1 and E5) are simulated in parallel, read also the L5/E5
%sheet and save it in the (:,:,2) of the pattern array
if geoSYSp.indexTxSignal==1 || geoSYSp.indexTxSignal==3
    optsReadTableZenith.Sheet = 'L5E5';
    zenithAntennaPatternArray_5_tab = readtable(filenameRxUPGainPattern, optsReadTableZenith, "UseExcel", false);
    zenithAntennaPatternArray_5     = table2array(zenithAntennaPatternArray_5_tab);
    zenithPattern_gain_linear(:,:,2)= db2pow(zenithAntennaPatternArray_5(2:length(zenithPattern_azimuth_0),2:length(zenithPattern_offPointing_0)));
end

TRAp.theta_RxUPpattern = zenithPattern_offPointing;
TRAp.phi_RxUPpattern   = zenithPattern_azimuth;
TRAp.GainRxUP          = zenithPattern_gain_linear;

%interpolation of the pattern
if tag.PattInterp == 1

    [ZenithPattern_offPointing, ZenithPattern_azimuth] = meshgrid(zenithPattern_offPointing, zenithPattern_azimuth);

    zenithPattern_offPointing_new =  0:configParam.dang:90;
    zenithPattern_azimuth_new     = -180:configParam.dang:180;
    [ZenithPattern_offPointing_new, ZenithPattern_azimuth_new] = meshgrid(zenithPattern_offPointing_new, zenithPattern_azimuth_new);

    zenithPattern_gain_linear_new = interp2(ZenithPattern_offPointing, ZenithPattern_azimuth, ...
                zenithPattern_gain_linear(:,:,1), ZenithPattern_offPointing_new, ZenithPattern_azimuth_new);

    %if L1 and L5 (or E1 and E5) are simulated in parallel, interpolate also the L5/E5
    %pattern and save it in the (:,:,2) of the pattern array
    if geoSYSp.indexTxSignal==1 || geoSYSp.indexTxSignal==3
        zenithPattern_gain_linear_new(:,:,2) = interp2(ZenithPattern_offPointing, ZenithPattern_azimuth, ...
                zenithPattern_gain_linear(:,:,2), ZenithPattern_offPointing_new, ZenithPattern_azimuth_new);
    end

    if tag.plotPatterns == 1
        figure
        imagesc(zenithPattern_offPointing,zenithPattern_azimuth,zenithPattern_gain_linear(:,:,1))
        ylabel('\phi'), xlabel('\theta'), colorbar, title('original Rx UP pattern for L1/E1')

        figure
        imagesc(ZenithPattern_offPointing_new(:),ZenithPattern_azimuth_new(:),zenithPattern_gain_linear_new(:,:,1))
        ylabel('\phi'), xlabel('\theta'), colorbar, title('interp Rx UP pattern for L1/E1')
    end

    %overwrite with interpolated variables
    TRAp.theta_RxUPpattern = zenithPattern_offPointing_new;
    TRAp.phi_RxUPpattern   = zenithPattern_azimuth_new;
    TRAp.GainRxUP          = zenithPattern_gain_linear_new;
end


%% ------------Rx NOISE INPUTS---------------------------

geoSYSp.cohIntTime           = SensorNoise(1);
geoSYSp.incoIntTime          = SensorNoise(2); 
geoSYSp.cohIntTime_sec       = 1e-3*SensorNoise(1); %must be in seconds in HSAVERS
geoSYSp.incoIntTime_sec      = 1e-3*SensorNoise(2); %must be in seconds in HSAVERS
geoSYSp.DeltaF_interf_dB     = SensorNoise(3);
geoSYSp.TA_down_K            = SensorNoise(4);
geoSYSp.TA_up_K              = SensorNoise(5);
geoSYSp.DeltaTRx_interf_K    = SensorNoise(6);
geoSYSp.DeltaG_sys_interf_dB = SensorNoise(7);

geoSYSp.CodeRepetitionFreq  = configParam.CodeRepetitionFreq_Hz;
geoSYSp.alfa                = degreeofcoherence;

% read two parameters from the instrumentData configuration file
% this is needed to write then in the metadata
geoSYSp.TRx_K         = instrumentConfParam.TRx_K;
geoSYSp.TL_K          = instrumentConfParam.TL_K;


%% -----simulation resolution and DDM settings--------------

geoSYSp.dx                 = SimulationSetting(1);
geoSYSp.dy                 = SimulationSetting(1);
geoSYSp.FRQmax             = SimulationSetting(2);
Dly_RNGmax_L1E1            = SimulationSetting(3);
Dly_RNGmax_L5E5            = SimulationSetting(4);
geoSYSp.nx_Rdelay          = SimulationSetting(5);
geoSYSp.ny_Dfreq           = SimulationSetting(6);
TRAp.TxBeamwidth           = SimulationSetting(7);
TRAp.RxBeamwidth           = SimulationSetting(8);
geoSYSp.sigPowerBooster_dB = SimulationSetting(9);

geoSYSp.Dly_RNGmax     = [Dly_RNGmax_L1E1 Dly_RNGmax_L5E5];

%Hertz freq. step for mapping [Hz] - Ratio
% geoSYSp.fpas = 2*geoSYSp.FRQmax/geoSYSp.ny_Dfreq;
geoSYSp.fpas = geoSYSp.FRQmax / geoSYSp.ny_Dfreq; % Added by amir for Rongowai
%timerange step for mapping [usec], two values array, for L1/E1 and L5/E5
geoSYSp.mpas = geoSYSp.Dly_RNGmax./double(geoSYSp.nx_Rdelay);

%compute delay offset from inputs and using the delay offset factor in the
%configuration file
geoSYSp.delayOffset_bin         = double(floor(geoSYSp.nx_Rdelay*configParam.delayOffsetFactor));
geoSYSp.delayOffset_microsec    = geoSYSp.delayOffset_bin.*geoSYSp.Dly_RNGmax./geoSYSp.nx_Rdelay; %for L1/E1 and for L5/E5

%% AMIR
% geoSYSp.chipTime = [1.0 0.1]; %[L1/E1 L5/E5] chip time in microsec [usec]
% int_mpas_L1 = max(1, round(geoSYSp.nx_Rdelay * geoSYSp.chipTime(1) /geoSYSp.Dly_RNGmax(1)));  % integer bins per chip
% int_mpas_L5 = max(1, round(geoSYSp.nx_Rdelay * geoSYSp.chipTime(2) /geoSYSp.Dly_RNGmax(2)));  % integer bins per chip
% mpas_L1_int = geoSYSp.chipTime(1)/int_mpas_L1;
% mpas_L5_int = geoSYSp.chipTime(2)/int_mpas_L5;
% 
% clear geoSYSp.mpas
% geoSYSp.mpas = [mpas_L1_int mpas_L5_int];
% geoSYSp.delayOffset_microsec = geoSYSp.delayOffset_bin * geoSYSp.mpas;
%% AMIR


%CodeDelaySpacingSamplesBetweenPixels
geoSYSp.SamplingFrequency_Hz                 = configParam.SamplingFrequency_Hz;
geoSYSp.CodeDelaySpacingSamplesBetweenPixels = round(geoSYSp.mpas.*(configParam.SamplingFrequency_Hz/10^6));

%width of the geobounded window around nominal SP where to search the peak
geoSYSp.geoboundedWindowWidthInDelay_numBins = floor(configParam.geoboundedWindowWidth(1)./geoSYSp.mpas); %two values array for L1/E1 and for L5/E5
if geoSYSp.Dly_RNGmax(2)==0.0
    geoSYSp.geoboundedWindowWidthInDelay_numBins(2) = 0.0;
end
geoSYSp.geoboundedWindowWidthInFreq_numBins  = floor(configParam.geoboundedWindowWidth(2)./geoSYSp.fpas);

%% -------- ionization parameters --------------------------
filenameIonoParam   = [data_directory, '\', configParam.filenameIonisationLevelParam];
%read the ionization parameters
optsReadTableIono=detectImportOptions(filenameIonoParam);
switch configParam.SolarActivityLevel
    case 'high'
        optsReadTableIono.Sheet = 'high';
    case 'medium'
        optsReadTableIono.Sheet = 'medium';
    case 'low'
        optsReadTableIono.Sheet = 'low';
end
ionoParam_tab        = readtable(filenameIonoParam, optsReadTableIono, "UseExcel", false);
ionisationLevelParam = table2array(ionoParam_tab);

%% --------- select the site and define the DEM/maps ---------------------
if tag.debug == 1   
    soilmoisture = 'Tibesti Mountains';
end

[filenameDEM30arcsec,filenameDEM,filenameDEM_flooded,SPcoodinatesLimits,tag.scenario] = siteDependentInputs(soilmoisture,forest,freeze_thaw,wetlands,geoSYSp.dx,geoSYSp.dy,configParam);

DEMinfo.filenameDEM30arcsec = [data_directory, '\DEM\', filenameDEM30arcsec];
DEMinfo.filenameDEM         = [data_directory, '\DEM\', filenameDEM];
DEMinfo.filenameDEM_flooded = [data_directory, '\DEM\', filenameDEM_flooded];
DEMinfo.filenameDEM5km      = [data_directory, '\DEM\', configParam.filenameDEM5km];

%SPcoodinatesLimits=[26.15 26.3 -80.7 -79.5]; %for test global mode on Everglades
DEMinfo.SPcoodinatesLimits  = SPcoodinatesLimits;
filenameCoverMap            = [data_directory, '\LandCoverMap\',configParam.filenameCoverMap];


%Values in degree
%Distance between the SP and the border of the map
%This is used to define the size of the window used to open and process only a part of the DEM,
%land cover map and other maps eventually given as input
%If one SP is simulated then a map window of 6° by 6° centered on the SP will be open.
%NOTE all SPs should be at least 3° far from all borders (north, east, and west
%south) of the map
%in coarse mode a smaller simulation area is used therefore a smaller map
%window is necessary (distance SP borders 1°)
%for high latitudes sites, the window map must be larger

if tag.CoarseMode == 0
    if SPcoodinatesLimits(1) < 60
        DEMinfo.sizeInputMapWindow = configParam.sizeInputMapWindow;
    else
        DEMinfo.sizeInputMapWindow = configParam.sizeInputMapWindowHighLatitudes;
    end
else
    if SPcoodinatesLimits(1) < 60
        DEMinfo.sizeInputMapWindow = configParam.sizeInputMapWindowCoarseMode;
    else
        DEMinfo.sizeInputMapWindow = configParam.sizeInputMapWindowHighLatitudesCoarseMode;
    end
end

%% --------- output directory ---------------------

%%%Definition of a destination directory taking current date and time
date_hour          = datestr(now);
date_hour([12,15,18]) = '_';
name_data_tag=[RUNID '_' date_hour ];
 satellite = string(satellite); % to solve for temporary
dest_dir           = [outfile_directory '\' name_data_tag '\' satellite '\'];
 mkdir(strjoin(dest_dir,''));
%  mkdir(dest_dir);
%%  read the orbit file and select the specular point(s) to simulate

if tag.debug == 0
    %%%Position and speeds are in ECEF
    [Tx, VTx, Rx, VRx, outputFilename, SPlatlonAlt, incAngle, PRN, GNSS_Week, GNSS_SoW,TrackID,SVN,goToSimulation,OrbitWithReflectivity,ReflectivitydB_selected,tag,Rx_roll_selected,Rx_pitch_selected,Rx_yaw_selected,Rx_Attitude_Variation_RPY]...
                        = selectAndReadGNSSRdata_newOrbitFile(filenameData,DEMinfo,filenameCoverMap, dest_dir, mode1, geoSYSp,Sampling,Sampling_all_Text,tag);
    %%added trackid -by aneesha

    %saving the vector within a structure to keep track of them              
    geoSYSp.SPlat_orbitFile = SPlatlonAlt(:,1);   
    geoSYSp.SPlon_orbitFile = SPlatlonAlt(:,2);
    geoSYSp.SPAlt_orbitFile = SPlatlonAlt(:,3);
    geoSYSp.incAngle_orbitFile = incAngle;

    geoSYSp.Rx_roll = Rx_roll_selected; 
    geoSYSp.Rx_pitch = Rx_pitch_selected; 
    geoSYSp.Rx_yaw = Rx_yaw_selected; 
    geoSYSp.Rx_Attitude_Variation_RPY = Rx_Attitude_Variation_RPY;

    if geoSYSp.Rx_Attitude_Variation_RPY(1) == 0
       geoSYSp.Rx_roll = repmat(TRAp.Rx_nomAtt_angles(1), 1, size(Rx, 1));
    end
    if geoSYSp.Rx_Attitude_Variation_RPY(2) == 0
        geoSYSp.Rx_pitch = repmat(TRAp.Rx_nomAtt_angles(2), 1, size(Rx, 1));
    end
    if geoSYSp.Rx_Attitude_Variation_RPY(3) == 0
        geoSYSp.Rx_yaw = repmat(TRAp.Rx_nomAtt_angles(3), 1, size(Rx, 1));
    end

else
    %DEFAULT SP FOR DEBUG PHASE
    Tx  = [-9.731997352920022e+06 -1.763080636031368e+07  2.170164458069115e+07];
    VTx = [ 2.061569403093482e+03   .104969783279554e+02  1.502141428956273e+03];
    Rx  = [-1.432660211498533e+06 -5.453296873849013e+06  4.028154566084063e+06];
    VRx = [-2.692565968365223e+03 -3.787609703848653e+03 -6.086552372852673e+03];
    PRN = 36;
    outputFilename       = cellstr(['PRN' num2str(PRN) '_SP' num2str(1)]);
    SPlatlon             = [36.872484669606050 -105.8682734932170];    
    goToSimulation       = 1;
end

if goToSimulation == 1
    %% ------- READ THE TX PATTERN OF THE SELECTED PRN ------------------------

    % NOTE:at this moment the file of the Tx antenna pattern is fixed.
    % This part should be changed once we have all the Tx patterns !
    filenameTxGainPattern = [data_directory,'\AntennaPattern\', TxGainPattern];

    %NOTE: all patterns files should be excel files with a predefined structure

    %The GPS antenna pattern file is an .xls file with two spreadsheets,
    % one for L1 frequency and one for L2 frequency. GPS L2 pattern is used
    % instead of GPS L5 pattern, they are very similar, but the latter is
    % not available.
    %It is a table of 37x92, where the first raw reports theta values and the first column
    %reports phi values. The cell 1,1 is text.
    %The data is stored as 2 degree steps in theta and 10 degree step in phi.
    %Theta ranges from -90° to 90° but the pattern is almost symmetric around 0°.
    %Phi ranges from 0° to 350°.
    %The gain is in dB.
    %The GPS antenna pattern file an .xls file with two spreadsheets,
    % one for E1 frequency and one for E5 frequency.
    %It is a table of 721x92,where the first raw reports theta values and the first column
    %reports phi values. The cell 1,1 is text.
    %The data is stored as 0.5 degree steps in theta and in phi.
    %Theta ranges from 0° to 45° for E1 and from 0° to 52° for E5
    %Phi ranges from 0° to 359.5°.
    %The gain is in dB.    

    optsReadTableTxPattern=detectImportOptions(filenameTxGainPattern);

    if geoSYSp.indexTxSignal == 0 || geoSYSp.indexTxSignal == 1
        %read GPS L1 pattern
        optsReadTableTxPattern.Sheet = 'GPS_IIR_SVN61L1';    
        TxAntennaPatternArray_tab = readtable(filenameTxGainPattern, optsReadTableTxPattern, "UseExcel", false);
        TxAntennaPatternArray     = table2array(TxAntennaPatternArray_tab);

        %reshape the variable in order to have (as it is for the Rx) the gain in a vector 37x91,
        %phi in a vector 37x1, theta in a vector 91x1
        TxPattern_phi_temp_0   = TxAntennaPatternArray(:,1);
        TxPattern_phi          = TxPattern_phi_temp_0(2:length(TxPattern_phi_temp_0));
        TxPattern_theta_temp_0 = TxAntennaPatternArray(1,:);
        TxPattern_theta_temp   = TxPattern_theta_temp_0(2:length(TxPattern_theta_temp_0));
        TxPattern_gain         = TxAntennaPatternArray(2:length(TxPattern_phi_temp_0), 2:length(TxPattern_theta_temp_0));
        %make the format of the pattern similar to that of Rx (to be
        %coherent with the used frame)
        %wrap phi to 180° and and sort to have phi ranging from -180° to 180°
        TxPattern_phi_180                     = wrapTo180(TxPattern_phi);
        [TxPattern_phi_180_sort,indexPhiSort] = sort(TxPattern_phi_180);
        TxPattern_gain_phiSort                   = TxPattern_gain(indexPhiSort,:);
        %add a row in phi and pattern corresponding to -180° (equal to the
        %row corresponding to 180°) to improve the performance of the
        %interpolation
        TxPattern_phi_180_sort_complete(1)      = -180.0;
        TxPattern_phi_180_sort_complete(2:length(TxPattern_phi_180_sort)+1) = TxPattern_phi_180_sort;
        TxPattern_gain_phiSort_complete(1,:) = TxPattern_gain_phiSort(length(TxPattern_phi_180_sort),:);
        TxPattern_gain_phiSort_complete(2:length(TxPattern_phi_180_sort)+1,:) = TxPattern_gain_phiSort;

        %NOTE only half of the pattern is used for the simulations
        %theta -90° to 0° (absolute value) ;
        indexResizedTheta       = find(TxPattern_theta_temp <= 0);
        TxPattern_theta_res     = abs(TxPattern_theta_temp(indexResizedTheta));
        TxPattern_gain_res = TxPattern_gain_phiSort_complete(:,indexResizedTheta);
        
        if geoSYSp.indexTxSignal == 1
            %read GPS L5 pattern
            optsReadTableTxPattern.Sheet = 'GPS_IIR_SVN61L2';    
            TxAntennaPatternArray_5_tab = readtable(filenameTxGainPattern, optsReadTableTxPattern, "UseExcel", false);
            TxAntennaPatternArray_5     = table2array(TxAntennaPatternArray_5_tab);
            TxPattern_gain_5  = TxAntennaPatternArray_5(2:length(TxPattern_phi_temp_0), 2:length(TxPattern_theta_temp_0));
            TxPattern_gain_5_phiSort                   = TxPattern_gain_5(indexPhiSort,:);
            TxPattern_gain_5_phiSort_complete(1,:) = TxPattern_gain_5_phiSort(length(TxPattern_phi_180_sort),:);
            TxPattern_gain_5_phiSort_complete(2:length(TxPattern_phi_180_sort)+1,:) = TxPattern_gain_5_phiSort;
            TxPattern_gain_res(:,:,2) = TxPattern_gain_5_phiSort_complete(:,indexResizedTheta);
        end

        %NOTE: pattern is converted to linear units!
        TxPattern_gain_linear_res = db2pow(TxPattern_gain_res);

        TRAp.theta_Txpattern  = TxPattern_theta_res;
        TRAp.phi_Txpattern    = TxPattern_phi_180_sort_complete';
        TRAp.GainTx_linear    = TxPattern_gain_linear_res;

        if tag.plotPatterns == 1
                figure
                imagesc(TxPattern_theta_res, TxPattern_phi_180_sort_complete, TxPattern_gain_res(:,:,1))
                ylabel('\phi'), xlabel('\theta'), colorbar, title('Original Tx pattern - L1')
                %caxis([max(max(TxPattern_gain_temp(:,:,1))) - 40, max(max(TxPattern_gain_temp(:,:,1)))])
        end

        if tag.PattInterp == 1

            [TxPattern_theta_res_grid, TxPattern_phi_sort_grid] = meshgrid(TxPattern_theta_res, TxPattern_phi_180_sort_complete);

            TxPattern_theta_interp  =  abs(-90:configParam.dang:0);
            TxPattern_phi_interp         = -180:configParam.dang:180;
            [TxPattern_theta_interp_grid, TxPattern_phi_interp_grid] = meshgrid(TxPattern_theta_interp, TxPattern_phi_interp);

            TxPattern_gain_res_interp = interp2(TxPattern_theta_res_grid, TxPattern_phi_sort_grid,...
                        TxPattern_gain_res(:,:,1), TxPattern_theta_interp_grid, TxPattern_phi_interp_grid,'makima');
            if geoSYSp.indexTxSignal == 1
                TxPattern_gain_res_interp(:,:,2) = interp2(TxPattern_theta_res_grid, TxPattern_phi_sort_grid,...
                        TxPattern_gain_res(:,:,2), TxPattern_theta_interp_grid, TxPattern_phi_interp_grid,'makima');
            end

            %NOTE: pattern is converted to linear units!
            TxPattern_gain_linear_interp = db2pow(TxPattern_gain_res_interp);

            %overwrite to interp
            TRAp.theta_Txpattern  = TxPattern_theta_interp;
            TRAp.phi_Txpattern    = TxPattern_phi_interp';
            TRAp.GainTx_linear    = TxPattern_gain_linear_interp;

            if tag.plotPatterns == 1
                figure
                imagesc(TxPattern_theta_interp, TxPattern_phi_interp, TxPattern_gain_res_interp(:,:,1))
                ylabel('\phi'), xlabel('\theta'), colorbar, title('Interp Tx pattern  - L1')
                %caxis([max(max(TxPattern_gain_temp_new(:,:,1))) - 40, max(max(TxPattern_gain_temp_new(:,:,1)))])
            end

        end

     elseif geoSYSp.indexTxSignal == 2 || geoSYSp.indexTxSignal == 3

        %read Galileo E1 pattern
        optsReadTableTxPattern.Sheet = 'Galileo_E1';    
        TxAntennaPatternArray_tab = readtable(filenameTxGainPattern, optsReadTableTxPattern, "UseExcel", false);
        TxAntennaPatternArray     = table2array(TxAntennaPatternArray_tab);

        %reshape the variable in order to have (as it is for the Rx) the gain in a vector 720x91,
        %phi in a vector 720x1, theta in a vector 91x1
        TxPattern_phi_temp_0   = TxAntennaPatternArray(:,1);
        TxPattern_phi          = TxPattern_phi_temp_0(2:length(TxPattern_phi_temp_0));       
        TxPattern_theta_temp_0 = TxAntennaPatternArray(1,:);
        TxPattern_theta   = TxPattern_theta_temp_0(2:length(TxPattern_theta_temp_0));
        TxPattern_gain_orig         = TxAntennaPatternArray(2:length(TxPattern_phi_temp_0), 2:length(TxPattern_theta_temp_0));

        %make the format of the pattern similar to that of Rx (to be
        %coherent with the used frame)
        %wrap phi to 180° and and sort to have phi ranging from -180° to 180°
        TxPattern_phi_180                     = wrapTo180(TxPattern_phi);
        [TxPattern_phi_180_sort,indexPhiSort] = sort(TxPattern_phi_180);
        TxPattern_gain_phiSort                   = TxPattern_gain_orig(indexPhiSort,:);
        %add a row in phi and pattern corresponding to -180° (equal to the
        %row corresponding to 180°) to improve the performance of the
        %interpolation
        TxPattern_phi_180_sort_complete(1)      = -180.0;
        TxPattern_phi_180_sort_complete(2:length(TxPattern_phi_180_sort)+1) = TxPattern_phi_180_sort;
        TxPattern_gain_phiSort_complete(1,:) = TxPattern_gain_phiSort(length(TxPattern_phi_180_sort),:);
        TxPattern_gain_phiSort_complete(2:length(TxPattern_phi_180_sort)+1,:) = TxPattern_gain_phiSort;
        TxPattern_gain=TxPattern_gain_phiSort_complete;

        if geoSYSp.indexTxSignal == 3
            %read Galileo E5 pattern and resize to the same size of E1
            optsReadTableTxPattern.Sheet = 'Galileo_E5';    
            TxAntennaPatternArray_5_tab = readtable(filenameTxGainPattern, optsReadTableTxPattern, "UseExcel", false);
            TxAntennaPatternArray_5     = table2array(TxAntennaPatternArray_5_tab);
            TxPattern_gain_5 = TxAntennaPatternArray_5(2:length(TxPattern_phi_temp_0), 2:length(TxPattern_theta_temp_0));
            TxPattern_gain_5_phiSort                   = TxPattern_gain_5(indexPhiSort,:);
            TxPattern_gain_5_phiSort_complete(1,:) = TxPattern_gain_5_phiSort(length(TxPattern_phi_180_sort),:);
            TxPattern_gain_5_phiSort_complete(2:length(TxPattern_phi_180_sort)+1,:) = TxPattern_gain_5_phiSort;
            TxPattern_gain(:,:,2) = TxPattern_gain_5_phiSort_complete;
        end

        %NOTE: pattern is converted to linear units!
        TxPattern_gain_linear = db2pow(TxPattern_gain);
        TRAp.theta_Txpattern  = fliplr(TxPattern_theta);
        TRAp.phi_Txpattern    = TxPattern_phi_180_sort_complete';
        TRAp.GainTx_linear    = fliplr(TxPattern_gain_linear);

        if tag.plotPatterns == 1
             figure
             imagesc(TxPattern_theta, TxPattern_phi_180_sort_complete, TxPattern_gain(:,:,1))
             ylabel('\phi'), xlabel('\theta'), colorbar, title('Original Tx pattern - E1')
             %caxis([max(TxPattern_gain_temp(:)) - 40, max(TxPattern_gain_temp(:))])
        end

        if tag.PattInterp == 1
            
            [TxPattern_theta_grid, TxPattern_phi_sort_grid] = meshgrid(TxPattern_theta, TxPattern_phi_180_sort_complete);

            TxPattern_theta_interp  =  0:configParam.dang:45;
            TxPattern_phi_interp         = -180:configParam.dang:180;
            [TxPattern_theta_interp_grid, TxPattern_phi_interp_grid] = meshgrid(TxPattern_theta_interp, TxPattern_phi_interp);

            TxPattern_gain_interp = interp2(TxPattern_theta_grid, TxPattern_phi_sort_grid,...
                        TxPattern_gain(:,:,1), TxPattern_theta_interp_grid, TxPattern_phi_interp_grid,'makima');
            if geoSYSp.indexTxSignal == 3
                TxPattern_gain_interp(:,:,2) = interp2(TxPattern_theta_grid, TxPattern_phi_sort_grid,...
                        TxPattern_gain(:,:,2), TxPattern_theta_interp_grid, TxPattern_phi_interp_grid,'makima');
            end

            %NOTE: pattern is converted to linear units!
            TxPattern_gain_linear_interp = db2pow(TxPattern_gain_interp);

            %overwrite to interp
            TRAp.theta_Txpattern  = fliplr(TxPattern_theta_interp);
            TRAp.phi_Txpattern    = TxPattern_phi_interp';
            TRAp.GainTx_linear    = fliplr(TxPattern_gain_linear_interp);             

            if tag.plotPatterns == 1
                figure
                imagesc(TxPattern_theta_interp, TxPattern_phi_interp, TxPattern_gain_interp(:,:,1))
                ylabel('\phi'), xlabel('\theta'), colorbar, title('Interp Tx pattern - E1')
                %caxis([max(TxPattern_gain_temp_new(:)) - 40, max(TxPattern_gain_temp_new(:))]
            end
        end
     end

    %% --- read the land cover map and define the classes parameters ---------------

    if tag.varibleBioGeoInputs==0
        [landCoverMap,coverMap_IDs,terrainParameters,vegetationParameters,landCoverMap_all] = readCoverClassesAndDefineParam(tag,filenameCoverMap, SPlatlonAlt, DEMinfo, bioGeoInputs, presentLogInputsFilename);
    else
        [landCoverMap,coverMap_IDs,terrainParameters,vegetationParameters,descreteVariableParam,randomVariableParameter,landCoverMap_all] = readCoverClassesAndDefineVariableParam(tag,filenameCoverMap, SPlatlonAlt, DEMinfo,tag.scenario, bioGeoInputsVariable, presentLogInputsFilename);
    end

    if tag.saveCoverMap==1
        save([dest_dir 'landCoverMap.mat'], 'landCoverMap');
    end


    %% ------ START THE SIMULATION-------------------
    if tag.varibleBioGeoInputs       == 0
        %simulations for a fixed value of input bio and geo variables
        [outputSimulations, infoAtSPforRefFile, UP_pattern] = HSAVERS_IHS_OENUonDEM(dest_dir, outputFilename, Tx, VTx, Rx, VRx, geoSYSp, TRAp,...
                                                                                    DEMinfo, landCoverMap, coverMap_IDs, instrumentConfParam, terrainParameters,...
                                                                                    vegetationParameters, ionisationLevelParam, GNSS_Week, GNSS_SoW, tag,landCoverMap_all);
        
        SPlatlonAlt_HSAVERS=[infoAtSPforRefFile.HSAVERS_SPlat' infoAtSPforRefFile.HSAVERS_SPlon' infoAtSPforRefFile.HSAVERS_SPalt_5km'];
        TrackID_cohChannel_0 = uint32(0);
        for i = 1:length(TrackID)
            TrackID_cohChannel_i = uint32(zeros(1,geoSYSp.incoIntTime/geoSYSp.cohIntTime))+uint32(TrackID(i));
            TrackID_cohChannel_0 = [TrackID_cohChannel_0 TrackID_cohChannel_i];% only valid in trackwise mode
            TrackID_cohChannel =TrackID_cohChannel_0(2:length(TrackID_cohChannel_0));
        end

        %create structure to take track of all inputs and outputs of each run
        inOutReferenceFile.TrackID=uint32(TrackID);
        inOutReferenceFile.orbitParameters.GNSS_Week = GNSS_Week;
        inOutReferenceFile.orbitParameters.GNSS_SoW  = GNSS_SoW;
        inOutReferenceFile.orbitParameters.PRN       = PRN;
        inOutReferenceFile.orbitParameters.SVN       = SVN;
        inOutReferenceFile.orbitParameters.Rx        = Rx;
        inOutReferenceFile.orbitParameters.VRx       = VRx;
        inOutReferenceFile.orbitParameters.Tx        = Tx;
        inOutReferenceFile.orbitParameters.VTx       = VTx;

        inOutReferenceFile.inputFilenames.filenameOrbits            = orbitdataHG1;
        inOutReferenceFile.inputFilenames.filenameCoverMap          = filenameCoverMap;
        inOutReferenceFile.inputFilenames.filenameTxGainPattern     = TxGainPattern;
        if satellite == "HydroGNSS-1"  % added by ansha
            inOutReferenceFile.inputFilenames.filenameRxDOWNGainPattern = nadirantennadataHG1;
            inOutReferenceFile.inputFilenames.filenameRxUPGainPattern   = zenithantennadataHG1;
        else
            inOutReferenceFile.inputFilenames.filenameRxDOWNGainPattern = nadirantennadataHG2;
            inOutReferenceFile.inputFilenames.filenameRxUPGainPattern   = zenithantennadataHG2;
        end
        inOutReferenceFile.outputDirectory                          = dest_dir;
        inOutReferenceFile.outputFilename                           = outputFilename;
        inOutReferenceFile.geoSYSp                                  = geoSYSp;
        inOutReferenceFile.geoSYSp.SPlat_series                     = infoAtSPforRefFile.HSAVERS_SPlat;
        inOutReferenceFile.geoSYSp.SPlon_series                     = infoAtSPforRefFile.HSAVERS_SPlon;
        inOutReferenceFile.geoSYSp.SPalt_5km_series                 = infoAtSPforRefFile.HSAVERS_SPalt_5km;
        inOutReferenceFile.geoSYSp.SPalt_origDEM_series             = infoAtSPforRefFile.HSAVERS_SPalt_origDEM;
        inOutReferenceFile.geoSYSp.TxBeamwidth                      = TRAp.TxBeamwidth;
        inOutReferenceFile.geoSYSp.RxBeamwidth                      = TRAp.RxBeamwidth;
        inOutReferenceFile.geoSYSp.Mode                             = mode1;
        inOutReferenceFile.DEMinfo                                  = DEMinfo;
        inOutReferenceFile.TxSettings.TxSignalPower                 = TRAp.TxPower_dBW;
        inOutReferenceFile.TxSettings.TxSignalPowerError            = TRAp.TxPowerError_dBW;
        inOutReferenceFile.TxSettings.TxAttitudeAngles_nominal      = TRAp.Tx_nomAtt_angles;
        inOutReferenceFile.TxSettings.TxAttitudeAngles_error        = TRAp.TxAttitude_angles_GUI;
        inOutReferenceFile.TxSettings.TxEuler_angles                = TRAp.TxEuler_angles;
        inOutReferenceFile.TxSettings.TxFrequency                   = TRAp.TxFrequency;
        inOutReferenceFile.TxSettings.TxMismatch                    = buttonTxMismatch;

        inOutReferenceFile.TxSettings.Yaw_TX                        = infoAtSPforRefFile.Yaw_TX;
        % with Error
        inOutReferenceFile.TxSettings.Yaw_TX_Error                        = infoAtSPforRefFile.Yaw_TX_Error;

        
        inOutReferenceFile.RxSettings.RxAttitudeAngles_error        = TRAp.RxAttitude_angles_GUI;
        inOutReferenceFile.RxSettings.RxDOWNmismatch    = buttonRxDOWNmismatch;
        inOutReferenceFile.RxSettings.RxUPmismatch      = buttonRxUPmismatch;
        inOutReferenceFile.RxSettings.RxEuler_angles    = TRAp.RxEuler_angles;

        if geoSYSp.Rx_Attitude_Variation_RPY(1) == 1 || geoSYSp.Rx_Attitude_Variation_RPY(2) == 1 || geoSYSp.Rx_Attitude_Variation_RPY(3) == 1
            inOutReferenceFile.RxSettings.Rx_roll = geoSYSp.Rx_roll;
            inOutReferenceFile.RxSettings.Rx_pitch = geoSYSp.Rx_pitch;
            inOutReferenceFile.RxSettings.Rx_yaw = geoSYSp.Rx_yaw;
            inOutReferenceFile.RxSettings.Rx_Attitude_Variation_RPY = geoSYSp.Rx_Attitude_Variation_RPY;
            inOutReferenceFile.RxSettings.Rx_Attitude_Variation_mode = 'Attitude variation from track file';
            clear inOutReferenceFile.geoSYSp.Rx_roll inOutReferenceFile.geoSYSp.Rx_pitch inOutReferenceFile.geoSYSp.Rx_yaw inOutReferenceFile.geoSYSp.Rx_Attitude_Variation_RPY
        else
            inOutReferenceFile.RxSettings.RxAttitudeAngles_nominal      = TRAp.Rx_nomAtt_angles;
            inOutReferenceFile.RxSettings.Rx_Attitude_Variation_RPY = geoSYSp.Rx_Attitude_Variation_RPY;
            inOutReferenceFile.RxSettings.Rx_Attitude_Variation_mode = 'Stable attitude from configuration file';
            clear inOutReferenceFile.geoSYSp.Rx_roll inOutReferenceFile.geoSYSp.Rx_pitch inOutReferenceFile.geoSYSp.Rx_yaw inOutReferenceFile.geoSYSp.Rx_Attitude_Variation_RPY
        
        end
        inOutReferenceFile.terrainParameters            = terrainParameters;
        inOutReferenceFile.vegetationParameters         = vegetationParameters;

        %look for the cover class where the SP falls and the bio/geo
        %parameters given as input for the SP class
        bioGeoParametersAtSP=biogeoparamAtSP(infoAtSPforRefFile,coverMap_IDs,terrainParameters,vegetationParameters);

        inOutReferenceFile.bioGeoParametersAtSP         = bioGeoParametersAtSP;

        inOutReferenceFile.Tx1GainAtSP                   = infoAtSPforRefFile.Tx1GainAtSP;
        inOutReferenceFile.Tx1GainAtRx                   = infoAtSPforRefFile.Tx1GainAtRx;
        inOutReferenceFile.RxGainLR1atSP                 = infoAtSPforRefFile.RxGainLR1atSP;
        inOutReferenceFile.RxGainRR1atSP                 = infoAtSPforRefFile.RxGainRR1atSP;

        inOutReferenceFile.squaredDistanceSum            = infoAtSPforRefFile.squaredDistanceSum;
        inOutReferenceFile.RxZenithGainTowardsTx1        = UP_pattern(:,1);
        inOutReferenceFile.absoluteDelayAtSP_1_microsec             = outputSimulations.absoluteDelayAtSP_1_microsec;
        inOutReferenceFile.delayAtSPReferencedToDirect_1_microsec   = outputSimulations.delayAtSPReferencedToDirect_1_microsec;

        inOutReferenceFile.max_sigPower_noNoise_LR1         = outputSimulations.max_sigPower_noNoise_LR1;
        inOutReferenceFile.max_sigPower_noNoise_RR1         = outputSimulations.max_sigPower_noNoise_RR1;
        inOutReferenceFile.max_sigPower_Noise_LR1           = outputSimulations.max_sigPower_Noise_LR1;
        inOutReferenceFile.max_sigPower_Noise_RR1           = outputSimulations.max_sigPower_Noise_RR1;
        inOutReferenceFile.max_reflectivity_noNoise_LR1_dB  = outputSimulations.max_reflectivity_noNoise_LR1_dB;
        inOutReferenceFile.max_reflectivity_noNoise_RR1_dB  = outputSimulations.max_reflectivity_noNoise_RR1_dB;
        inOutReferenceFile.max_reflectivity_NoisyCal_LR1_dB = outputSimulations.max_reflectivity_Noisy_LR1_dB;
        inOutReferenceFile.max_reflectivity_NoisyCal_RR1_dB = outputSimulations.max_reflectivity_Noisy_RR1_dB;
        inOutReferenceFile.max_geobReflec_noNoise_LR1_dB    = outputSimulations.max_geobReflec_noNoise_LR1_dB;
        inOutReferenceFile.max_geobReflec_noNoise_RR1_dB    = outputSimulations.max_geobReflec_noNoise_RR1_dB;
        inOutReferenceFile.max_geobPower_Noisy_LR1   = outputSimulations.max_geobPower_Noisy_LR1;
        inOutReferenceFile.max_geobPower_Noisy_RR1   = outputSimulations.max_geobPower_Noisy_RR1;

%         inOutReferenceFile.IdealAntenna = outputSimulations.IdealAntenna_all;
%         inOutReferenceFile.MismatchAntenna = outputSimulations.MismatchAntenna_all;

% inOutReferenceFile.reflcomponent.max_sigPower_noNoise_LR      = outputSimulations.reflcomponentx.max_sigPower_noNoise_LR;
% inOutReferenceFile.reflcomponent.max_sigPower_Noise_LR        = outputSimulations.reflcomponentx.max_sigPower_Noise_LR;
% inOutReferenceFile.reflcomponent.TxGainAtRx                   = outputSimulations.reflcomponentx.TxGainAtRx;
% inOutReferenceFile.reflcomponent.PonSurf                      = outputSimulations.reflcomponentx.PonSurf;
% inOutReferenceFile.reflcomponent.LR_UP_Power                  = outputSimulations.reflcomponentx.LR_UP_Power;
% inOutReferenceFile.reflcomponent.Gain_LHCP_onsurf_lin         = outputSimulations.reflcomponentx.Gain_LHCP_onsurf_lin;
% inOutReferenceFile.reflcomponent.PATH_factor                  = outputSimulations.reflcomponentx.PATH_factor;
% 
% % New reflectivity components
% % inOutReferenceFile.reflcomponent.sigPower_noNoise_LR_lin      = outputSimulations.reflcomponentx.sigPower_noNoise_LR_lin;
% inOutReferenceFile.reflcomponent.fourPi_squared               = outputSimulations.reflcomponentx.fourPi_squared;
% inOutReferenceFile.reflcomponent.squaredDistanceSum           = outputSimulations.reflcomponentx.squaredDistanceSum;
% inOutReferenceFile.reflcomponent.lambda_squared               = outputSimulations.reflcomponentx.lambda_squared;
% inOutReferenceFile.reflcomponent.Gain_LHCP_onsurf_atSP        = outputSimulations.reflcomponentx.Gain_LHCP_onsurf_atSP;
% inOutReferenceFile.reflcomponent.TxSignalPower_lin            = outputSimulations.reflcomponentx.TxSignalPower_lin;
% inOutReferenceFile.reflcomponent.TxGainAtSP                   = outputSimulations.reflcomponentx.TxGainAtSP;
% inOutReferenceFile.reflcomponent.sigPower_noNoise_LR_lin      = outputSimulations.reflcomponentx.sigPower_noNoise_LR_lin;
% inOutReferenceFile.reflcomponent.reflectivity_noNoise_LR_tot      = outputSimulations.reflcomponentx.reflectivity_noNoise_LR_tot;


% inOutReferenceFile.reflcomponent.reflectivity_noNoise_LR_tot  = outputSimulations.reflcomponentx.reflectivity_noNoise_LR_tot;


       %%Amir
%         inOutReferenceFile.AMIR_Sig_postI = outputSimulations.Sig_postI;
%         inOutReferenceFile.AMIR_Noise_postI_cnt = outputSimulations.Noise_postI_cnt;
%         inOutReferenceFile.AMIR_Sum_Sig_Noise_int_amir = outputSimulations.Sum_Sig_Noise_int_amir;
%         inOutReferenceFile.AMIR_RR_UP_Power_Flu_def = outputSimulations.RR_UP_Power_Flu_def;
%         inOutReferenceFile.AMIR_NoiseSamples_def = outputSimulations.NoiseSamples_def;
%         inOutReferenceFile.AMIR_RR_UP_Power_def = outputSimulations.RR_UP_Power_def;

       %%Amir
        if tag.SNR==1
            inOutReferenceFile.SNR_LR1                   = outputSimulations.SNR_LR1;
            inOutReferenceFile.SNR_RR1                   = outputSimulations.SNR_RR1;

%             inOutReferenceFile.geobSNR_LR1               = outputSimulations.geobSNR_LR1;
%             inOutReferenceFile.geobSNR_RR1               = outputSimulations.geobSNR_RR1;
        end
        if geoSYSp.indexTxSignal==1 || geoSYSp.indexTxSignal==3
            inOutReferenceFile.Tx5GainAtSP                   = infoAtSPforRefFile.Tx5GainAtSP;
            inOutReferenceFile.Tx5GainAtRx                   = infoAtSPforRefFile.Tx5GainAtRx;
            inOutReferenceFile.RxGainLR5atSP                 = infoAtSPforRefFile.RxGainLR5atSP;
            inOutReferenceFile.RxGainRR5atSP                 = infoAtSPforRefFile.RxGainRR5atSP;
            inOutReferenceFile.RxZenithGainTowardsTx5        = UP_pattern(:,2);
            inOutReferenceFile.absoluteDelayAtSP_5_microsec  = outputSimulations.absoluteDelayAtSP_5_microsec;
            inOutReferenceFile.delayAtSPReferencedToDirect_5_microsec   = outputSimulations.delayAtSPReferencedToDirect_5_microsec;
            inOutReferenceFile.max_sigPower_noNoise_LR5         = outputSimulations.max_sigPower_noNoise_LR5;
            inOutReferenceFile.max_sigPower_noNoise_RR5         = outputSimulations.max_sigPower_noNoise_RR5;
            inOutReferenceFile.max_sigPower_Noise_LR5           = outputSimulations.max_sigPower_Noise_LR5;
            inOutReferenceFile.max_sigPower_Noise_RR5           = outputSimulations.max_sigPower_Noise_RR5;
            inOutReferenceFile.max_reflectivity_noNoise_LR5_dB  = outputSimulations.max_reflectivity_noNoise_LR5_dB;
            inOutReferenceFile.max_reflectivity_noNoise_RR5_dB  = outputSimulations.max_reflectivity_noNoise_RR5_dB;
            inOutReferenceFile.max_reflectivity_NoisyCal_LR5_dB = outputSimulations.max_reflectivity_Noisy_LR5_dB;
            inOutReferenceFile.max_reflectivity_NoisyCal_RR5_dB = outputSimulations.max_reflectivity_Noisy_RR5_dB;
            inOutReferenceFile.max_geobReflec_noNoise_LR5_dB    = outputSimulations.max_geobReflec_noNoise_LR5_dB;
            inOutReferenceFile.max_geobReflec_noNoise_RR5_dB    = outputSimulations.max_geobReflec_noNoise_RR5_dB;
            inOutReferenceFile.max_geobPower_Noisy_LR5   = outputSimulations.max_geobPower_Noisy_LR5;
            inOutReferenceFile.max_geobPower_Noisy_RR5   = outputSimulations.max_geobPower_Noisy_RR5;
            
            if tag.SNR==1
                inOutReferenceFile.SNR_LR5                   = outputSimulations.SNR_LR5;
                inOutReferenceFile.SNR_RR5                   = outputSimulations.SNR_RR5;

%                 inOutReferenceFile.geobSNR_LR5               = outputSimulations.geobSNR_LR5;
%                 inOutReferenceFile.geobSNR_RR5               = outputSimulations.geobSNR_RR5;
            end
        end
    else
        %variable that varies in a descrete way (soil moisture for scenario of
        %soil moisture and biomass for scenario forest). Note that all the
        %other bio/geo inputs will vary randomly
        varPam             = descreteVariableParam(1):descreteVariableParam(2):descreteVariableParam(3);
        %define the output variables (this should be large enougth to
        %include all outputs for all runs
        inOutReferenceFile = struct([]);
        TrackID            = uint32(TrackID);
        numSP              = length(TrackID);
        TrackID_cohChannel_0 = uint32(0);
        for i = 1:length(TrackID)
            TrackID_cohChannel_i = uint32(zeros(1,geoSYSp.incoIntTime/geoSYSp.cohIntTime))+uint32(TrackID(i));
            TrackID_cohChannel_0 = [TrackID_cohChannel_0 TrackID_cohChannel_i];% only valid in trackwise mode
            TrackID_cohChannel =TrackID_cohChannel_0(2:length(TrackID_cohChannel_0));
        end

        outputSimulations.max_sigPower_noNoise_LR1        = zeros(1,numSP*length(varPam)*descreteVariableParam(4));
        outputSimulations.max_sigPower_noNoise_RR1        = outputSimulations.max_sigPower_noNoise_LR1;
        outputSimulations.max_sigPower_Noise_LR1          = outputSimulations.max_sigPower_noNoise_LR1;
        outputSimulations.max_sigPower_Noise_RR1          = outputSimulations.max_sigPower_noNoise_LR1;
        
        outputSimulations.delayAtSPReferencedToDirect_1_microsec = outputSimulations.max_sigPower_noNoise_LR1;
        outputSimulations.Atotal_1                        = zeros(1,numSP*length(varPam)*descreteVariableParam(4)*geoSYSp.incoIntTime/geoSYSp.cohIntTime);        

        outputSimulations.LR1_UP_Power                    = zeros(numSP*length(varPam)*descreteVariableParam(4),1);
        outputSimulations.RR1_UP_Power                    = outputSimulations.LR1_UP_Power;
        outputSimulations.LR1_UP_Power_Fluc               = outputSimulations.LR1_UP_Power;
        outputSimulations.RR1_UP_Power_Fluc               = outputSimulations.LR1_UP_Power;
        outputSimulations.PN_star_1                       = outputSimulations.LR1_UP_Power;
        outputSimulations.sigPower_noNoise_LR1_lin        = zeros(geoSYSp.ny_Dfreq,geoSYSp.nx_Rdelay,numSP*length(varPam)*descreteVariableParam(4));
        outputSimulations.sigPower_noNoise_RR1_lin        = outputSimulations.sigPower_noNoise_LR1_lin;
        outputSimulations.sigPower_Noise_LR1_lin          = outputSimulations.sigPower_noNoise_LR1_lin;
        outputSimulations.sigPower_Noise_RR1_lin          = outputSimulations.sigPower_noNoise_LR1_lin;


        outputSimulations.NLblackbody_DOWN_LR_1              = outputSimulations.sigPower_noNoise_LR1_lin;
        outputSimulations.NLblackbody_DOWN_RR_1              = outputSimulations.sigPower_noNoise_LR1_lin;
        outputSimulations.NLblackbody_UP_1                = outputSimulations.sigPower_noNoise_LR1_lin;
        outputSimulations.NLblackbody_DOWNmean_LR_1          = outputSimulations.max_sigPower_noNoise_LR1;
        outputSimulations.NLblackbody_DOWNmax_LR_1           = outputSimulations.max_sigPower_noNoise_LR1;
        outputSimulations.NLblackbody_DOWNmean_RR_1          = outputSimulations.max_sigPower_noNoise_LR1;
        outputSimulations.NLblackbody_DOWNmax_RR_1           = outputSimulations.max_sigPower_noNoise_LR1;
        outputSimulations.NLblackbody_UPmean_1            = outputSimulations.max_sigPower_noNoise_LR1;
        outputSimulations.NLblackbody_UPmax_1             = outputSimulations.max_sigPower_noNoise_LR1;

        outputSimulations.max_sigPower_noNoise_LR5        = zeros(1,numSP*length(varPam)*descreteVariableParam(4));
        outputSimulations.max_sigPower_noNoise_RR5        = outputSimulations.max_sigPower_noNoise_LR5;
        outputSimulations.max_sigPower_Noise_LR5          = outputSimulations.max_sigPower_noNoise_LR5;
        outputSimulations.max_sigPower_Noise_RR5          = outputSimulations.max_sigPower_noNoise_LR5;
        
        outputSimulations.delayAtSPReferencedToDirect_5_microsec = outputSimulations.max_sigPower_noNoise_LR5;
        outputSimulations.Atotal_5                        = zeros(1,numSP*length(varPam)*descreteVariableParam(4)*geoSYSp.incoIntTime/geoSYSp.cohIntTime);        

        outputSimulations.LR5_UP_Power                    = zeros(numSP*length(varPam)*descreteVariableParam(4),1);
        outputSimulations.RR5_UP_Power                    = outputSimulations.LR5_UP_Power;
        outputSimulations.LR5_UP_Power_Fluc               = outputSimulations.LR5_UP_Power;
        outputSimulations.RR5_UP_Power_Fluc               = outputSimulations.LR5_UP_Power;
        outputSimulations.PN_star_5                       = outputSimulations.LR5_UP_Power;
        outputSimulations.sigPower_noNoise_LR5_lin        = zeros(geoSYSp.ny_Dfreq,geoSYSp.nx_Rdelay,numSP*length(varPam)*descreteVariableParam(4));
        outputSimulations.sigPower_noNoise_RR5_lin        = outputSimulations.sigPower_noNoise_LR5_lin;
        outputSimulations.sigPower_Noise_LR5_lin          = outputSimulations.sigPower_noNoise_LR5_lin;
        outputSimulations.sigPower_Noise_RR5_lin          = outputSimulations.sigPower_noNoise_LR5_lin;

        outputSimulations.NLblackbody_DOWN_LR_5           = outputSimulations.sigPower_noNoise_LR5_lin;
        outputSimulations.NLblackbody_DOWN_RR_5           = outputSimulations.sigPower_noNoise_LR5_lin;
        outputSimulations.NLblackbody_UP_5                = outputSimulations.sigPower_noNoise_LR5_lin;
        outputSimulations.NLblackbody_DOWNmean_LR_5       = outputSimulations.max_sigPower_noNoise_LR5;
        outputSimulations.NLblackbody_DOWNmax_LR_5        = outputSimulations.max_sigPower_noNoise_LR5;
        outputSimulations.NLblackbody_DOWNmean_RR_5       = outputSimulations.max_sigPower_noNoise_LR5;
        outputSimulations.NLblackbody_DOWNmax_RR_5        = outputSimulations.max_sigPower_noNoise_LR5;
        outputSimulations.NLblackbody_UPmean_5            = outputSimulations.max_sigPower_noNoise_LR5;
        outputSimulations.NLblackbody_UPmax_5             = outputSimulations.max_sigPower_noNoise_LR5;
        

        %repeat the run for all descrete values and for the number of
        %times given as input
        for m=1:length(varPam)
            for l=1:descreteVariableParam(4)
                %define the bio/geo parameters that will be used for each
                %run
                switch tag.scenario
                    case 1 %soil moisture
                        indexNoWater=find(terrainParameters(:,1)~=2);
                        %descrete variation of the soil moisture (this is
                        %homogeneous over the whole simulation area)
                        terrainParameters(indexNoWater,5)=varPam(m);
                        %random variation of the biomass for the forest classes
                        indexForest=find(vegetationParameters(:,2)==1);
                            a=randomVariableParameter(indexForest,3) ; 
                            b=randomVariableParameter(indexForest,4) ; 
                        vegetationParameters(indexForest,7)=random('Normal',(b+a)/2,(b-a)/configParam.biogeofact) ;
                        disp(['Soil moisture [%]: ', num2str(varPam(m)), ', RUN ' , num2str(l)])

                    case 2 %forest
                        indexNoWater=find(terrainParameters(:,1)~=2);
                        %random variation of the soil moisture (each class has
                        %a different soil moisture)
                            a=randomVariableParameter(indexNoWater,3) ;
                            b=randomVariableParameter(indexNoWater,4) ;
                        terrainParameters(indexNoWater,5)=random('Normal',(b+a)/2,(b-a)/configParam.biogeofact);
                        %descrete variation of the biomass (all forest classes
                        %will have the same biomass)
                        indexForest=find(vegetationParameters(:,2)==1);
                        vegetationParameters(indexForest,7)=varPam(m);
                        disp(['Biomass [t/ha]: ', num2str(varPam(m)) , ', RUN ' , num2str(l)])

                    case 3 %freeze-thaw
                        terrainParameters(:,1) = varPam(m); %0 all classes have unfrozen soil, 1 all classes have frozen soil
                        %if frozen, soil temperature must be below zero and
                        %water must have soil moisture 10%
                        indexWater=find(terrainParameters(:,1)==2);
                        if varPam(m)==1
                            terrainParameters(:,7) = -10.0;
                            terrainParameters(indexWater,5)=10.0;
                        end
                        indexNoWater=find(terrainParameters(:,1)~=2);
                        %random variation of the soil moisture (each class has
                        %a different soil moisture)
                            a=randomVariableParameter(indexNoWater,3) ; 
                            b=randomVariableParameter(indexNoWater,4) ;                       
                        terrainParameters(indexNoWater,5)=random('Normal',(b+a)/2,(b-a)/configParam.biogeofact);                       
                        %random variation of the biomass for the forest classes
                        indexForest=find(vegetationParameters(:,2)==1);
                            a=randomVariableParameter(indexForest,7) ;
                            b=randomVariableParameter(indexForest,8) ;
                        vegetationParameters(indexForest,7)=random('Normal',(b+a)/2,(b-a)/configParam.biogeofact);
                        disp(['Frozen [0 = No; 1 = Yes] ', num2str(varPam(m)) , ', RUN ' , num2str(l)])

                    case 4 %wetland
                        %only vegetated classes defined flooded in the CCI map
                        indexFloodedVegetation=find(coverMap_IDs(:,1)==5 | coverMap_IDs(:,1)==6);
                        terrainParameters(indexFloodedVegetation,1) = varPam(m); % 0 not flooded, 2 flooded
                        indexNoWater=find(terrainParameters(:,1)~=2);
                        %random variation of the soil moisture (each class has
                        %a different soil moisture)
                            a=randomVariableParameter(indexNoWater,3) ;
                            b=randomVariableParameter(indexNoWater,4) ;
                        terrainParameters(indexNoWater,5)=random('Normal',(b+a)/2,(b-a)/configParam.biogeofact);
                        indexWater=find(terrainParameters(:,1)==2);
                        terrainParameters(indexWater,5)=80.0;
                        terrainParameters(indexWater,2)=0.1;
                        %random variation of the biomass for the forest classes
                        indexForest=find(vegetationParameters(:,2)==1);
                            a=randomVariableParameter(indexForest,7) ;
                            b=randomVariableParameter(indexForest,8) ;
                        vegetationParameters(indexForest,7)=random('Normal',(b+a)/2,(b-a)/configParam.biogeofact);
                        disp(['Flooded [0 = No; 2 = Yes] ', num2str(varPam(m)) , ', RUN ' , num2str(l)])
                     
                end
                indexNoWater=find(terrainParameters(:,1)~=2);
                %random soil roughness
                    a=randomVariableParameter(indexNoWater,1) ;
                    b=randomVariableParameter(indexNoWater,2) ;
                terrainParameters(indexNoWater,2)   = random('Normal',(b+a)/2,(b-a)/configParam.biogeofact);
                %random LAI (for the low vegetation class)
                indexLowVeg                         = find(vegetationParameters(:,1)==1 & vegetationParameters(:,2)==0);
                    a=randomVariableParameter(indexLowVeg,5) ; 
                    b=randomVariableParameter(indexLowVeg,6) ;            
                vegetationParameters(indexLowVeg,3) =random('Normal',(b+a)/2,(b-a)/configParam.biogeofact);

                %output filename that will be used for each run
                outputFilename_variable      = cell(1,numSP);
                outputFilename_variable(1,:) = {['_' num2str(m) num2str(100000*l)]};
                outputFilename_lm            = cellstr([char(outputFilename) char(outputFilename_variable)]);

                [outputSimulations_lm, infoAtSPforRefFile_lm, UP_pattern_lm] = HSAVERS_IHS_OENUonDEM(dest_dir,outputFilename_lm,Tx, VTx, Rx, VRx, geoSYSp, TRAp,...
                                                                                                        DEMinfo,landCoverMap,coverMap_IDs, instrumentConfParam,terrainParameters,...
                                                                                                        vegetationParameters, ionisationLevelParam, GNSS_Week, GNSS_SoW, tag);
            
                %put the outputs of each run in a common output variable
                startIndex        = 1+numSP*((m-1)*descreteVariableParam(4)+(l-1));
                endIndex          = numSP*((m-1)*descreteVariableParam(4)+l);
                startIndexAtotal  = 1+numSP*geoSYSp.incoIntTime./geoSYSp.cohIntTime*((m-1)*descreteVariableParam(4)+(l-1));
                endIndexAtotal    = numSP*geoSYSp.incoIntTime./geoSYSp.cohIntTime*((m-1)*descreteVariableParam(4)+l);

                outputSimulations.max_sigPower_noNoise_LR1(startIndex:endIndex)               = outputSimulations_lm.max_sigPower_noNoise_LR1;
                outputSimulations.max_sigPower_noNoise_RR1(startIndex:endIndex)               = outputSimulations_lm.max_sigPower_noNoise_RR1;
                outputSimulations.max_sigPower_Noise_LR1(startIndex:endIndex)                 = outputSimulations_lm.max_sigPower_Noise_LR1;
                outputSimulations.max_sigPower_Noise_RR1(startIndex:endIndex)                 = outputSimulations_lm.max_sigPower_Noise_RR1;

 
                outputSimulations.delayAtSPReferencedToDirect_1_microsec(startIndex:endIndex) = outputSimulations_lm.delayAtSPReferencedToDirect_1_microsec;
                outputSimulations.Atotal_1(startIndexAtotal:endIndexAtotal)                   = outputSimulations_lm.Atotal_1;

                outputSimulations.LR1_UP_Power(startIndex:endIndex)                    = outputSimulations_lm.LR1_UP_Power;
                outputSimulations.RR1_UP_Power(startIndex:endIndex)                    = outputSimulations_lm.RR1_UP_Power;
                outputSimulations.LR1_UP_Power_Fluc(startIndex:endIndex)               = outputSimulations_lm.LR1_UP_Power_Fluc;
                outputSimulations.RR1_UP_Power_Fluc(startIndex:endIndex)               = outputSimulations_lm.RR1_UP_Power_Fluc;
                outputSimulations.PN_star_1(startIndex:endIndex)                       = outputSimulations_lm.PN_star_1;
                outputSimulations.sigPower_noNoise_LR1_lin(:,:,startIndex:endIndex)    = outputSimulations_lm.sigPower_noNoise_LR1_lin;
                outputSimulations.sigPower_noNoise_RR1_lin(:,:,startIndex:endIndex)    = outputSimulations_lm.sigPower_noNoise_RR1_lin;
                outputSimulations.sigPower_Noise_LR1_lin(:,:,startIndex:endIndex)      = outputSimulations_lm.sigPower_Noise_LR1_lin;
                outputSimulations.sigPower_Noise_RR1_lin(:,:,startIndex:endIndex)      = outputSimulations_lm.sigPower_Noise_RR1_lin;

                outputSimulations.NLblackbody_DOWN_LR_1(:,:,startIndex:endIndex)       = outputSimulations_lm.NLblackbody_DOWN_LR_1;
                outputSimulations.NLblackbody_DOWNmean_LR_1(startIndex:endIndex)       = outputSimulations_lm.NLblackbody_DOWNmean_LR_1;
                outputSimulations.NLblackbody_DOWNmax_LR_1(startIndex:endIndex)        = outputSimulations_lm.NLblackbody_DOWNmax_LR_1;
                outputSimulations.NLblackbody_DOWN_RR_1(:,:,startIndex:endIndex)       = outputSimulations_lm.NLblackbody_DOWN_RR_1;
                outputSimulations.NLblackbody_DOWNmean_RR_1(startIndex:endIndex)       = outputSimulations_lm.NLblackbody_DOWNmean_RR_1;
                outputSimulations.NLblackbody_DOWNmax_RR_1(startIndex:endIndex)        = outputSimulations_lm.NLblackbody_DOWNmax_RR_1;
                outputSimulations.NLblackbody_UP_1(:,:,startIndex:endIndex)            = outputSimulations_lm.NLblackbody_UP_1;
                outputSimulations.NLblackbody_UPmean_1(startIndex:endIndex)            = outputSimulations_lm.NLblackbody_UPmean_1;
                outputSimulations.NLblackbody_UPmax_1(startIndex:endIndex)             = outputSimulations_lm.NLblackbody_UPmax_1;

                if geoSYSp.indexTxSignal==1 || geoSYSp.indexTxSignal==3
                    outputSimulations.max_sigPower_noNoise_LR5(startIndex:endIndex)               = outputSimulations_lm.max_sigPower_noNoise_LR5;
                    outputSimulations.max_sigPower_noNoise_RR5(startIndex:endIndex)               = outputSimulations_lm.max_sigPower_noNoise_RR5;
                    outputSimulations.max_sigPower_Noise_LR5(startIndex:endIndex)                 = outputSimulations_lm.max_sigPower_Noise_LR5;
                    outputSimulations.max_sigPower_Noise_RR5(startIndex:endIndex)                 = outputSimulations_lm.max_sigPower_Noise_RR5;

                    outputSimulations.delayAtSPReferencedToDirect_5_microsec(startIndex:endIndex) = outputSimulations_lm.delayAtSPReferencedToDirect_5_microsec;
                    outputSimulations.Atotal_5(startIndexAtotal:endIndexAtotal)                   = outputSimulations_lm.Atotal_5;

                    outputSimulations.LR5_UP_Power(startIndex:endIndex)                    = outputSimulations_lm.LR5_UP_Power;
                    outputSimulations.RR5_UP_Power(startIndex:endIndex)                    = outputSimulations_lm.RR5_UP_Power;
                    outputSimulations.LR5_UP_Power_Fluc(startIndex:endIndex)               = outputSimulations_lm.LR5_UP_Power_Fluc;
                    outputSimulations.RR5_UP_Power_Fluc(startIndex:endIndex)               = outputSimulations_lm.RR5_UP_Power_Fluc;
                    outputSimulations.PN_star_5(startIndex:endIndex)                       = outputSimulations_lm.PN_star_5;
                    outputSimulations.sigPower_noNoise_LR5_lin(:,:,startIndex:endIndex)    = outputSimulations_lm.sigPower_noNoise_LR5_lin;
                    outputSimulations.sigPower_noNoise_RR5_lin(:,:,startIndex:endIndex)    = outputSimulations_lm.sigPower_noNoise_RR5_lin;
                    outputSimulations.sigPower_Noise_LR5_lin(:,:,startIndex:endIndex)      = outputSimulations_lm.sigPower_Noise_LR5_lin;
                    outputSimulations.sigPower_Noise_RR5_lin(:,:,startIndex:endIndex)      = outputSimulations_lm.sigPower_Noise_RR5_lin;

                    outputSimulations.NLblackbody_DOWN_LR_5(:,:,startIndex:endIndex)       = outputSimulations_lm.NLblackbody_DOWN_LR_5;
                    outputSimulations.NLblackbody_DOWNmean_LR_5(startIndex:endIndex)       = outputSimulations_lm.NLblackbody_DOWNmean_LR_5;
                    outputSimulations.NLblackbody_DOWNmax_LR_5(startIndex:endIndex)        = outputSimulations_lm.NLblackbody_DOWNmax_LR_5;
                    outputSimulations.NLblackbody_DOWN_RR_5(:,:,startIndex:endIndex)       = outputSimulations_lm.NLblackbody_DOWN_RR_5;
                    outputSimulations.NLblackbody_DOWNmean_RR_5(startIndex:endIndex)       = outputSimulations_lm.NLblackbody_DOWNmean_RR_5;
                    outputSimulations.NLblackbody_DOWNmax_RR_5(startIndex:endIndex)        = outputSimulations_lm.NLblackbody_DOWNmax_RR_5;
                    outputSimulations.NLblackbody_UP_5(:,:,startIndex:endIndex)            = outputSimulations_lm.NLblackbody_UP_5;
                    outputSimulations.NLblackbody_UPmean_5(startIndex:endIndex)            = outputSimulations_lm.NLblackbody_UPmean_5;
                    outputSimulations.NLblackbody_UPmax_5(startIndex:endIndex)             = outputSimulations_lm.NLblackbody_UPmax_5;
                end


                %increase the Track ID to recognize each run in the L1A,
                %L1B and L2 products
                TrackID_varPar(startIndex:endIndex)       = TrackID+1000000*m+100000*l;
                TrackID_varPar_cohChannel(startIndexAtotal:endIndexAtotal) = TrackID_cohChannel+1000000*m+100000*l; %valid only in Trackwise mode
                %repeat some of the input variables for the number of
                %runs (note that this step is necessary to provide
                %correct inputs to the routine writeSixHours
                GNSS_Week_varPar(startIndex:endIndex)     = GNSS_Week;
                GNSS_SoW_varPar(startIndex:endIndex)      = GNSS_SoW;
                Rx_varPar(startIndex:endIndex,:)          = Rx;
                VRx_varPar(startIndex:endIndex,:)         = VRx;
                Tx_varPar(startIndex:endIndex,:)          = Tx;
                VTx_varPar(startIndex:endIndex,:)         = VTx;
                SPlatlonAlt_HSAVERS=[infoAtSPforRefFile_lm.HSAVERS_SPlat' infoAtSPforRefFile_lm.HSAVERS_SPlon' infoAtSPforRefFile_lm.HSAVERS_SPalt_5km'];
                SPlatlonAlt_HSAVERS_varPar(startIndex:endIndex,:) = SPlatlonAlt_HSAVERS;
                PRN_varPar(startIndex:endIndex)           = PRN;
                SVN_varPar(startIndex:endIndex)           = SVN;

                %create structure to take track of all inputs and outputs of
                %each run
                inOutReferenceFileTemp.TrackID                       = TrackID+1000000*m+100000*l;
                inOutReferenceFileTemp.orbitParameters.GNSS_Week     = GNSS_Week;
                inOutReferenceFileTemp.orbitParameters.GNSS_SoW      = GNSS_SoW;
                inOutReferenceFileTemp.orbitParameters.PRN           = PRN;
                inOutReferenceFileTemp.orbitParameters.SVN           = SVN;
                inOutReferenceFileTemp.orbitParameters.Rx            = Rx;
                inOutReferenceFileTemp.orbitParameters.VRx           = VRx;
                inOutReferenceFileTemp.orbitParameters.Tx            = Tx;
                inOutReferenceFileTemp.orbitParameters.VTx           = VTx;
                
                inOutReferenceFileTemp.inputFilenames.filenameOrbits            = orbitdataHG1;
                inOutReferenceFileTemp.inputFilenames.filenameCoverMap          = filenameCoverMap;
                inOutReferenceFileTemp.inputFilenames.filenameTxGainPattern     = TxGainPattern;
                if satellite == "HydroGNSS-1"  % added by ansha
                    inOutReferenceFileTemp.inputFilenames.filenameRxDOWNGainPattern = nadirantennadataHG1;
                    inOutReferenceFileTemp.inputFilenames.filenameRxUPGainPattern   = zenithantennadataHG1;
                else
                    inOutReferenceFileTemp.inputFilenames.filenameRxDOWNGainPattern = nadirantennadataHG2;
                    inOutReferenceFileTemp.inputFilenames.filenameRxUPGainPattern   = zenithantennadataHG2;
                end
                inOutReferenceFileTemp.outputDirectory                          = dest_dir;
                inOutReferenceFileTemp.outputFilename                           = outputFilename_lm;
                inOutReferenceFileTemp.geoSYSp                                  = geoSYSp;
                inOutReferenceFileTemp.geoSYSp.SPlat_series                     = infoAtSPforRefFile_lm.HSAVERS_SPlat;
                inOutReferenceFileTemp.geoSYSp.SPlon_series                     = infoAtSPforRefFile_lm.HSAVERS_SPlon;
                inOutReferenceFileTemp.geoSYSp.SPalt_5km_series                 = infoAtSPforRefFile_lm.HSAVERS_SPalt_5km;
                inOutReferenceFileTemp.geoSYSp.SPalt_origDEM_series             = infoAtSPforRefFile_lm.HSAVERS_SPalt_origDEM;
                inOutReferenceFileTemp.geoSYSp.TxBeamwidth                      = TRAp.TxBeamwidth;
                inOutReferenceFileTemp.geoSYSp.RxBeamwidth                      = TRAp.RxBeamwidth;
                inOutReferenceFileTemp.geoSYSp.Mode                             = mode1;
                inOutReferenceFileTemp.DEMinfo                                  = DEMinfo;
                inOutReferenceFileTemp.TxSettings.TxSignalPower                 = TRAp.TxPower_dBW;
                inOutReferenceFileTemp.TxSettings.TxSignalPowerError            = TRAp.TxPowerError_dBW;
                inOutReferenceFileTemp.TxSettings.TxAttitudeAngles_nominal      = TRAp.Tx_nomAtt_angles;
                inOutReferenceFileTemp.TxSettings.TxAttitudeAngles_error        = TRAp.TxAttitude_angles_GUI;
                inOutReferenceFileTemp.TxSettings.TxEuler_angles                = TRAp.TxEuler_angles;
                inOutReferenceFileTemp.TxSettings.TxFrequency                   = TRAp.TxFrequency;
                inOutReferenceFileTemp.TxSettings.TxMismatch                    = buttonTxMismatch;

                inOutReferenceFileTemp.TxSettings.Yaw_TX                        = infoAtSPforRefFile_lm.Yaw_TX;

                inOutReferenceFileTemp.RxSettings.RxAttitudeAngles_nominal      = TRAp.Rx_nomAtt_angles;
                inOutReferenceFileTemp.RxSettings.RxAttitudeAngles_error        = TRAp.RxAttitude_angles_GUI;
                inOutReferenceFileTemp.RxSettings.RxEuler_angles                = TRAp.RxEuler_angles;
                inOutReferenceFileTemp.RxSettings.RxDOWNmismatch                = buttonRxDOWNmismatch;
                inOutReferenceFileTemp.RxSettings.RxUPmismatch                  = buttonRxUPmismatch;

                inOutReferenceFileTemp.terrainParameters                        = terrainParameters;
                inOutReferenceFileTemp.vegetationParameters                     = vegetationParameters;

                %look for the cover class where the SP falls and the bio/geo
                %parameters given as input for the SP class
                bioGeoParametersAtSP=biogeoparamAtSP(infoAtSPforRefFile_lm,coverMap_IDs,terrainParameters,vegetationParameters);

                inOutReferenceFileTemp.bioGeoParametersAtSP                     = bioGeoParametersAtSP;

                inOutReferenceFileTemp.Tx1GainAtSP                               = infoAtSPforRefFile_lm.Tx1GainAtSP;
                inOutReferenceFileTemp.Tx1GainAtRx                               = infoAtSPforRefFile_lm.Tx1GainAtRx;
                inOutReferenceFileTemp.RxGainLR1atSP                             = infoAtSPforRefFile_lm.RxGainLR1atSP;
                inOutReferenceFileTemp.RxGainRR1atSP                             = infoAtSPforRefFile_lm.RxGainRR1atSP;
                inOutReferenceFileTemp.squaredDistanceSum                        = infoAtSPforRefFile_lm.squaredDistanceSum;
                inOutReferenceFileTemp.RxZenithGainTowardsTx1                    = UP_pattern_lm(:,1);
                inOutReferenceFileTemp.absoluteDelayAtSP_1_microsec              = outputSimulations_lm.absoluteDelayAtSP_1_microsec;
                inOutReferenceFileTemp.delayAtSPReferencedToDirect_1_microsec    = outputSimulations_lm.delayAtSPReferencedToDirect_1_microsec;

                inOutReferenceFileTemp.max_sigPower_noNoise_LR1                  = outputSimulations_lm.max_sigPower_noNoise_LR1;
                inOutReferenceFileTemp.max_sigPower_noNoise_RR1                  = outputSimulations_lm.max_sigPower_noNoise_RR1;
                inOutReferenceFileTemp.max_sigPower_Noise_LR1                    = outputSimulations_lm.max_sigPower_Noise_LR1;
                inOutReferenceFileTemp.max_sigPower_Noise_RR1                    = outputSimulations_lm.max_sigPower_Noise_RR1;
                inOutReferenceFileTemp.max_reflectivity_noNoise_LR1_dB           = outputSimulations_lm.max_reflectivity_noNoise_LR1_dB;
                inOutReferenceFileTemp.max_reflectivity_noNoise_RR1_dB           = outputSimulations_lm.max_reflectivity_noNoise_RR1_dB;
                inOutReferenceFileTemp.max_reflectivity_NoisyCal_LR1_dB          = outputSimulations_lm.max_reflectivity_Noisy_LR1_dB;
                inOutReferenceFileTemp.max_reflectivity_NoisyCal_RR1_dB          = outputSimulations_lm.max_reflectivity_Noisy_RR1_dB;
                inOutReferenceFileTemp.max_geobReflec_noNoise_LR1_dB             = outputSimulations_lm.max_geobReflec_noNoise_LR1_dB;
                inOutReferenceFileTemp.max_geobReflec_noNoise_RR1_dB             = outputSimulations_lm.max_geobReflec_noNoise_RR1_dB;
                inOutReferenceFileTemp.max_geobPower_Noisy_LR1            = outputSimulations_lm.max_geobPower_Noisy_LR1;
                inOutReferenceFileTemp.max_geobPower_Noisy_RR1            = outputSimulations_lm.max_geobPower_Noisy_RR1;       


%                 inOutReferenceFileTemp.IdealAntenna = outputSimulations_lm.IdealAntenna_all;
%                 inOutReferenceFileTemp.MismatchAntenna = outputSimulations_lm.MismatchAntenna_all;

                if tag.SNR==1
                    inOutReferenceFileTemp.SNR_LR1             = outputSimulations_lm.SNR_LR1;
                    inOutReferenceFileTemp.SNR_RR1             = outputSimulations_lm.SNR_RR1;

%                     inOutReferenceFileTemp.geobSNR_LR1         = outputSimulations_lm.geobSNR_LR1;
%                     inOutReferenceFileTemp.geobSNR_RR1         = outputSimulations_lm.geobSNR_RR1;
                    if geoSYSp.indexTxSignal==1 || geoSYSp.indexTxSignal==3
                        inOutReferenceFileTemp.SNR_LR5             = outputSimulations_lm.SNR_LR5;
                        inOutReferenceFileTemp.SNR_RR5             = outputSimulations_lm.SNR_RR5;

%                         inOutReferenceFileTemp.geobSNR_LR5         = outputSimulations_lm.geobSNR_LR5;
%                         inOutReferenceFileTemp.geobSNR_RR5         = outputSimulations_lm.geobSNR_RR5;
                    end
                end
                if geoSYSp.indexTxSignal==1 || geoSYSp.indexTxSignal==3
                    inOutReferenceFileTemp.Tx5GainAtSP                               = infoAtSPforRefFile_lm.Tx5GainAtSP;
                    inOutReferenceFileTemp.Tx5GainAtRx                               = infoAtSPforRefFile_lm.Tx5GainAtRx;
                    inOutReferenceFileTemp.RxGainLR5atSP                             = infoAtSPforRefFile_lm.RxGainLR5atSP;
                    inOutReferenceFileTemp.RxGainRR5atSP                             = infoAtSPforRefFile_lm.RxGainRR5atSP;
                    inOutReferenceFileTemp.squaredDistanceSum                        = infoAtSPforRefFile_lm.squaredDistanceSum;
                    inOutReferenceFileTemp.RxZenithGainTowardsTx5                    = UP_pattern_lm(:,2);
                    inOutReferenceFileTemp.absoluteDelayAtSP_5_microsec               = outputSimulations_lm.absoluteDelayAtSP_5_microsec;
                    inOutReferenceFileTemp.delayAtSPReferencedToDirect_5_microsec     = outputSimulations_lm.delayAtSPReferencedToDirect_5_microsec;

                    inOutReferenceFileTemp.max_sigPower_noNoise_LR5                  = outputSimulations_lm.max_sigPower_noNoise_LR5;
                    inOutReferenceFileTemp.max_sigPower_noNoise_RR5                  = outputSimulations_lm.max_sigPower_noNoise_RR5;
                    inOutReferenceFileTemp.max_sigPower_Noise_LR5                    = outputSimulations_lm.max_sigPower_Noise_LR5;
                    inOutReferenceFileTemp.max_sigPower_Noise_RR5                    = outputSimulations_lm.max_sigPower_Noise_RR5;
                    inOutReferenceFileTemp.max_reflectivity_noNoise_LR5_dB           = outputSimulations_lm.max_reflectivity_noNoise_LR5_dB;
                    inOutReferenceFileTemp.max_reflectivity_noNoise_RR5_dB           = outputSimulations_lm.max_reflectivity_noNoise_RR5_dB;
                    inOutReferenceFileTemp.max_reflectivity_NoisyCal_LR5_dB          = outputSimulations_lm.max_reflectivity_Noisy_LR5_dB;
                    inOutReferenceFileTemp.max_reflectivity_NoisyCal_RR5_dB          = outputSimulations_lm.max_reflectivity_Noisy_RR5_dB;
                    inOutReferenceFileTemp.max_geobReflec_noNoise_LR5_dB             = outputSimulations_lm.max_geobReflec_noNoise_LR5_dB;
                    inOutReferenceFileTemp.max_geobReflec_noNoise_RR5_dB             = outputSimulations_lm.max_geobReflec_noNoise_RR5_dB;
                    inOutReferenceFileTemp.max_geobPower_Noisy_LR5            = outputSimulations_lm.max_geobPower_Noisy_LR5;
                    inOutReferenceFileTemp.max_geobPower_Noisy_RR5            = outputSimulations_lm.max_geobPower_Noisy_RR5;

                   
                end

                inOutReferenceFile = [inOutReferenceFile inOutReferenceFileTemp];
            end
        end

        % Changed by Mauro since the 'original'0 TrackID value is needed to
        % segment the L1A product in six hour blocks
        TrackID    = uint32(TrackID_varPar);
        TrackID_cohChannel = TrackID_varPar_cohChannel;
        % end of Mauro's change
        GNSS_Week   = GNSS_Week_varPar;
        GNSS_SoW    = GNSS_SoW_varPar;
        Rx          = Rx_varPar;
        VRx         = VRx_varPar;
        Tx          = Tx_varPar;
        VTx         = VTx_varPar;
        SPlatlonAlt_HSAVERS = SPlatlonAlt_HSAVERS_varPar;
        PRN         = PRN_varPar;
        SVN         = SVN_varPar;
        
    end

Gsys_dB = instrumentConfParam.Gsys_dB_LHCP(1);
[UTCDateMidPoint]=WriteSixHours(satellite,TRAp,Gsys_dB,GNSS_Week,GNSS_SoW,Rx,VRx,Tx,VTx,SPlatlonAlt_HSAVERS,PRN,TrackID,TrackID_cohChannel,SVN,outputSimulations,name_data_tag,geoSYSp,outfile_directory);

% Save the Reference File
if OrbitWithReflectivity
    save(strjoin([dest_dir RUNID '_inOutReferenceFile.mat'],''), 'inOutReferenceFile', 'ReflectivitydB_selected');
else 
    save(strjoin([dest_dir RUNID '_inOutReferenceFile.mat'],''), 'inOutReferenceFile');
    save(strjoin([dest_dir  'DataRelease\' RUNID  '_inOutReferenceFile.mat'],''), 'inOutReferenceFile');
end

Path=fullfile('..\conf\');
abspath=fopen([Path 'AbsoluteFilePath.txt'],'wt');
fprintf(abspath,'%s',strjoin(dest_dir,''));%%added by aneesha
fclose(abspath);

disp('.......................................')
disp(' ### SIMULATION COMPLETED ### ')
disp('........................................')    
new=append('RunOutPutPath|',strjoin(dest_dir,''));
disp(new)

else
    disp('.......................................')
    disp(' ### SIMULATION STOPPED BECAUSE OF AN ERROR OR BY THE USER ### ')
    disp('........................................')
    error('Fatal error')
    %return
end

%% starting counter
toc
 end
