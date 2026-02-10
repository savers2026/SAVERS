%% ************************************************************************
% MODULE NAME:      UI_HSAVERS_FirstGUI.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      Aneesha Musunuri
% @copyright:   Tor Vergata University of Rome and Sapienza University of Rome
% Original version:      Sep 2021 - Sep 2022 by Aneesha Musunuri
% Main updates:          Oct 2022 - Dec 2022 by L. Dente; Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel)
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% This routine opens the main GUI of HSAVERS, restore and save the inputs
% in the UserInput.mat file and in the log file
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function closeIniGUI=UI_HSAVERS_FirstGUI(tag)

    global mode0;
    global tag_1;
    global OperationMode;  
    global SensorNoise;
    global RxAnswers;
    global TxAnswers;
    global SimulationSetting;
    global Sampling;
    global Sampling_all;
    global Sampling_all_Text;
    global OutputPlotText;
    global degreeofcoherence;
    global Txswitch;
    global Rxswitch_Zenith;
    global Rxswitch_Nadir;
    global OutputPlot;
    global signal;
    global mode1;
    
    % added by ansha
    global satellite;
    %added by Laura
    global presentLogInputsFilename
    global bioGeoInputs
    global bioGeoInputsVariable
    %
    closeIniGUI    = 0;
    OutputPlot     = [0,0,0,0,0,0];
    OutputPlotText = {0,0,0,0,0,0};
    OperationMode  = {0,0,0};
    mode0          = {0,0,0};
    signal         = "GPS L1";
    mode1           = "Single Point";
    Sampling       = {0,0,0};
    Sampling_all   = {0,0};
    TxAnswers      = [0,0,0,0,0,0,0];
    RxAnswers      = [0,0,0,0,0,0];
    SensorNoise    = [0,0,0,0,0,0,0];
    %SimulationSetting = [0,0,0,0,0,0,0,0,0,0];
    SimulationSetting = [0,0,0,0,0,0,0,0,0]; % 10 values
    Sampling_all_Text = {0,0};
    degreeofcoherence = (0);
    Txswitch          = 'Mismatch';
    Rxswitch_Zenith   = 'ZenithMismatch';
    Rxswitch_Nadir    = 'NadirMismatch';
    satellite         = {0};
    % added by Laura
    bioGeoInputs.water  = {'Normal','0.1', '20.0'};
    bioGeoInputs.BLForest   = {'Normal','1.0','5','30','1.18', '20.0', '5.','300.','25.','10.', '150.'};
    bioGeoInputs.NLForest   = {'Normal','1.0','5','15.0','1.18','20.0','5.','600.','15.','10.', '150.'};
    bioGeoInputs.LowVeg     = {'Normal','1.0','5','15','1.18','20.0','2.','0.75', '10.', '4.'};
    bioGeoInputs.bare       = {'Normal','1.0','5','15','1.18','20.'};
    bioGeoInputs.floodedForest  = {'Flooded','0.1','5','80','1.18', '20.0', '5.','300.','25.','10.', '150.'};
    bioGeoInputs.floodedLowVeg  = {'Flooded','0.1','5','80','1.18','20.0','2.','0.75', '10.', '4.'};
    bioGeoInputsVariable.sm.inputsVar   = {'1.0','41.0','10','3'};
    bioGeoInputsVariable.sm.BLForest    = {'0.3','2.0','25.','275.','5','1.18','5.','300.','25.','10.'};
    bioGeoInputsVariable.sm.NLForest    = {'0.3','2.0', '25.','275','5','1.18','5.','600.','15.','10.'};
    bioGeoInputsVariable.sm.LowVeg      = {'0.3','3.0','0.','6.','5','1.18','0.75', '10.', '4.'};
    bioGeoInputsVariable.sm.bare        = {'0.3','3.0','5','1.18'};
    bioGeoInputsVariable.sm.floodedForest  = {'0.3','2.0','25.','275.','5','1.18','5.','300.','25.','10.'};
    bioGeoInputsVariable.sm.floodedLowVeg  = {'0.3','3.0','0.','6.','5','1.18','0.75', '10.', '4.'};
    bioGeoInputsVariable.bio.inputsVar  = {'25.0','275.0','50','3'};
    bioGeoInputsVariable.bio.BLForest   = {'0.3','2.0','10.0','40.','5','1.18','5.','300.','25.','10.'};
    bioGeoInputsVariable.bio.NLForest   = {'0.3','2.0','10.0','40.','5','1.18','5.','600.','15.','10.'};
    bioGeoInputsVariable.bio.LowVeg     = {'0.3','3.0','1.0','40.','0.','6.','5','1.18','0.75', '10.', '4.'};
    bioGeoInputsVariable.bio.bare       = {'0.3','3.0','1.0','40.','5','1.18'};
    bioGeoInputsVariable.bio.floodedForest  = {'0.3','2.0','10.0','40.','5','1.18','5.','300.','25.','10.'};
    bioGeoInputsVariable.bio.floodedLowVeg  = {'0.3','3.0','1.0','40.','0.','6.','5','1.18','0.75', '10.', '4.'};
    bioGeoInputsVariable.ft.inputsVar   = {'3'};
    bioGeoInputsVariable.ft.BLForest    = {'0.3','2.0','10.0','40.','25.','275.','5','1.18','5.','300.','25.','10.'};
    bioGeoInputsVariable.ft.NLForest    = {'0.3','2.0','10.0','40.','25.','275.','5','1.18','5.','600.','15.','10.'};
    bioGeoInputsVariable.ft.LowVeg      = {'0.3','3.0','1.0','40.','0.','6.','5','1.18','0.75', '10.', '4.'};
    bioGeoInputsVariable.ft.bare        = {'0.3','3.0','1.0','40.','5','1.18'};
    bioGeoInputsVariable.ft.floodedForest  = {'0.3','2.0','10.0','40.','25.','275.','5','1.18','5.','600.','15.','10.'};
    bioGeoInputsVariable.ft.floodedLowVeg  = {'0.3','3.0','1.0','40.','0.','6.','5','1.18','0.75', '10.', '4.'};
    bioGeoInputsVariable.wet.inputsVar  = {'3'};
    bioGeoInputsVariable.wet.BLForest   = {'0.3','2.0','10.0','40.','25.','275','5','1.18','5.','300.','25.','10.'};
    bioGeoInputsVariable.wet.NLForest   = {'0.3','2.0','10.0','40.','25.','275.','5','1.18','5.','600.','15.','10.'};
    bioGeoInputsVariable.wet.LowVeg     = {'0.3','3.0','1.0','40.','0.','6.','5','1.18','0.75','10.', '4.'};
    bioGeoInputsVariable.wet.bare       = {'0.3','3.0','1.0','40.','5','1.18'};
    bioGeoInputsVariable.wet.floodedForest  = {'0.3','2.0','1.0','40.','25.','275.','5','1.18','5.','300.','25.','10.'};
    bioGeoInputsVariable.wet.floodedLowVeg  = {'0.3','3.0','1.0','40.','0.','6.','5','1.18','0.75','10.', '4.'};
    %

    fig        = uifigure('Name','HSAVERS by Sapienza University of Rome and Tor Vergata University of Rome'); %the main GUI window
    fig.Units  = 'pixels';
    fig.Resize = 'off';
    %fig.Position = [150 303 1200 500];
    fig.Position  = [150 303 1200 550];
    movegui(fig,'center');
    
    bg1     = uibuttongroup(fig,'Title', 'Mode','Position',[10 420 210 120],'SelectionChangedFcn',@(bg1,event) modeselection(bg1)); %Mode group
    tb1_bg1 = uitogglebutton(bg1,'Text','Single Point','Position',[10 64 180 22],'Tag','1');
    tb2_bg1 = uitogglebutton(bg1,'Text','Trackwise','Position',[10 44 180 22],'Tag','2');
    tb3_bg1 = uitogglebutton(bg1,'Text','Monte Carlo','Position',[10 24 180 22],'Tag','3');
    tb4_bg1 = uitogglebutton(bg1,'Text','Global','Position',[10 04 180 22],'Tag','4');
    
    bg2     = uibuttongroup(fig,'Title', 'Signal','Position',[230 420 210 120],'SelectionChangedFcn',@(bg2,event) signalselection(bg2)); %signal group
    tb1_bg2 = uitogglebutton(bg2,'Text','GPS L1','Position',[10 64 180 22],'Tag','1');
    tb2_bg2 = uitogglebutton(bg2,'Text','GPS L1&L5','Position',[10 44 180 22],'Tag','2');
    tb3_bg2 = uitogglebutton(bg2,'Text','Galileo E1c','Position',[10 24 180 22],'Tag','3');
    tb4_bg2 = uitogglebutton(bg2,'Text','Galileo E1c&E5a','Position',[10 04 180 22],'Tag','4');
    
    cd1 = uicheckbox(fig, 'Text','Coarse mode','Position',[465 490 180 15],'ValueChangedFcn',@(cd1,event) changedmodeoperation(cd1,1));     %coarse mode
    cd2 = uicheckbox(fig, 'Text','Variable bio/geo params','Position',[465 470 180 15],'ValueChangedFcn',@(cd2,event) changedmodeoperation(cd2,2));      %Verbose mode
    cd3 = uicheckbox(fig, 'Text','Flat Topography','Position',[465 450 180 15],'ValueChangedFcn',@(cd3,event) changedmodeoperation(cd3,3));    %Flat Topography
    
    
    m       = uimenu(fig,'Text','&Import');   %Import
    mitem1  = uimenu(m,'Text','&Platform Files');
    mchild1 = uimenu(mitem1,'Text','&Orbit Data HG1','MenuSelectedFcn',@(mchild1,event) MenuSelected_orbitdata(mchild1));
    mchild2 = uimenu(mitem1,'Text','&Orbit Data HG2','MenuSelectedFcn',@(mchild2,event) MenuSelected_metadata(mchild2));
    mitem2  = uimenu(m,'Text','&Antenna Patterns');
    mchild3 = uimenu(mitem2,'Text','&Nadir Antenna HG1','MenuSelectedFcn',@(mchild3,event) MenuSelected_Nadirhg1(mchild3));
    mchild4 = uimenu(mitem2,'Text','&Zenith Antenna HG1','MenuSelectedFcn',@(mchild4,event) MenuSelected_Zenithhg1(mchild4));
    mchild5 = uimenu(mitem2,'Text','&Nadir Antenna HG2','MenuSelectedFcn',@(mchild5,event) MenuSelected_Nadirhg2(mchild5));
    mchild6 = uimenu(mitem2,'Text','&Zenith Antenna HG2','MenuSelectedFcn',@(mchild6,event) MenuSelected_Zenithhg2(mchild6));
    
    % mitem1_m2.MenuSelectedFcn = @MenuSelected;
    % m3 = uimenu(fig,'Text','&Help');
    % mitem1_m3 = uimenu(m3,'Text','&Documentation');
    
    pn = uipanel(fig,'Title','Site Selection','Position',[640 420 280 120]);
    dlabel1_pn = uilabel (pn,'Text','Soil Moisture','HorizontalAlignment','right','Position',[5 80 80 15]);
    dd_pn = uidropdown(pn,'Position',[90 80 150 20],'Items',{'Not Selected','Tibesti Mountains','San Luis Valley','Site 3'});
    dd_pn.ValueChangedFcn = @(dd_pn,event) dropdownmenuchanged_soilmoisture(dd_pn);
    dlabel2_pn = uilabel (pn,'Text','Forest','HorizontalAlignment','right','Position',[5 55 80 15]);
    dd1_pn = uidropdown(pn,'Position',[90 55 150 20],'Items',{'Not Selected','Gabon','Finland','Site 3'});
    dd1_pn.ValueChangedFcn = @(dd1_pn,event) dropdownmenuchanged_forest(dd1_pn);
    dlabel3_pn = uilabel (pn,'Text','Freeze/Thaw','HorizontalAlignment','right','Position',[5 30 80 15]);
    dd2_pn = uidropdown(pn,'Position',[90 30 150 20],'Items',{'Not Selected','Finland','Chersky','Site 3'});
    dd2_pn.ValueChangedFcn = @(dd2_pn,event) dropdownmenuchanged_freeze_thaw(dd2_pn);
    dlabel4_pn = uilabel (pn,'Text','Wetlands','HorizontalAlignment','right','Position',[5 05 80 15]);
    dd3_pn = uidropdown(pn,'Position',[90 05 150 20],'Items',{'Not Selected','Yucatan Lake','Everglades','Site 3'});
    dd3_pn.ValueChangedFcn = @(dd3_pn,event) dropdownmenuchanged_wetlands(dd3_pn);
    
    pn1 = uipanel(fig,'Title','Simulation Settings','Position',[10 5 210 390]);
    dlabel1_pn1 = uilabel (pn1,'Text','Simulation Resolution  [m] ','Position',[10 343 200 15]);
    edt8_pn1 = uieditfield(pn1,'numeric','Tag','Simulation Resolution ','Position',[10 326 150 19]);
    edt8_pn1.ValueChangedFcn=@(edt8_pn1,event) textChangedss(edt8_pn1,1);
%     dlabel2_pn1 = uilabel (pn1,'Text','Simulation Resolution Y [m]','Position',[10 309 200 15]); [10 309 200 15]
%     edt1_pn1 = uieditfield(pn1,'numeric','Tag','Simulation Resolution Y','Position',[10 292 150 19]);[10 292 150 19]
%     edt1_pn1.ValueChangedFcn=@(edt1_pn1,event) textChangedss(edt1_pn1,2);
    dlabel3_pn1 = uilabel (pn1,'Text','Doppler Extent [Hz]','Position',[10 309 200 15]);
    edt2_pn1 = uieditfield(pn1,'numeric','Tag','Doppler Extent [Hz]','Position',[10 292 150 19]);
    edt2_pn1.ValueChangedFcn=@(edt2_pn1,event) textChangedss(edt2_pn1,2);%added by ansha - 12052033
    dlabel4_pn1 = uilabel (pn1,'Text','Delay Extent L1\E1 [µs]','Position',[10 275 200 15]);
    edt3_pn1 = uieditfield(pn1,'numeric','Tag','Delay Extent L1\E1 [µs]','Position',[10 258 150 19]);
    edt3_pn1.ValueChangedFcn=@(edt3_pn1,event) textChangedss(edt3_pn1,3);%added by ansha - 12052033
    dlabel10_pn1 = uilabel (pn1,'Text','Delay Extent L5\E5 [µs]','Position',[10 241 200 15]);
    edt10_pn1 = uieditfield(pn1,'numeric','Tag','Delay Extent L5\E5 [µs]','Position',[10 224 150 19]);
    edt10_pn1.ValueChangedFcn=@(edt10_pn1,event) textChangedss(edt10_pn1,4);
    dlabel5_pn1 = uilabel (pn1,'Text','Delay Bins [#]','Position',[10 207 200 15]);
    edt4_pn1 = uieditfield(pn1,'numeric','Tag','Delay Bins [#]','Position',[10 190 150 19]);
    edt4_pn1.ValueChangedFcn=@(edt4_pn1,event) textChangedss(edt4_pn1,5);
    dlabel6_pn1 = uilabel (pn1,'Text','Doppler Bins [#]','Position',[10 173 200 15]);
    edt5_pn1 = uieditfield(pn1,'numeric','Tag','Doppler Bins [#]','Position',[10 156 150 19]);
    edt5_pn1.ValueChangedFcn=@(edt5_pn1,event) textChangedss(edt5_pn1,6);
    dlabel7_pn1 = uilabel (pn1,'Text','Tx near-specular beamwidth [deg]','Position',[10 139 200 15]);
    edt6_pn1 = uieditfield(pn1,'numeric','Tag','Tx near-specular beamwidth','Position',[10 122 150 19]);
    edt6_pn1.ValueChangedFcn=@(edt6_pn1,event) textChangedss(edt6_pn1,7);
%     dlabel8_pn1 = uilabel (pn1,'Text','Rx near-specular beamwidth [deg] (0:Auto calculation)','Position',[10 105 200 15]);
    dlabel8_pn1 = uilabel (pn1,'Text','Beta_c[deg](0:Auto)[amir test]','Position',[10 105 200 15]);
    edt7_pn1 = uieditfield(pn1,'numeric','Tag','Rx near-specular beamwidth','Position',[10 88 150 19] );
    edt7_pn1.ValueChangedFcn=@(edt7_pn1,event) textChangedss(edt7_pn1,8);
    dlabel9_pn1 = uilabel (pn1,'Text','SNR Booster Δpower [dB]','Position',[10 71 200 15]);
    edt9_pn1 = uieditfield(pn1,'numeric','Tag','SNR Booster','Position',[10 54 150 19]);
    edt9_pn1.ValueChangedFcn=@(edt9_pn1,event) textChangedss(edt9_pn1,9);
%     dlabel10_pn1 = uilabel (pn1,'Text','Atmo Error','Position',[10 37 200 15]);
%     edt10_pn1 = uieditfield(pn1,'numeric','Tag','Atmo Error','Position',[10 20 150 19]);
%     edt10_pn1.ValueChangedFcn=@(edt10_pn1,event) textChangedss(edt10_pn1,10);
    
    pn2 = uipanel(fig,'Title','Tx Settings','Position',[230 5 210 390]);
    dlabel1_pn2 = uilabel (pn2,'Text','Tx Power error [dBW] ','Position',[10 343 200 15]);
    edt8_pn2 = uieditfield(pn2,'numeric','Tag','Tx Power error','Position',[10 326 150 19]);
    edt8_pn2.ValueChangedFcn=@(edt8_pn2,event) textChangedtx(edt8_pn2,1);
    dlabel2_pn2 = uilabel (pn2,'Text','Tx Roll error [deg]','Position',[10 309 200 15]);
    edt1_pn2 = uieditfield(pn2,'numeric','Tag','Tx Roll error ','Position',[10 292 150 19]);
    edt1_pn2.ValueChangedFcn=@(edt1_pn2,event) textChangedtx(edt1_pn2,2);
    dlabel3_pn2 = uilabel (pn2,'Text','Tx Pitch error [deg]','Position',[10 275 200 15]);
    edt2_pn2 = uieditfield(pn2,'numeric','Tag','Tx Pitch error','Position',[10 258 150 19]);
    edt2_pn2.ValueChangedFcn=@(edt2_pn2,event) textChangedtx(edt2_pn2,3);
    dlabel4_pn2 = uilabel (pn2,'Text','Tx Yaw error [deg]','Position',[10 241 200 15]);
    edt3_pn2 = uieditfield(pn2,'numeric','Tag','Tx Yaw error','Position',[10 224 150 19]);
    edt3_pn2.ValueChangedFcn=@(edt3_pn2,event) textChangedtx(edt3_pn2,4);
    dlabel5_pn2 = uilabel (pn2,'Text','Tx EulerAngle phi [deg]','Position',[10 207 200 15]);
    edt4_pn2 = uieditfield(pn2,'numeric','Tag','Tx EulerAngle phi','Position',[10 190 150 19]);
    edt4_pn2.ValueChangedFcn=@(edt4_pn2,event) textChangedtx(edt4_pn2,5);
    dlabel6_pn2 = uilabel (pn2,'Text','Tx EulerAngle theta [deg]','Position',[10 173 200 15]);
    edt5_pn2 = uieditfield(pn2,'numeric','Tag','Tx EulerAngle theta','Position',[10 156 150 19]);
    edt5_pn2.ValueChangedFcn=@(edt5_pn2,event) textChangedtx(edt5_pn2,6);
    dlabel7_pn2 = uilabel (pn2,'Text','Tx EulerAngle psi [deg]','Position',[10 139 200 15]);
    edt6_pn2 = uieditfield(pn2,'numeric','Tag','Tx EulerAngle psi','Position',[10 122 150 19]);
    edt6_pn2.ValueChangedFcn=@(edt6_pn2,event) textChangedtx(edt6_pn2,7);
    sw1_pn2 = uiswitch(pn2,'Items',{'Mismatch','Ideal'},'Position',[80 90 150 15],'ValueChangedFcn',@(sw1_pn2,event) Tx_switchvalue(sw1_pn2));
    
    pn3 = uipanel(fig,'Title','Rx Settings','Position',[450 5 210 390]);
    dlabel1_pn3 = uilabel (pn3,'Text','Rx Roll error [deg]','Position',[10 343 200 15]);
    edt8_pn3 = uieditfield(pn3,'numeric','Tag','Rx Roll error','Position',[10 326 150 19]);
    edt8_pn3.ValueChangedFcn=@(edt8_pn3,event) textChangedrx(edt8_pn3,1);
    dlabel2_pn3 = uilabel (pn3,'Text','Rx Pitch error [deg] ','Position',[10 309 200 15]);
    edt1_pn3 = uieditfield(pn3,'numeric','Tag','Rx Pitch error','Position',[10 292 150 19]);
    edt1_pn3.ValueChangedFcn=@(edt1_pn3,event) textChangedrx(edt1_pn3,2);
    dlabel3_pn3 = uilabel (pn3,'Text','Rx Yaw error [deg]','Position',[10 275 200 15]);
    edt2_pn3 = uieditfield(pn3,'numeric','Tag','Rx Yaw','Position',[10 258 150 19]);
    edt2_pn3.ValueChangedFcn=@(edt2_pn3,event) textChangedrx(edt2_pn3,3);
    dlabel4_pn3 = uilabel (pn3,'Text','Rx nominal Roll [deg]','Position',[10 241 200 15]);
    edt3_pn3 = uieditfield(pn3,'numeric','Tag','Rx nominal Roll','Position',[10 224 150 19]);
    edt3_pn3.ValueChangedFcn=@(edt3_pn3,event) textChangedrx(edt3_pn3,4);
    dlabel5_pn3 = uilabel (pn3,'Text','Rx nominal Pitch [deg] ','Position',[10 207 200 15]);
    edt4_pn3 = uieditfield(pn3,'numeric','Tag','Rx nominal Pitch','Position',[10 190 150 19]);
    edt4_pn3.ValueChangedFcn=@(edt4_pn3,event) textChangedrx(edt4_pn3,5);
    dlabel6_pn3 = uilabel (pn3,'Text','Rx nominal Yaw [deg]','Position',[10 173 200 15]);
    edt5_pn3 = uieditfield(pn3,'numeric','Tag','Rx nominal Yaw','Position',[10 156 150 19]);
    edt5_pn3.ValueChangedFcn=@(edt5_pn3,event) textChangedrx(edt5_pn3,6);
    sw2_pn3 = uiswitch(pn3,'Items',{'NadirMismatch','Ideal'},'Position',[110 126 50 15],'ValueChangedFcn',@(sw2_pn3,event) Rx_switchvalue_1(sw2_pn3));
    sw3_pn3 = uiswitch(pn3,'Items',{'ZenithMismatch','Ideal'},'Position',[110 103 130 15],'ValueChangedFcn',@(sw3_pn3,event) Rx_switchvalue_2(sw3_pn3));
    
    pn4 = uipanel(fig,'Title','Output Plots','Position',[900 50 200 160]);
    cb1_pn4 = uicheckbox(fig, 'Text','Map of selected SPs','Position',[920 170 150 20],'ValueChangedFcn',@(cb1_pn4,event) outputplotchange(cb1_pn4,1));
    cb2_pn4 = uicheckbox(fig, 'Text','DEM','Position',[920 150 150 20],'ValueChangedFcn',@(cb2_pn4,event) outputplotchange(cb2_pn4,2));
    cb3_pn4 = uicheckbox(fig, 'Text','Cover map','Position',[920 130 150 20],'ValueChangedFcn',@(cb3_pn4,event) outputplotchange(cb3_pn4,3));
    cb4_pn4 = uicheckbox(fig, 'Text','Angles & Isolines','Position',[920 110 150 20],'ValueChangedFcn',@(cb4_pn4,event) outputplotchange(cb4_pn4,4));
    cb5_pn4 = uicheckbox(fig, 'Text','Sigma0','Position',[920 90 150 20],'ValueChangedFcn',@(cb5_pn4,event) outputplotchange(cb5_pn4,5));
    cb6_pn4 = uicheckbox(fig, 'Text','DDM','Position',[920 70 150 20],'ValueChangedFcn',@(cb6_pn4,event) outputplotchange(cb6_pn4,6));
    
    dlabel_ri = uilabel (fig,'Text','Run ID ','Position',[910 225 200 25]);
    edt8_ri = uieditfield(fig,'Text','Tag','RUN ID','Position',[950 226 150 20]);
    edt8_ri.ValueChangedFcn=@(edt8_ri,event) textchangedrunid(edt8_ri);
    
    pn5 = uipanel(fig,'Title','Sensor Settings','Position',[670 5 210 390]);
    dlabel1_pn5 = uilabel (pn5,'Text','Coherent Integration Time [msec] ','Position',[10 343 200 15]);
    edt8_pn5 = uieditfield(pn5,'numeric','Tag','Rx Coherent Integration Time','Position',[10 326 150 19]);
    edt8_pn5.ValueChangedFcn=@(edt8_pn5,event) textChangedsn(edt8_pn5,1);
    dlabel2_pn5 = uilabel (pn5,'Text','Incoherent Integration Time [msec] ','Position',[10 309 200 15]);
    edt1_pn5 = uieditfield(pn5,'numeric','Tag','Rx Incoherent Integration Time','Position',[10 292 150 19]);
    edt1_pn5.ValueChangedFcn=@(edt1_pn5,event) textChangedsn(edt1_pn5,2);
    dlabel3_pn5 = uilabel (pn5,'Text','Noise Figure error [dB]','Position',[10 275 200 15]);
    edt2_pn5 = uieditfield(pn5,'numeric','Tag','Noise Figure error [dB]', 'Position',[10 258 150 19]);
    edt2_pn5.ValueChangedFcn=@(edt2_pn5,event) textChangedsn(edt2_pn5,3);
    dlabel4_pn5 = uilabel (pn5,'Text','Nadir Antenna Temp TA [K]','Position',[10 241 200 15]);
    edt3_pn5 = uieditfield(pn5,'numeric','Tag','Nadir Antenna Temp TA','Position',[10 224 150 19]);
    edt3_pn5.ValueChangedFcn=@(edt3_pn5,event) textChangedsn(edt3_pn5,4);
    dlabel3a_pn5 = uilabel(pn5,'Text','Zenith Antenna Temp Tsky [K]','Position',[10 207 200 15]);
    edt3a_pn5=uieditfield(pn5,'numeric','Tag','Zenith Antenna Temp Tski [K]','Position',[10 190 150 19]);
    edt3a_pn5.ValueChangedFcn=@(edt3a_pn5,event) textChangedsn(edt3a_pn5,5);
    dlabel5_pn5 = uilabel (pn5,'Text','Receiver Temperature error [K]','Position',[10 173 200 15]);
    edt4_pn5 = uieditfield(pn5,'numeric','Tag','Receiver Temperature error [K]','Position',[10 156 150 19]);
    edt4_pn5.ValueChangedFcn=@(edt4_pn5,event) textChangedsn(edt4_pn5,6);
    dlabel6_pn5 = uilabel (pn5,'Text','Gain error [dB]','Position',[10 139 200 15]);
    edt5_pn5 = uieditfield(pn5,'numeric','Tag','Gain error [dB]','Position',[10 122 150 19]);
    edt5_pn5.ValueChangedFcn=@(edt5_pn5,event) textChangedsn(edt5_pn5,7);
     dlabelsl_pn5 = uilabel (pn5,'Text','Degree of Coherence','Position',[10 92 200 15]);
    sl1_pn5 = uislider(pn5,'Limits',[0 1],'MajorTicks',[0 0.5 1],'Value',0,'Position',[40 80 60 3],'ValueChangedFcn',@(sl1_pn5,event) slidervaluechange(sl1_pn5));


    m2 = uimenu(fig,'Text','&Default');%Default
    mitem1_m2 = uimenu(m2,'Text','&Previous GUI Input File','Accelerator','P');
    mitem1_m2.MenuSelectedFcn = @(mitem1_m2,event) loadinputfile(mitem1_m2);
    mitem2_m2 = uimenu(m2,'Text','&Auxiliary Folder','Accelerator','A');
    mitem2_m2.MenuSelectedFcn= @(mitem2_m2,event) loadinputfilespath(mitem2_m2);
    mitem3_m2 = uimenu(m2,'Text','&Output Folder ','Accelerator','O');
    mitem3_m2.MenuSelectedFcn = @(mitem3_m2,event) loadoutfilespath(mitem3_m2);


    %added by ansha - 12052033
    bg3     = uibuttongroup(fig,'Title', 'Satellite Selection','Position',[950 450 150 90],'SelectionChangedFcn',@(bg3,event) satelliteselection(bg3));
    tb1_bg3 = uitogglebutton(bg3,'Text','HydroGNSS-1','Position',[30 42 100 22],'Tag','1');
    tb2_bg3 = uitogglebutton(bg3,'Text','HydroGNSS-2','Position',[30 22 100 22],'Tag','2');

    continuebutton = uibutton(fig,'Push','Text','Continue','Position',[980 415 100 20],'ButtonPushedFcn',@(continuebutton,event) saveduserinput(fig));
    load_val_last = uibutton(fig,'Push','Text','Previous Input','Position',[980 390 100 20],'ButtonPushedFcn',@(load_val_last,event) loadinputfile_previous(load_val_last));
    
     pn7 = uipanel(fig,'Title','Global Orbitdata Sub Sampling ','Position',[900 265 200 105]);
     dlabel1_pn7 = uilabel (pn7,'Text','First_Day','Position',[10 60 190 15]);
     edt1_pn7=uieditfield(pn7,'numeric','Tag','First_Day','Position',[20 38 40 20]);
     edt1_pn7.ValueChangedFcn=@(edt1_pn7,event) textChangedoss(edt1_pn7,1);
     dlabel2_pn7 = uilabel(pn7,'Text','Last_Day','Position',[70 60 190 15]);
     edt2_pn7=uieditfield(pn7,'numeric','Tag','Last_Day','Position',[80 38 40 20]);
     edt2_pn7.ValueChangedFcn=@(edt2_pn7,event) textChangedoss(edt2_pn7,2);
     dlabel3_pn7=uilabel(pn7,'Text','Num_SP','Position',[130 60 190 15]);
     edt3_pn7=uieditfield(pn7,'numeric','Tag','Num_SP','Position',[140 38 40 20 ]);
     edt3_pn7.ValueChangedFcn=@(edt3_pn7,event) textChangedoss(edt3_pn7,3);
     cb1_pn7=uicheckbox(fig, 'Text','All Days','Position',[925 275 100 20],'ValueChangedFcn',@(cb1_pn7,event) change_day_sp(cb1_pn7,edt1_pn7,edt2_pn7,edt3_pn7,1));
     cb2_pn7=uicheckbox(fig, 'Text','All SPs','Position',[998 275 100 20],'ValueChangedFcn',@(cb2_pn7,event) change_day_sp(cb2_pn7,edt1_pn7,edt2_pn7,edt3_pn7,2));
    
    %tag_1 =[edt8_pn1,edt1_pn1,edt2_pn1,edt3_pn1,edt4_pn1,edt5_pn1,edt6_pn1,edt7_pn1,edt8_pn2,edt1_pn2,edt2_pn2,edt3_pn2,edt4_pn2,edt5_pn2,edt6_pn2,edt8_pn3,edt1_pn3,edt2_pn3,edt3_pn3,edt4_pn3,edt5_pn3,edt8_pn5,edt1_pn5,edt2_pn5,edt3_pn5,edt4_pn5,bg1,bg2,sw1_pn2,sw2_pn3,sw3_pn3,sl1_pn5,tb1_bg1,tb2_bg1,tb3_bg1,tb4_bg1,tb1_bg2,tb2_bg2,tb3_bg2,tb4_bg2,cd1,cd2,cd3,edt8_ri,dd_pn,dd1_pn,dd2_pn,dd3_pn,mchild3,mchild4,mchild1,mchild2,mitem2_m2,mitem3_m2,cb1_pn4,cb2_pn4,cb3_pn4,cb4_pn4,cb5_pn4,cb6_pn4,edt1_pn7,edt2_pn7,edt3_pn7,cb1_pn7,cb2_pn7,edt3a_pn5,edt5_pn5,edt9_pn1,edt10_pn1];%orbitdata,z_m,data_directory,outfile_directory];
%     tag_1 =[edt8_pn1,edt1_pn1,edt2_pn1,edt3_pn1,edt4_pn1,edt5_pn1,edt6_pn1,edt7_pn1,edt8_pn2,edt1_pn2,edt2_pn2,edt3_pn2,edt4_pn2,edt5_pn2,edt6_pn2,edt8_pn3,edt1_pn3,edt2_pn3,edt3_pn3,edt4_pn3,edt5_pn3,edt8_pn5,edt1_pn5,edt2_pn5,edt3_pn5,edt4_pn5,bg1,bg2,sw1_pn2,sw2_pn3,sw3_pn3,sl1_pn5,tb1_bg1,tb2_bg1,tb3_bg1,tb4_bg1,tb1_bg2,tb2_bg2,tb3_bg2,tb4_bg2,cd1,cd2,cd3,edt8_ri,dd_pn,dd1_pn,dd2_pn,dd3_pn,mchild3,mchild4,mchild1,mchild2,mitem2_m2,mitem3_m2,cb1_pn4,cb2_pn4,cb3_pn4,cb4_pn4,cb5_pn4,cb6_pn4,edt1_pn7,edt2_pn7,edt3_pn7,cb1_pn7,cb2_pn7,edt3a_pn5,edt5_pn5,edt9_pn1,edt10_pn1,bg3,tb1_bg3,tb2_bg3];%orbitdata,z_m,data_directory,outfile_directory];
      tag_1 =[edt8_pn1,edt2_pn1,edt3_pn1,edt4_pn1,edt5_pn1,edt6_pn1,edt7_pn1,edt8_pn2,edt1_pn2,edt2_pn2,edt3_pn2,edt4_pn2,edt5_pn2,edt6_pn2,edt8_pn3,edt1_pn3,edt2_pn3,edt3_pn3,edt4_pn3,edt5_pn3,edt8_pn5,edt1_pn5,edt2_pn5,edt3_pn5,edt4_pn5,bg1,bg2,sw1_pn2,sw2_pn3,sw3_pn3,sl1_pn5,tb1_bg1,tb2_bg1,tb3_bg1,tb4_bg1,tb1_bg2,tb2_bg2,tb3_bg2,tb4_bg2,cd1,cd2,cd3,edt8_ri,dd_pn,dd1_pn,dd2_pn,dd3_pn,mchild3,mchild4,mchild1,mchild2,mitem2_m2,mitem3_m2,cb1_pn4,cb2_pn4,cb3_pn4,cb4_pn4,cb5_pn4,cb6_pn4,edt1_pn7,edt2_pn7,edt3_pn7,cb1_pn7,cb2_pn7,edt3a_pn5,edt5_pn5,edt9_pn1,edt10_pn1,bg3,tb1_bg3,tb2_bg3,mchild5,mchild6];%orbitdata,z_m,data_directory,outfile_directory];

    fig.CloseRequestFcn = @(closeBut,event)my_closereq(closeBut);

    uiwait(fig);
    
    function my_closereq(fig)
        delete(fig)
        closeIniGUI=1;
        return
    end

    function MenuSelected_orbitdata(txt)
        global orbitdataHG1;
        global data_directory;
        pathname_orbitfiles =  data_directory;
        titleDataFileSelection = 'Select the Orbitdata file for HG1';
        [filename,pathname_orbitfiles] = uigetfile({'*.nc';'*.mat'}, titleDataFileSelection,(pathname_orbitfiles));
        orbitdataHG1 = filename;
        name= append('OrbitFile HG1:',orbitdataHG1);
        txt.Text= name;
        end
    function MenuSelected_metadata(txt)
        global orbitdataHG2;
        global data_directory;
        pathname_metafiles = data_directory;
        titleDataFileSelection = 'Select the Orbitdata file for HG2';
        [filename, pathname_metafiles] =  uigetfile({'*.nc';'*.mat'}, titleDataFileSelection,(pathname_metafiles));
        orbitdataHG2 =  filename;
        name= append('OrbitFile HG2:',orbitdataHG2);
        txt.Text= name;
        
    end
    function MenuSelected_Nadirhg1(txt)
        global data_directory;
        global nadirantennadataHG1;
        path_nadir = data_directory;
        titleDataFileSelection = 'Select Nadir Antenna Pattern File for HG1';
        [filename,path_nadir] = uigetfile({'*.xlsx';'*.mat';'*.nc'}, titleDataFileSelection,(path_nadir));
        nadirantennadataHG1 = filename;
        name=append('Nadir HG1:',nadirantennadataHG1);
        txt.Text=name;
         
    end
    function MenuSelected_Nadirhg2(txt)
        global data_directory;
        global nadirantennadataHG2;
        path_nadir = data_directory;
        titleDataFileSelection = 'Select Nadir Antenna Pattern File for HG2';
        [filename,path_nadir] = uigetfile({'*.xlsx';'*.mat';'*.nc'}, titleDataFileSelection,(path_nadir));
        nadirantennadataHG2 = filename;
        name=append('Nadir HG2:',nadirantennadataHG2);
        txt.Text=name;
         
    end
    
    function MenuSelected_Zenithhg1(txt)
        global data_directory;
        global zenithantennadataHG1;
        path_zenith = data_directory;
        titleDataFileSelection = 'Select Zenith Antenna Pattern File';
        [filename,path_zenith] = uigetfile({'*.xlsx';'*.mat';'*.nc'}, titleDataFileSelection,(path_zenith));
        zenithantennadataHG1 = filename;
        name=append('Zenith HG1:',zenithantennadataHG1);
        txt.Text=name; 
    end
    function MenuSelected_Zenithhg2(txt)
        global data_directory;
        global zenithantennadataHG2;
        path_zenith = data_directory;
        titleDataFileSelection = 'Select Zenith Antenna Pattern File';
        [filename,path_zenith] = uigetfile({'*.xlsx';'*.mat';'*.nc'}, titleDataFileSelection,(path_zenith));
        zenithantennadataHG2 = filename;
        name=append('Zenith HG2:',zenithantennadataHG2);
        txt.Text=name; 
    end
    function loadinputfilespath(txt)
        global data_directory;
        title = 'Select the Auxiliary data directory';
        path ='../';
        data_directory = uigetdir(path,title);
        name_aux=append('AuxiliaryPath:',data_directory);
        txt.Text=name_aux;
    end
    
    function loadoutfilespath(txt)
        global outfile_directory;
        title = 'Select Output data directory';
        path ='../';
        outfile_directory = uigetdir(path,title);
        name_out=append('OutputPath:',outfile_directory);
        txt.Text=name_out;
    end

    function textChangedoss(txt,number_inside_array)
%         global Sampling; %only in function
%         if bg1.SelectedObject.Text ~= "Global" || bg1.SelectedObject.Text~= "Monte Carlo"
%             uiwait(msgbox({'under sampling is only for monte carlo or global mode'}))
%         end
        Sampling{number_inside_array}=txt.Value;
        if Sampling{1} || Sampling{2} ~= 0
         cb1_pn7.Value=0; 
        Sampling_all={0,0};
        Sampling_all_Text{1}=0;
        end
        if Sampling{3} ~=0
          cb2_pn7.Value=0;
          Sampling_all_Text{2}=0;
        end
    end
    function change_day_sp(txt,txt1,txt2,txt3,number_inside_array)
%         global Sampling_all;
%         global Sampling; %only in function

        Sampling_all{number_inside_array}=txt.Value;

          if Sampling_all{1}==1

              txt1.Value=0;
              txt2.Value=0;
             Sampling{1}=0;
             Sampling{2}=0;
         end
         if Sampling_all{2}==1
             txt3.Value=0;
             Sampling{3}=0;
         end
      if cb1_pn7.Value==1
          Sampling_all_Text{1}=cb1_pn7.Text;
      else
          Sampling_all_Text{1}=0;
      end

      if cb2_pn7.Value==1
          Sampling_all_Text{2}=cb2_pn7.Text;
      else
          Sampling_all_Text{2}=0;
      end
    end

    function textChangedtx(txt,number_inside_array)
%        global TxAnswers; %only in function
        TxAnswers(number_inside_array)=txt.Value;

    end

    function textChangedrx(txt,number_inside_array)
%          global RxAnswers; %only in function
        RxAnswers(number_inside_array)=txt.Value;
    end
    function textChangedss(txt,number_inside_array)
%          global SimulationSetting; %only in function
         SimulationSetting(number_inside_array)=txt.Value;
    end
    function textChangedsn(txt,number_inside_array)
%          global SensorNoise; %only in function
        SensorNoise(number_inside_array) = txt.Value;
    end
    function textchangedrunid(txt)
        global RUNID;
        RUNID = txt.Value;
    end
    
    function outputplotchange(txt,number)
%          global OutputPlot;
%          global OutputPlotText;  %only in function
        OutputPlot(number)=txt.Value;
        if txt.Value==1
          OutputPlotText{number}=txt.Text;
        else
          OutputPlotText{number}=0;
        end
    end

    function dropdownmenuchanged_soilmoisture(txt)
        global soilmoisture;
        global forest;
        global freeze_thaw;
        global wetlands;
        soilmoisture = get(txt,'Value');
        dd1_pn.Value='Not Selected';
        forest = dd1_pn.Value;
        dd2_pn.Value='Not Selected';
        freeze_thaw=dd2_pn.Value;
        dd3_pn.Value='Not Selected';
        wetlands=dd3_pn.Value;

    end
    function dropdownmenuchanged_forest(txt)
        global soilmoisture;
        global forest;
        global freeze_thaw;
        global wetlands;

        forest = get(txt,'Value');
        dd_pn.Value='Not Selected';
        soilmoisture=dd_pn.Value;
        dd2_pn.Value='Not Selected';
        freeze_thaw=dd2_pn.Value;
        dd3_pn.Value='Not Selected';
        wetlands=dd3_pn.Value;
    end
    function dropdownmenuchanged_freeze_thaw(txt)
        global soilmoisture;
        global forest;
        global freeze_thaw;
        global wetlands;
        freeze_thaw = get(txt,'Value');
        dd1_pn.Value='Not Selected';
        forest=dd1_pn.Value;
        dd_pn.Value='Not Selected';
        soilmoisture=dd_pn.Value;
        dd3_pn.Value='Not Selected';
        wetlands=dd3_pn.Value;
    end
    
    function dropdownmenuchanged_wetlands(txt)
        global soilmoisture;
        global forest;
        global freeze_thaw;
        global wetlands;
        wetlands = get(txt,'Value');
         dd1_pn.Value='Not Selected';
         forest=dd1_pn.Value;
         dd2_pn.Value='Not Selected';
         freeze_thaw=dd2_pn.Value;
         dd_pn.Value='Not Selected';
         soilmoisture=dd_pn.Value;

    end
    function changedmodeoperation(txt,number)
%          global OperationMode;
%             global mode0; %only in function
          mode0{number}=txt.Value;

       if txt.Value==1
            OperationMode{number}=txt.Text;
        else
            OperationMode{number}=0;
        end
    end
  
    function satelliteselection(txt)
%        global satellite
        satellite = get(txt.SelectedObject,'Text');

    end

    function signalselection(txt)
%         global signal;
        signal = get(txt.SelectedObject,'Text');
    end
    function  modeselection(txt)
%         global mode1;
        mode1 = get(txt.SelectedObject,'Text');
    end
    function slidervaluechange(txt)
%         global degreeofcoherence;

        degreeofcoherence = get(txt,'Value');
    end
    function Tx_switchvalue(sw)
%         global Txswitch;
        Txswitch = sw.Value;
    end
    
    function Rx_switchvalue_1(sw)
%         global Rxswitch_Nadir;
        Rxswitch_Nadir = sw.Value;
    end
    function Rx_switchvalue_2(sw)
%         global Rxswitch_Zenith;
        Rxswitch_Zenith = sw.Value;
    end
    
      function saveduserinput(fig)

%         global mode1;
%         global signal;
        global RUNID;
        global orbitdataHG1;
        global orbitdataHG2;
        global data_directory;
        global outfile_directory;
        global soilmoisture;
        global forest;
        global freeze_thaw;
        global wetlands;
        global nadirantennadataHG1;
        global nadirantennadataHG2;
        global zenithantennadataHG1;
        global zenithantennadataHG2;
%         global satellite;
        filename_input =sprintf('%s_%s.mat',RUNID,datestr(now,'ddmmyyyyHHMMSS'));
        if tag.exe == 1
            path_input = '../log/';
        elseif tag.exe ==0
            path_input = '../log/';
        end
        file_input = fullfile(path_input,filename_input);
                %added by Laura
        presentLogInputsFilename=file_input;
        %
        if tag.exe == 1
           save('../conf/UserInputs.mat','SimulationSetting','TxAnswers','RxAnswers','SensorNoise','Rxswitch_Zenith',...
                    'Rxswitch_Nadir','Txswitch','degreeofcoherence','signal','mode1','mode0','OperationMode','RUNID',...
                    'orbitdataHG1','orbitdataHG2','nadirantennadataHG1','nadirantennadataHG2','zenithantennadataHG1','zenithantennadataHG2','data_directory','outfile_directory',...
                    'soilmoisture','forest','freeze_thaw','wetlands','OutputPlot','OutputPlotText','Sampling','Sampling_all',...
                    'Sampling_all_Text','bioGeoInputs','bioGeoInputsVariable','satellite','presentLogInputsFilename','-mat');
        elseif tag.exe ==0
           save('../conf/UserInputs.mat','SimulationSetting','TxAnswers','RxAnswers','SensorNoise','Rxswitch_Zenith',...
                    'Rxswitch_Nadir','Txswitch','degreeofcoherence','signal','mode1','mode0','OperationMode','RUNID',...
                    'orbitdataHG1','orbitdataHG2','nadirantennadataHG1','nadirantennadataHG2','zenithantennadataHG1','zenithantennadataHG2','data_directory','outfile_directory',...
                    'soilmoisture','forest','freeze_thaw','wetlands','OutputPlot','OutputPlotText','Sampling','Sampling_all',...
                    'Sampling_all_Text','bioGeoInputs','bioGeoInputsVariable','satellite','presentLogInputsFilename','-mat');
        end

        save(file_input,'SimulationSetting','TxAnswers','RxAnswers','SensorNoise','Rxswitch_Zenith',...
                    'Rxswitch_Nadir','Txswitch','degreeofcoherence','signal','mode1','mode0','OperationMode','RUNID',...
                    'orbitdataHG1','orbitdataHG2','nadirantennadataHG1','nadirantennadataHG2','zenithantennadataHG1','zenithantennadataHG2','data_directory','outfile_directory',...
                    'soilmoisture','forest','freeze_thaw','wetlands','OutputPlot','OutputPlotText','Sampling','Sampling_all',...
                    'Sampling_all_Text','bioGeoInputs','bioGeoInputsVariable','satellite','presentLogInputsFilename','-mat');
        

        group = 'Exit';
        pref =  'Exit';
        title = "ContinueandExit";
        quest = {'Are you sure you want to continue?','Recheck entered Values.'};
        pbtns = {'Yes','No'};
        [pval,tf] = uigetpref(group,pref,title,quest,pbtns);
        drawnow;
        if pval == "yes"
            close(fig)
            closeIniGUI=0;
        end
    end
    function loadinputfile(~)
        global inputprevious;
%         global mode1;
%         global signal;
        global RUNID;
        global orbitdataHG1;
        global orbitdataHG2;
        global data_directory;
        global outfile_directory;
        global soilmoisture;
        global forest;
        global freeze_thaw;
        global wetlands;
        global nadirantennadataHG1;
        global zenithantennadataHG1;
        global nadirantennadataHG2;
        global zenithantennadataHG2;
        if tag.exe == 1
            Path = '../log/';
        elseif tag.exe ==0
             Path = '../log/';
        end
        file = 'values.mat';
        titleDataFileSelection   = 'Select previous input file';
        [file,Path] = uigetfile({'*.mat'}, titleDataFileSelection,[Path file]);
        file=fullfile(Path,file);
        inputprevious = load(file);
       tag_1(1).Value =inputprevious.SimulationSetting(1);
        SimulationSetting(1,1) = tag_1(1).Value;
%         tag_1(2).Value = inputprevious_last .SimulationSetting(2);
%         SimulationSetting(1,2) = tag_1(2).Value;
        tag_1(2).Value = inputprevious.SimulationSetting(2);
        SimulationSetting(1,2) = tag_1(2).Value;
        tag_1(3).Value = inputprevious.SimulationSetting(3);
        SimulationSetting(1,3) = tag_1(3).Value;
        tag_1(4).Value = inputprevious.SimulationSetting(5);
        SimulationSetting(1,5) = tag_1(4).Value;
        tag_1(5).Value = inputprevious.SimulationSetting(6);
        SimulationSetting(1,6) = tag_1(5).Value;
        tag_1(6).Value = inputprevious.SimulationSetting(7);
        SimulationSetting(1,7) = tag_1(6).Value;
        tag_1(7).Value = inputprevious.SimulationSetting(8);
        SimulationSetting(1,8) = tag_1(7).Value;
        tag_1(67).Value = inputprevious.SimulationSetting(9);
        SimulationSetting(1,9) = tag_1(67).Value;
        tag_1(68).Value = inputprevious.SimulationSetting(4);
        SimulationSetting(1,4) = tag_1(68).Value;
%         tag_1(69).Value = inputprevious.SimulationSetting(10);
%         SimulationSetting(1,10) = tag_1(69).Value;
        
        tag_1(8).Value = inputprevious.TxAnswers(1);
        TxAnswers(1,1) = tag_1(8).Value;
        tag_1(9).Value = inputprevious.TxAnswers(2);
        TxAnswers(1,2) = tag_1(9).Value;
        tag_1(10).Value = inputprevious.TxAnswers(3);
        TxAnswers(1,3) = tag_1(10).Value;
        tag_1(11).Value = inputprevious.TxAnswers(4);
        TxAnswers(1,4) = tag_1(11).Value;
        tag_1(12).Value =inputprevious.TxAnswers(5);
        TxAnswers(1,5) = tag_1(12).Value;
        tag_1(13).Value = inputprevious.TxAnswers(6);
        TxAnswers(1,6) = tag_1(13).Value;
        tag_1(14).Value =inputprevious.TxAnswers(7);
        TxAnswers(1,7) = tag_1(14).Value;
        
        tag_1(15).Value = inputprevious.RxAnswers(1);
        RxAnswers(1,1) = tag_1(15).Value;
        tag_1(16).Value = inputprevious.RxAnswers(2);
        RxAnswers(1,2) = tag_1(16).Value;
        tag_1(17).Value = inputprevious.RxAnswers(3);
        RxAnswers(1,3) = tag_1(17).Value;
        tag_1(18).Value = inputprevious.RxAnswers(4);
        RxAnswers(1,4) = tag_1(18).Value;
        tag_1(19).Value = inputprevious.RxAnswers(5);
        RxAnswers(1,5) = tag_1(19).Value;
        tag_1(20).Value = inputprevious.RxAnswers(6);
        RxAnswers(1,6) = tag_1(20).Value;
     
        tag_1(21).Value = inputprevious.SensorNoise(1);
        SensorNoise(1,1) = tag_1(21).Value;
        tag_1(22).Value =  inputprevious.SensorNoise(2);
        SensorNoise(1,2) = tag_1(22).Value;
        tag_1(23).Value =  inputprevious.SensorNoise(3);
        SensorNoise(1,3) = tag_1(23).Value;
        tag_1(24).Value =  inputprevious.SensorNoise(4);
        SensorNoise(1,4) = tag_1(24).Value;
        tag_1(65).Value=inputprevious.SensorNoise(5);
        SensorNoise(1,5)=tag_1(65).Value;
        tag_1(25).Value = inputprevious.SensorNoise(6);
        SensorNoise(1,6) = tag_1(25).Value;
        tag_1(66).Value=inputprevious.SensorNoise(7);
        SensorNoise(1,7)=tag_1(66).Value;
   
        if inputprevious.mode1 == "Single Point"
            tag_1(26).SelectedObject = tag_1(32);
            mode1 = "Single Point";
        end
        
        if inputprevious.mode1 == "Trackwise"
            tag_1(26).SelectedObject= tag_1(33);
            mode1 = "Trackwise" ;
        end
        
        if inputprevious.mode1 == "Monte Carlo"
            tag_1(26).SelectedObject= tag_1(34);
            mode1 = "Monte Carlo";
        end
        
        if inputprevious.mode1 == "Global"
            tag_1(26).SelectedObject=tag_1(35);
            mode1 = "Global";
        end

        if inputprevious.signal == "GPS L1"
            tag_1(27).SelectedObject = tag_1(36);
            signal = "GPS L1";
        end
        
        if inputprevious.signal == "GPS L1&L5"
            tag_1(27).SelectedObject = tag_1(37);
            signal = "GPS L1&L5";
        end
        
        if inputprevious.signal == "Galileo E1c"
            tag_1(27).SelectedObject= tag_1(38);
            signal = "Galileo E1c";
        end
        
        if inputprevious.signal == "Galileo E1c&E5a"
            tag_1(27).SelectedObject= tag_1(39);
            signal = "Galileo E1c&E5a";
        end
        
        tag_1(28).Value = inputprevious.Txswitch;
        Txswitch = tag_1(28).Value;
        tag_1(29).Value = inputprevious.Rxswitch_Nadir;
        Rxswitch_Nadir = tag_1(29).Value;
        tag_1(30).Value = inputprevious.Rxswitch_Zenith;
        Rxswitch_Zenith = tag_1(30).Value;
        tag_1(31).Value=inputprevious.degreeofcoherence;
        degreeofcoherence = tag_1(31).Value;
        tag_1(40).Value=inputprevious.mode0{1};
        mode0{1}=tag_1(40).Value;
        if tag_1(40).Value==1
            OperationMode{1}=tag_1(40).Text;
        else 
            OperationMode{1}=0;
        end

         tag_1(41).Value=inputprevious.mode0{2};
        mode0{2}=tag_1(41).Value;
        if tag_1(41).Value==1
            OperationMode{2}=tag_1(41).Text;
        else 
            OperationMode{2}=0;
        end
         tag_1(42).Value=inputprevious.mode0{3};
        mode0{3}=tag_1(42).Value;
        if tag_1(43).Value==1
            OperationMode{3}=tag_1(42).Text;
        else 
            OperationMode{3}=0;
        end
        tag_1(43).Value = inputprevious.RUNID;
        RUNID = tag_1(43).Value;
        orbitdataHG1 =inputprevious.orbitdataHG1 ;
        name_orbit=append('OrbitFile HG1:',orbitdataHG1);
        tag_1(50).Text=name_orbit;
        orbitdataHG2 = inputprevious.orbitdataHG2;
        name_meta=append('OrbitFile HG2:',orbitdataHG2);
        tag_1(51).Text=name_meta;
        zenithantennadataHG1=inputprevious.zenithantennadataHG1;
        name_zenith=append('Zenith HG1:',zenithantennadataHG1);
        tag_1(49).Text=name_zenith;
        nadirantennadataHG1 =inputprevious.nadirantennadataHG1;
        name_nadir=append('Nadir HG1:',nadirantennadataHG1);
        tag_1(48).Text=name_nadir;
        zenithantennadataHG2=inputprevious.zenithantennadataHG2;
        name_zenith=append('Zenith HG2:',zenithantennadataHG2);
        tag_1(72).Text=name_zenith;
        nadirantennadataHG2 =inputprevious.nadirantennadataHG2;
        name_nadir=append('Nadir HG2:',nadirantennadataHG2);
        tag_1(73).Text=name_nadir;
        data_directory = inputprevious.data_directory;
        name_aux=append('AuxiliaryPath:',data_directory);
        tag_1(52).Text=name_aux;
        outfile_directory =inputprevious.outfile_directory;
        name_out=append('OutputPath:',outfile_directory);
        tag_1(53).Text=name_out;
        tag_1(44).Value = inputprevious.soilmoisture;
        soilmoisture = tag_1(44).Value;
        tag_1(45).Value = inputprevious.forest;
        forest = tag_1(45).Value;
        tag_1(46).Value = inputprevious.freeze_thaw;
        freeze_thaw = tag_1(46).Value;
        tag_1(47).Value = inputprevious.wetlands;
        wetlands = tag_1(47).Value;
        tag_1(54).Value=inputprevious.OutputPlot(1);
        OutputPlot(1)=tag_1(54).Value;
        if tag_1(54).Value==1
          OutputPlotText{1}=tag_1(54).Text;
        else 
            OutputPlotText{1}=0;
        end
        tag_1(55).Value=inputprevious.OutputPlot(2);
        OutputPlot(2)=tag_1(55).Value;
          if tag_1(55).Value==1
          OutputPlotText{2}=tag_1(55).Text;
          else
          OutputPlotText{2}=0;
          end
         tag_1(56).Value=inputprevious.OutputPlot(3);
         OutputPlot(3)=tag_1(56).Value;
         if tag_1(56).Value==1
            OutputPlotText{3}=tag_1(56).Text;
         else 
             OutputPlotText{3}=0;
         end
         tag_1(57).Value=inputprevious.OutputPlot(4);
        OutputPlot(4)=tag_1(57).Value;
        if tag_1(57).Value==1   
             OutputPlotText{4}=tag_1(57).Text;
        else 
             OutputPlotText{4}=0;
        end
         tag_1(58).Value=inputprevious.OutputPlot(5);
         OutputPlot(5)=tag_1(58).Value;
        if tag_1(58).Value==1   
              OutputPlotText{5}=tag_1(58).Text;
        else 
             OutputPlotText{5}=0;
        end
        tag_1(59).Value=inputprevious.OutputPlot(6);
        OutputPlot(6)=tag_1(59).Value;
        if tag_1(59).Value==1   
              OutputPlotText{6}=tag_1(59).Text;
        else 
             OutputPlotText{6}=0;
        end
       tag_1(60).Value=inputprevious.Sampling{1};
       Sampling{1}=tag_1(60).Value;
        tag_1(61).Value=inputprevious.Sampling{2};
       Sampling{2}=tag_1(61).Value;
        tag_1(62).Value=inputprevious.Sampling{3};
       Sampling{3}=tag_1(62).Value;
       tag_1(63).Value=inputprevious.Sampling_all{1};
       Sampling_all{1}=tag_1(63).Value;
        if tag_1(63).Value==1
             Sampling_all_Text{1}=tag_1(63).Text;
        else
             Sampling_all_Text{1}=0;
        end
       tag_1(64).Value=inputprevious.Sampling_all{2};
       Sampling_all{2}=tag_1(64).Value;
        if tag_1(64).Value==1
                Sampling_all_Text{2}=tag_1(64).Text;
        else
             Sampling_all_Text{2}=0;
        end
        if inputprevious.satellite == "HydroGNSS-1"
            tag_1(69).SelectedObject = tag_1(70);
            satellite = "HydroGNSS-1";
        end
        if inputprevious.satellite == "HydroGNSS-2"
            tag_1(69).SelectedObject = tag_1(71);
            satellite = "HydroGNSS-2";
        end
       %added by Laura
       if isfield(inputprevious, 'bioGeoInputs')==0
            inputprevious.bioGeoInputs=bioGeoInputs;
            inputprevious.bioGeoInputsVariable=bioGeoInputsVariable;
       end
       bioGeoInputs=inputprevious.bioGeoInputs;
       bioGeoInputsVariable=inputprevious.bioGeoInputsVariable;
       %
    end

    function loadinputfile_previous(~)
        global inputprevious_last;
%         global mode1;
%         global signal;
        global RUNID;
        global orbitdataHG1;
        global orbitdataHG2;
        global data_directory;
        global outfile_directory;
        global soilmoisture;
        global forest;
        global freeze_thaw;
        global wetlands;
        global nadirantennadataHG1;
        global zenithantennadataHG1;
        global nadirantennadataHG2;
        global zenithantennadataHG2;

        if tag.exe == 1
                Path = '../conf/';
        elseif tag.exe ==0
                Path = '../conf/';
        end

        filename = 'UserInputs.mat';
        file=fullfile(Path,filename);     
        inputprevious_last = load(file);

        tag_1(1).Value =inputprevious_last .SimulationSetting(1);
        SimulationSetting(1,1) = tag_1(1).Value;
%         tag_1(2).Value = inputprevious_last .SimulationSetting(2);
%         SimulationSetting(1,2) = tag_1(2).Value;
        tag_1(2).Value = inputprevious_last .SimulationSetting(2);
        SimulationSetting(1,2) = tag_1(2).Value;
        tag_1(3).Value = inputprevious_last .SimulationSetting(3);
        SimulationSetting(1,3) = tag_1(3).Value;
        tag_1(4).Value = inputprevious_last .SimulationSetting(5);
        SimulationSetting(1,5) = tag_1(4).Value;
        tag_1(5).Value = inputprevious_last .SimulationSetting(6);
        SimulationSetting(1,6) = tag_1(5).Value;
        tag_1(6).Value = inputprevious_last .SimulationSetting(7);
        SimulationSetting(1,7) = tag_1(6).Value;
        tag_1(7).Value = inputprevious_last .SimulationSetting(8);
        SimulationSetting(1,8) = tag_1(7).Value;
        tag_1(67).Value = inputprevious_last .SimulationSetting(9);
        SimulationSetting(1,9) = tag_1(67).Value;
        tag_1(68).Value = inputprevious_last.SimulationSetting(4);
        SimulationSetting(1,4) = tag_1(68).Value;
%         tag_1(69).Value = inputprevious_last .SimulationSetting(10);
%         SimulationSetting(1,10) = tag_1(69).Value;
        
        tag_1(8).Value = inputprevious_last.TxAnswers(1);
        TxAnswers(1,1) = tag_1(8).Value;
        tag_1(9).Value = inputprevious_last.TxAnswers(2);
        TxAnswers(1,2) = tag_1(9).Value;
        tag_1(10).Value = inputprevious_last.TxAnswers(3);
        TxAnswers(1,3) = tag_1(10).Value;
        tag_1(11).Value = inputprevious_last.TxAnswers(4);
        TxAnswers(1,4) = tag_1(11).Value;
        tag_1(12).Value =inputprevious_last.TxAnswers(5);
        TxAnswers(1,5) = tag_1(12).Value;
        tag_1(13).Value = inputprevious_last.TxAnswers(6);
        TxAnswers(1,6) = tag_1(13).Value;
        tag_1(14).Value =inputprevious_last.TxAnswers(7);
        TxAnswers(1,7) = tag_1(14).Value;
        
        tag_1(15).Value = inputprevious_last.RxAnswers(1);
        RxAnswers(1,1) = tag_1(15).Value;
        tag_1(16).Value = inputprevious_last.RxAnswers(2);
        RxAnswers(1,2) = tag_1(16).Value;
        tag_1(17).Value = inputprevious_last.RxAnswers(3);
        RxAnswers(1,3) = tag_1(17).Value;
        tag_1(18).Value = inputprevious_last.RxAnswers(4);
        RxAnswers(1,4) = tag_1(18).Value;
        tag_1(19).Value = inputprevious_last.RxAnswers(5);
        RxAnswers(1,5) = tag_1(19).Value;
        tag_1(20).Value = inputprevious_last.RxAnswers(6);
        RxAnswers(1,6) = tag_1(20).Value;
     
        tag_1(21).Value = inputprevious_last.SensorNoise(1);
        SensorNoise(1,1) = tag_1(21).Value;
        tag_1(22).Value =  inputprevious_last.SensorNoise(2);
        SensorNoise(1,2) = tag_1(22).Value;
        tag_1(23).Value =  inputprevious_last.SensorNoise(3);
        SensorNoise(1,3) = tag_1(23).Value;
        tag_1(24).Value =  inputprevious_last.SensorNoise(4);
        SensorNoise(1,4) = tag_1(24).Value;
        tag_1(65).Value=inputprevious_last.SensorNoise(5);
        SensorNoise(1,5)=tag_1(65).Value;
        tag_1(25).Value = inputprevious_last.SensorNoise(6);
        SensorNoise(1,6) = tag_1(25).Value;
        tag_1(66).Value=inputprevious_last.SensorNoise(7);
        SensorNoise(1,7)=tag_1(66).Value;
   
        if inputprevious_last.mode1 == "Single Point"
            tag_1(26).SelectedObject = tag_1(32);
            mode1 = "Single Point";
        end
        
        if inputprevious_last.mode1 == "Trackwise"
            tag_1(26).SelectedObject= tag_1(33);
            mode1 = "Trackwise" ;
        end
        
        if inputprevious_last.mode1 == "Monte Carlo"
            tag_1(26).SelectedObject= tag_1(34);
            mode1 = "Monte Carlo";
        end
        
        if inputprevious_last.mode1 == "Global"
            tag_1(26).SelectedObject=tag_1(35);
            mode1 = "Global";
        end

        if inputprevious_last.signal == "GPS L1"
            tag_1(27).SelectedObject = tag_1(36);
            signal = "GPS L1";
        end
        
        if inputprevious_last.signal == "GPS L1&L5"
            tag_1(27).SelectedObject = tag_1(37);
            signal = "GPS L1&L5";
        end
        
        if inputprevious_last.signal == "Galileo E1c"
            tag_1(27).SelectedObject= tag_1(38);
            signal = "Galileo E1c";
        end
        
        if inputprevious_last.signal == "Galileo E1c&E5a"
            tag_1(27).SelectedObject= tag_1(39);
            signal = "Galileo E1c&E5a";
        end
        
        tag_1(28).Value = inputprevious_last.Txswitch;
        Txswitch = tag_1(28).Value;
        tag_1(29).Value = inputprevious_last.Rxswitch_Nadir;
        Rxswitch_Nadir = tag_1(29).Value;
        tag_1(30).Value = inputprevious_last.Rxswitch_Zenith;
        Rxswitch_Zenith = tag_1(30).Value;
        tag_1(31).Value=inputprevious_last.degreeofcoherence;
        degreeofcoherence = tag_1(31).Value;
        tag_1(40).Value=inputprevious_last.mode0{1};
        mode0{1}=tag_1(40).Value;
        if tag_1(40).Value==1
            OperationMode{1}=tag_1(40).Text;
        else 
            OperationMode{1}=0;
        end

         tag_1(41).Value=inputprevious_last.mode0{2};
        mode0{2}=tag_1(41).Value;
        if tag_1(41).Value==1
            OperationMode{2}=tag_1(41).Text;
        else 
            OperationMode{2}=0;
        end
         tag_1(42).Value=inputprevious_last.mode0{3};
        mode0{3}=tag_1(42).Value;
        if tag_1(43).Value==1
            OperationMode{3}=tag_1(42).Text;
        else 
            OperationMode{3}=0;
        end
        tag_1(43).Value = inputprevious_last.RUNID;
        RUNID = tag_1(43).Value;
        orbitdataHG1 =inputprevious_last.orbitdataHG1 ;
        name_orbit=append('OrbitFile HG1:',orbitdataHG1);
        tag_1(50).Text=name_orbit;
        orbitdataHG2 = inputprevious_last.orbitdataHG2;
        name_meta=append('OrbitFile HG2:',orbitdataHG2);
        tag_1(51).Text=name_meta;
        zenithantennadataHG1=inputprevious_last.zenithantennadataHG1;
        name_zenith=append('Zenith HG1:',zenithantennadataHG1);
        tag_1(49).Text=name_zenith;
        nadirantennadataHG1 =inputprevious_last.nadirantennadataHG1;
        name_nadir=append('Nadir HG1:',nadirantennadataHG1);
        tag_1(48).Text=name_nadir;
        zenithantennadataHG2=inputprevious_last.zenithantennadataHG2;
        name_zenith=append('Zenith HG2:',zenithantennadataHG2);
        tag_1(72).Text=name_zenith;
        nadirantennadataHG2 =inputprevious_last.nadirantennadataHG2;
        name_nadir=append('Nadir HG2:',nadirantennadataHG2);
        tag_1(73).Text=name_nadir;
        data_directory = inputprevious_last.data_directory;
        name_aux=append('AuxiliaryPath:',data_directory);
        tag_1(52).Text=name_aux;
        outfile_directory =inputprevious_last.outfile_directory;
        name_out=append('OutputPath:',outfile_directory);
        tag_1(53).Text=name_out;
        tag_1(44).Value = inputprevious_last.soilmoisture;
        soilmoisture = tag_1(44).Value;
        tag_1(45).Value = inputprevious_last.forest;
        forest = tag_1(45).Value;
        tag_1(46).Value = inputprevious_last.freeze_thaw;
        freeze_thaw = tag_1(46).Value;
        tag_1(47).Value = inputprevious_last.wetlands;
        wetlands = tag_1(47).Value;
        tag_1(54).Value=inputprevious_last.OutputPlot(1);
        OutputPlot(1)=tag_1(54).Value;
        if tag_1(54).Value==1
          OutputPlotText{1}=tag_1(54).Text;
        else 
            OutputPlotText{1}=0;
        end
        tag_1(55).Value=inputprevious_last.OutputPlot(2);
        OutputPlot(2)=tag_1(55).Value;
          if tag_1(55).Value==1
          OutputPlotText{2}=tag_1(55).Text;
          else
          OutputPlotText{2}=0;
          end
         tag_1(56).Value=inputprevious_last.OutputPlot(3);
         OutputPlot(3)=tag_1(56).Value;
         if tag_1(56).Value==1
            OutputPlotText{3}=tag_1(56).Text;
         else 
             OutputPlotText{3}=0;
         end
         tag_1(57).Value=inputprevious_last.OutputPlot(4);
        OutputPlot(4)=tag_1(57).Value;
        if tag_1(57).Value==1   
             OutputPlotText{4}=tag_1(57).Text;
        else 
             OutputPlotText{4}=0;
        end
         tag_1(58).Value=inputprevious_last.OutputPlot(5);
         OutputPlot(5)=tag_1(58).Value;
        if tag_1(58).Value==1   
              OutputPlotText{5}=tag_1(58).Text;
        else 
             OutputPlotText{5}=0;
        end
        tag_1(59).Value=inputprevious_last.OutputPlot(6);
        OutputPlot(6)=tag_1(59).Value;
        if tag_1(59).Value==1   
              OutputPlotText{6}=tag_1(59).Text;
        else 
             OutputPlotText{6}=0;
        end
       tag_1(60).Value=inputprevious_last.Sampling{1};
       Sampling{1}=tag_1(60).Value;
        tag_1(61).Value=inputprevious_last.Sampling{2};
       Sampling{2}=tag_1(61).Value;
        tag_1(62).Value=inputprevious_last.Sampling{3};
       Sampling{3}=tag_1(62).Value;
       tag_1(63).Value=inputprevious_last.Sampling_all{1};
       Sampling_all{1}=tag_1(63).Value;
        if tag_1(63).Value==1
             Sampling_all_Text{1}=tag_1(63).Text;
        else
             Sampling_all_Text{1}=0;
        end
       tag_1(64).Value=inputprevious_last.Sampling_all{2};
       Sampling_all{2}=tag_1(64).Value;
        if tag_1(64).Value==1
                Sampling_all_Text{2}=tag_1(64).Text;
        else
             Sampling_all_Text{2}=0;
        end
        if inputprevious_last.satellite == "HydroGNSS-1"
            tag_1(69).SelectedObject = tag_1(70);
            satellite = "HydroGNSS-1";
        end
        if inputprevious_last.satellite == "HydroGNSS-2"
            tag_1(69).SelectedObject = tag_1(71);
            satellite = "HydroGNSS-2";
        end

       %added by Laura
       if isfield(inputprevious_last, 'bioGeoInputs')==0
            inputprevious_last.bioGeoInputs=bioGeoInputs;
            inputprevious_last.bioGeoInputsVariable=bioGeoInputsVariable;
       end

       bioGeoInputs=inputprevious_last.bioGeoInputs;
       bioGeoInputsVariable=inputprevious_last.bioGeoInputsVariable;
       %

        drawnow;
    end
     end
