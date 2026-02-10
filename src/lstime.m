% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% The Intellectual Property of the contents of this directory, code, data
% and files belong to SSTL and are @copyright SSTL 2018 unless otherwise
% specifically indicated within the code or file.
%
% The contents are being transferred to ESA under the agreement within
% Contract 4000122597//17/NL/CRS/hh - Design and development of a tool for 
% antenna carrier characterisation from flying space GNSS receiver
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Function Description: 
% This function computes the thetaGMST in radians
% 
% INPUT: Julian Date 
% 
% OUTPUT: thetaGMST in radians
% -------------------------------------------------------------------------

function thetaGMST = lstime(jd)

tut1 = (jd - 2451545.0) / 36525.0;
theta = 67310.54841 + (((876600.0 * 3600.0) + 8640184.812866) * tut1) + (0.093104 * (tut1 * tut1)) - (6.2e-6 * (tut1 * tut1 * tut1));
thetaGMST = mod(theta, 86400.0) / 240.0;
thetaGMST = thetaGMST * pi / 180.0; %return thetaGMST in radians.