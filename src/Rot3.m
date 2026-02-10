% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% The Intellectual Property of the contents of this directory, code, data
% and files belong to SSTL and are @copyright SSTL 2018 unless otherwise
% specifically indicated within the code or file.
%
% The contents are being transferred to ESA under the agreement within
% Contract 4000122597//17/NL/CRS/hh - Design and development of a tool for 
% antenna carrier characterisation from flying space GNSS receiver
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% =========================================================================
% Description:This function computes a rotation matrix R
%
% INPUT: 
%    - angle: Greenwich Sideral Time hour angle (GST)
% 
% OUTPUT: 
%    - R: rotation matrix to go from ECEF to ECI frame
% =========================================================================

function R = Rot3(angle)
R = eye(3);
R(1, 1) =  cos(angle);
R(1, 2) =  sin(angle);
R(2, 1) = -sin(angle);
R(2, 2) =  cos(angle);
end %function