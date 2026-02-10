
function [Rx, Ry, Rz] = rotation_matrices(phi, theta, psi)

%%%%Rotation matrices to determines the attitude

%Generate rotation matrix with angle alpha (radians), with axis x
%See Eq 13.5 in the MERRByS reference
Rx   = [      1                  0                       0;...
              0            cosd(phi)                sind(phi);...
              0           -sind(phi)                cosd(phi)];

%Generate rotation matrix with angle alpha (radians), with axis y
%See Eq 13.5 in the MERRByS reference          
Ry   = [cosd(theta)             0                -sind(theta);...
              0                 1                        0;...
        sind(theta)             0                 cosd(theta)];

    
%Generate rotation matrix with angle alpha (radians), with axis z
%See Eq 13.5 in the MERRByS reference    
Rz   = [ cosd(psi)         sind(psi)                    0;...
        -sind(psi)         cosd(psi)                    0;...
              0               0                         1];
        
end