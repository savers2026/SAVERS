
function [Amat_ECEF2OF, Amat_OF2BF, Amat_BFtoAF] = matrixFrameGen(pHydroG_ECEF, vHydroG_ECEF, nom_attitude_angles, attitude_angles, Euler_angles)


   %%The function takes as input column vector
   %%Double check and fixin format
   
   if isrow(pHydroG_ECEF), pHydroG_ECEF = pHydroG_ECEF.'; end
   if isrow(vHydroG_ECEF), vHydroG_ECEF = vHydroG_ECEF.'; end
   
   %% DEFINING ATTITUDE ANGLES
  
   %OLD
   %phi_roll        = attitude_angles(1);
   %theta_pitch     = attitude_angles(2);
   %psi_yaw         = attitude_angles(3);        
   
   %NEW: ADDING nominal attitude to eventual error
   phi_roll        = attitude_angles(1) + nom_attitude_angles(1);
   theta_pitch     = attitude_angles(2) + nom_attitude_angles(2);
   psi_yaw         = attitude_angles(3) + nom_attitude_angles(3);    

   %% AMIR: check height for exluding coriolis considration for lower than 15km height, for airborne it provides problem for antenna footprint projection
  [~, ~, pHydroG_Height] = xyz2llh(pHydroG_ECEF(1),pHydroG_ECEF(2),pHydroG_ECEF(3));
   if pHydroG_Height > 15e3
  %% DC: correzione velocità - equivalent to Coriolis
   vHydroG_ECEF  = vHydroG_ECEF + (cross([0,0,7.2921158553*1.e-5],pHydroG_ECEF)).';
   end

   %% Orbit reference system case 1 - Rius approach

   %unit vector z
   zhat_orbit = -pHydroG_ECEF./norm(pHydroG_ECEF);

   %unit vector y orthogonal to the velocity and toward left wrt the movement
   yhat_orbit = -cross(pHydroG_ECEF, vHydroG_ECEF)./norm(cross(pHydroG_ECEF, vHydroG_ECEF));
    
   %unit vector x parallel to velocity  
   xhat_orbit = cross(yhat_orbit, zhat_orbit);

   orbitFrame = [xhat_orbit.'; yhat_orbit.'; zhat_orbit.'];

   %matrix for conversion from ECEF to OF
   Amat_ECEF2OF = orbitFrame;

   %% Body Frame

   % nota bene: orbit diretto lungo la prua (ram, asse x del body), ma la prua non è necessariamente
   % lungo la velocità

   % ypr disallineano il body dall'orbit
   [Rx, Ry, Rz] = rotation_matrices(phi_roll, theta_pitch, psi_yaw);
       
   %% Matrix for conversion from OF to BF

    %OLD
    %rotation 213 - (theta, phi, psi) --> (pitch, roll, yaw)
    %(remind: read from right to left)
%     Amat_OF2BF = Rz*Rx*Ry;
%     Amat_OF2BF = Ry*Rx*Rz;

    %NEW
    %rotation 123 -> (phi, theta, psi) --> (roll, pitch, yaw) 
    %(remind: read from right to left)
    Amat_OF2BF = Rz*Ry*Rx;  

   %% Antenna Frame - if Euler_angles == 0 it is the same as body

   [RxE, RyE, RzE] = rotation_matrices(Euler_angles(1), Euler_angles(2), Euler_angles(3));

    %rotation 123 - yrp
    Amat_BFtoAF = RzE*RyE*RxE;

end