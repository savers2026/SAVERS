function [SP_lla, SP_ECEF] = Specular_Point_Position_Calculation(Pos_GPS, Pos_Rec)
%--------------------------------------------------------------------------
% [SP_lla, SP_ECEF] = Specular_Point_Position_Calculation(Pos_Rec, Pos_GPS)
%--------------------------------------------------------------------------
% Function Name:  Specular_Point_Position_Calculation.m
% Author:         Starlab Barcelona SL
% Last Update:    2010/10/21
%                 March 2023 by Laura Dente (the minor axis is computed
%                 from the major axis and the flattening)
%--------------------------------------------------------------------------
%
% The Specular Point Positions are calculated in latitude, longitude, and
% altitude as well as in ECEF coordinates.
% The positions of the specular points are calculated based on the receiver
% position and the GPS satellite position using the WGS84 ellipsoid model 
% The number of emitter and receiver coordinates should be the same. If the
% receiver is static, then input just one row in the receiver coordinates
% matrix
%
% Input
% -----
% Pos_Rec   -   [N x 3] ECEF coordinates of the emitter satellite (m)
% Pos_GPS   -   [N x 3] ECEF coordinates of the receiver satellite (m)
%
% Output
% ------
% SP_lla    -   [N x 3] ECEF coordinates of the specular point position
% SP_ECEF   -   [N x 3] lla coordinates of the specular point position
%
%--------------------------------------------------------------------------

[SP_lla,SP_ECEF,inds] = findSpecularPoint(Pos_GPS,Pos_Rec);
    

return;



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%------------------------       FUNCTIONS         ------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


function [R_lla,R_ECEF,inds] = findSpecularPoint(emitter,receiver,model)
%--------------------------------------------------------------------------
% [R_lla,R_ECEF,inds] = findSpecularPoint(emitter,receiver,model)
%--------------------------------------------------------------------------
% Function Name:  trueAnomaly.m
% Author:         Starlab Barcelona SL
% Last Update:    2010/10/21
%--------------------------------------------------------------------------
%
% Finds the specular reflection point on earth surface for a given emitter and 
% receiver location. Assumes an ellipsoidal Earth.
% This function might return empty arrays if no specular points are found.
% The inds array is the array of indices of the emitter/receiver arrays
% corresponding to the discovered points.
%
% Input
% -----
% emitter   -   [N x 3] ECEF coordinates of the emitter satellite (m)
% receiver  -   [N x 3] ECEF coordinates of the receiver satellite (m)
% model     -   WGS84 / JASON ellipsoid model
%
% Output
% ------
% R_lla     -   [N x 3] ECEF coordinates of the specular point position
% lon       -   [N x 3] lla coordinates of the specular point position
% inds      -   [N]     Indices of the emitter/receiver arrays corresponding
%                       to the points discovered
%--------------------------------------------------------------------------

% Set WGS84 ellipsoid as default
if nargin==2
    model='WGS84';
end

if prod(double(model=='WGS84'))
    earthMajax = 6.37813700000000e6;
    f = 1/298.257223563;
    earthMinax = (1-f)*earthMajax;
    %earthMinax = 6.356752314140356e6;
    %f = (earthMajax-earthMinax)/earthMajax;
else
    % JASON Ellipsoid
    earthMajax = 6.37813630000000e6;
    f          = 1/298.257;
    earthMinax = (1-f)*earthMajax;
end

% Check input array sizes etc.
se = size(emitter);
sr = size(receiver);

if se(2)~=3 || sr(2)~=3
	error('Input arrays must have 3 columns');
end

if sr(1)==1
    receiver = repmat(receiver,sr(2),1);
end

num = se(1);

% Now find specular point for each input row

% Make a first guess
v = (emitter(1,:)+receiver(1,:))/2;
%v = receiver(j,:);
p = rect2sphere(v);
firstGuess = [p(2) p(3)];
res = firstGuess;

lats=[];
longs=[];
alts=[];
inds=[];
r=[];

at = 1;

for j=1:num
    % Check if specular point is impossible
    if rayIntersectsEarth(emitter(j,:),receiver(j,:))
        if j<num-1
            v = (emitter(j+1,:)+receiver(1,:))/2;
            p = rect2sphere(v);
            firstGuess = [p(2) p(3)];
            res = firstGuess;
        else
            if isempty(lats)
                R_lla   = NaN;
                R_ECEF  = NaN;
                inds    = NaN;
            end
            return
        end
    else
    
        firstGuess = res;

        % Set optimisation options
        options = optimset('TolFun',10e-8,'TolX',10e-8,'Display','Off',...
            'Diagnostics','Off','MaxIter',10000,'MaxFunEvals',10000);

        % Search for the specular point - see cost function for more details
        [res,feval,exitflag,output] = fminsearch(@costFunction,firstGuess,...
            options,emitter(j,:),receiver(1,:));


        if exitflag>0
            if nargout==5
                inds(at)=j;
            end
            % Convert to latitude and longitude
            r = earthSurfacePoint(res(1),res(2));
            R(at,:)=r;
            %[lats(at) longs(at) alts(at)] = ECEFtoLatLonAlt(r(1),r(2),r(3));
            at = at+1;
        end
        
    end
    
end

% Convert to Latitude Longitude Altitude
[lats longs alts] = ECEFtoLatLonAlt(R(:,1),R(:,2),R(:,3));

R_lla   = [lats longs alts];
R_ECEF  = R;


if isempty(lats)
    R_lla   = NaN;
    R_ECEF  = NaN;
    inds    = NaN;
end

return


function [lat,lon,alt] = ECEFtoLatLonAlt(x,y,z)
%--------------------------------------------------------------------------
% [lat,lon,alt] = ECEF2LatLonAlt(x,y,z)
%--------------------------------------------------------------------------
% Function Name:  trueAnomaly.m
% Author:         Starlab Barcelona SL
% Last Update:    2010/10/21
%--------------------------------------------------------------------------
%
% Converts from ECEF cartesian coordinates to latitude, longitude and altitude
%
% Input
% -----
% x     -   ECEF x coordinate (m)
% y     -   ECEF y coordinate (m)
% z     -   ECEF z coordinate (m)
%
% Output
% ------
% lat   -   Latitude (deg)
% lon   -   Longitude (deg)
% alt   -   Altitude (m)
%--------------------------------------------------------------------------

% WGS84 Ellipsoid

if nargin==1
  a_=x;
  clear x
  z=a_(:,3)';
  y=a_(:,2)';
  x=a_(:,1)';
end

earthMajax = 6.37813700000000e6;
f = 1/298.257223563;
earthMinax = (1-f)*earthMajax;
%earthMinax = 6.356752314140356e6;
%f = 1-earthMinax/earthMajax;
tol = 10e-10;

if ~isempty(x(x==0))
  x(x==0)=1e-30;
end
if ~isempty(y(y==0))
  y(y==0)=1e-30;
end
if ~isempty(z(z==0))
  z(z==0)=1e-30;
end

lon=atan2(y,x)*180/pi;
r = sqrt(x.*x+y.*y+z.*z);
delta = atan2(z,sqrt(x.*x+y.*y));
phidk = delta;
while 1
    rck = earthMajax*sqrt((1-2*f+f*f)./(1-(2*f-f*f)*cos(phidk).^2));
    phik = atan(tan(phidk)./((1-f)*(1-f)));
    hk = sqrt(r.*r-rck.*rck.*sin(phik-phidk).^2)-rck.*cos(phik-phidk);
    phiddk = delta-asin(hk.*sin(phik-phidk)./r);
    if phiddk==0;keyboard;end

    if all(abs(phiddk-phidk)./abs(phiddk))<tol
        break;
    end
    phidk = phiddk;
end
alt= hk;
lat = 180*atan(tan(phiddk)./((1-f)*(1-f)))/pi;

return




function s = rect2sphere(p)
%--------------------------------------------------------------------------
% s = rect2sphere(p)
%--------------------------------------------------------------------------
% Function Name:  rect2sphere.m
% Author:         Starlab Barcelona SL
% Last Update:    2010/10/21
%--------------------------------------------------------------------------
%
% Conversts ECI rectangular coordinates to [radius ascension declination]
%
% Input
% -----
% p     -   ECI rectangular coordinates (m)
%
% Output
% ------
% s     -   Spherical coordinates (deg)
%--------------------------------------------------------------------------

s(1) = sqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3));
if s(1)==0
    s(2)=0;
    s(3)=0;
else
    r = sqrt(p(1)*p(1)+p(2)*p(2));
    s(3) = acos(r/s(1));
    if (p(3)<0)
        s(3)=-s(3);
    end
    if r==0
        s(2) = 0;
    else
        s(2) = acos(p(1)/r);
        if p(2)<0
            s(2) = 2*pi-s(2);
        end
    end
end

return



function x = sphere2rect(p)
%--------------------------------------------------------------------------
% x = sphere2rect(p)
%--------------------------------------------------------------------------
% Function Name:  sphere2rect.m
% Author:         Starlab Barcelona SL
% Last Update:    2010/10/21
%--------------------------------------------------------------------------
%
% Converts from spherical coordinates to rectangular
%
% Input
% -----
% p     -   [radius ascension declination] (in radians)
%
% Output
% ------
% x     -   [x y z] ECEF coordinates (m)
%--------------------------------------------------------------------------

x(1) = p(1)*cos(p(3))*cos(p(2));
x(2) = p(1)*cos(p(3))*sin(p(2));
x(3) = p(1)*sin(p(3));

return



function b = rayIntersectsEarth(p,q)
%--------------------------------------------------------------------------
% b = rayIntersectsEarth(p,q)
%--------------------------------------------------------------------------
% Function Name:  rayIntersectsEarth.m
% Author:         Starlab Barcelona SL
% Last Update:    2010/10/21
%--------------------------------------------------------------------------
%
% Returns 1 if the line segment from p to q (in ECI/ECEF rectangular coords)
% intersects with the earth ellipsoid
%
% Input
% -----
% p     -   [x y z] ECI/ECEF coordinates (m)
% q     -   [x y z] ECI/ECEF coordinates (m)
%
% Output
% ------
% b     -   [1] Boolean
%--------------------------------------------------------------------------

earthMajax = 6.37813700000000e6;
f = 1/298.257223563;
earthMinax = (1-f)*earthMajax;
%earthMinax = 6.356752314140356e6;
c = earthMajax/earthMinax;
c = c*c;
d = earthMajax*earthMajax;
i = q(1)-p(1);
j = q(2)-p(2);
k = q(3)-p(3);
cc = p(1)*p(1)+p(2)*p(2)+c*p(3)*p(3)-d;
bb = 2*p(1)*i + 2*p(2)*j+2*c*p(3)*k;
aa = i*i+j*j+c*k*k;
det = bb*bb-4*aa*cc;
b = 1;
if (det<=0)
    b = 0;
else
    % check if line segment actually passes through the earth
    sdet = sqrt(det);
    t1 = (-bb+sdet)/(2*aa);
    t2 = (-bb-sdet)/(2*aa);
    if t1<0 | t1>1 | t2<0 | t2>1
        b = 0;
    end
end

return




function c = costFunction(pos,emitter,receiver)
%--------------------------------------------------------------------------
% c = costFunction(pos,emitter,receiver)
%--------------------------------------------------------------------------
% Function Name:  costFunction.m
% Author:         Starlab Barcelona SL
% Last Update:    2010/10/21
%--------------------------------------------------------------------------
%
% Returns the cost value to be minimised in finding the specular point on
% the surface of the earth.
%
% Input
% -----
% pos       -   [ascension declination] in radians
% emitter   -   [x y z] rect. coords of emitter
% receiver  -   [x y z] rect. coords of receiver
%
% Output
% ------
% c         -   [1] Cost function value
%--------------------------------------------------------------------------

ascension=pos(1);
declination=pos(2);
% Get xyz of point on earth surface
x=earthSurfacePoint(ascension,declination);
d1 = emitter-x;
d2 = receiver-x;
c = norm(d1)+norm(d2);

return




function v = earthSurfacePoint(ascension,declination)
%--------------------------------------------------------------------------
% v = earthSurfacePoint(ascension,declination)
%--------------------------------------------------------------------------
% Function Name:  earthSurfacePoint.m
% Author:         Starlab Barcelona SL
% Last Update:    2010/10/21
%--------------------------------------------------------------------------
%
% Returns the ECI rectangular coords of the earth's surface point at the
% supplied ascension,declination (in radians)
%
% Input
% -----
% ascension   -   [1] in radians
% declination -   [1] in radians
%
% Output
% ------
% v           -   [1 x 3] ECEF Earth surface coordinates
%--------------------------------------------------------------------------

p = [1 ascension declination];
s = sphere2rect(p);
r = earthRadiusInDirection(s);
v = r*s;

return




function r = earthRadiusInDirection(v)
%--------------------------------------------------------------------------
% r = earthRadiusInDirection(v)
%--------------------------------------------------------------------------
% Function Name:  earthRadiusInDirection.m
% Author:         Starlab Barcelona SL
% Last Update:    2010/10/21
%--------------------------------------------------------------------------
%
% Returns the earth radius in the direction v assuming an ellipsoidal earth
%
% Input
% -----
% v           -   [1 x 3] direction from centre of the Earth in ECEF coord.
%
% Output
% ------
% r           -   [1] Earth radius
%--------------------------------------------------------------------------

earthMajax = 6.37813700000000e6;
f = 1/298.257223563;
earthMinax = (1-f)*earthMajax;
%earthMinax = 6.356752314140356e6;
r = norm(v)*earthMajax/sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3)*earthMajax*earthMajax/(earthMinax*earthMinax));

return
