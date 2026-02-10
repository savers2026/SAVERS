
% % clear all; close all;
% theta=10; %deg
% phi=0;%deg
% theta_s=10; %deg
% phi_s=0;%deg
% f=1.41;  %GHz
% 
% t=deg2rad(theta);
% p=deg2rad(phi);
% ts=deg2rad(theta_s);
% ps=deg2rad(phi_s);
% SM=15; %[5 20 40];
% sstd=1; %[0.5 1.5 3];  %cm
% R0=22000000./cosd(theta); %m
% R2=1000./cosd(theta); %m
% betac=1; %[deg] %half beamwidth
% betac=deg2rad(betac);
% beta0rad=betac;
% beta2rad=betac;
% coeff=coerente(f*10^9,SM,sstd/100,R0,R2,beta0rad,beta2rad,t,p,ts,ps);
% db(coeff,'power')
function coeff=coerente(f,SM,sstd,R0,R2,beta0rad,beta2rad,t,p,ts,ps,g0,gs)

         ht=R0*cos(t);
         hr=R2*cos(ts);
cd=costantedielettrica(SM./100,f/10^9);
refr_index=sqrt(cd);

         if hr==ht && ps==pi+p && t==ts  % condition for the backward direction. In this case the specular angle is null

              angolo_speculare=0;

         else % in this case the specular angle has computed based on R(2.30) and R(2.31)
 															
      xt=ht*tan(t)*cos(p+pi);
	  yt=ht*tan(t)*sin(p+pi);

      x_spec=ht*hr*(tan(t)*cos(p+pi)+tan(ts)*cos(ps))/(ht+hr);
	  y_spec=ht*hr*(tan(t)*sin(p+pi)+tan(ts)*sin(ps))/(ht+hr);

      distance=sqrt((xt-x_spec)^2+(yt-y_spec)^2);
 
      alpha1=atan(ht/distance);
	  angolo_speculare=0.5*pi-alpha1;

         end

[rH,rV]=fresnel_Matlab(1,refr_index,rad2deg(angolo_speculare));

% r=rV;
r=(rV+rH)/2;

% r=rH;

 lam=3.d8/f;
 k=2.d0*pi/lam;


beta0rad=beta0rad*sqrt(2./log(2.));
beta2rad=beta2rad*sqrt(2./log(2.)); 

 g0=1/(R0*beta0rad);
 gs=1/(R2*beta2rad);


%%


 kx=real(k)*(sin(ts)*cos(ps)-sin(t));
 if ps==pi

     ky=0;

 else

     ky=real(k)*sin(ts)*sin(ps);

 end
 kz=k*(cos(t)+cos(ts));
%  kz=k*(cos(t)-cos(ts)); %Amir test, change + to - to change the trend 

 ax0=cos(t)^2;
 ax2=1-(sin(ts)^2)*(cos(ps)^2);
 ay2=1-(sin(ts)^2)*(sin(ps)^2);
 axy2=2*sin(ps)*cos(ps)*sin(ts)^2; 

 alphas=(cos(ts)^2)*cos(ps)^2+sin(ps)^2;
 betas=(cos(ts)^2)*sin(ps)^2+cos(ps)^2;
 gammas=-2*sin(ps)*cos(ps)*sin(ts)^2;
 alpha=(g0^2)*cos(t)^2+(gs^2)*alphas;
 beta=g0^2+(gs^2)*betas;
 gamma=(gs^2)*gammas;
 a1=2*beta;
 a2=2*alpha-(gamma^2)/a1;

 a3=beta/2+((k^2)/(4*a1)+(k^2)*(gamma^2)/(4*a2*a1^2))*(1/(R0^2)+(ay2^2)/(R2^2)+2*ay2/(R0*R2))+...
    (k^2)/(4*a2)*(axy2^2)/(4*R2^2)+((k^2)*gamma/(2*a2*a1))*(axy2/(2*R0*R2)+axy2*ay2/(2*R2^2));
 delta3=gamma/2-((k^2)/(4*a1)+(k^2)*(gamma^2)/(4*a2*a1^2))*(axy2/(R0*R2)+ay2*axy2/(R2^2))-...
        (k^2)/(4*a2)*(ax0*axy2/(R0*R2)+ax2*axy2/(R2^2))-((k^2)*gamma/(2*a2*a1))*(ax0/R0^2+ax2/(R0*R2)+ax0*ay2/(R0*R2)+ax2*ay2/(R2^2)+(axy2^2)/(4*R2^2));
 a4=alpha/2+((k^2)/(4*a1)+(k^2)*(gamma^2)/(4*a2*a1^2))*(axy2^2)/(4*R2^2)+...
    (k^2)/(4*a2)*((ax0^2)/(R0^2)+(ax2^2)/(R2^2)+2*ax0*ax2/(R0*R2))+...
    (((k^2)*gamma)/(2*a2*a1))*(ax0*axy2/(2*R0*R2)+ax2*axy2/(2*R2^2))-delta3^2/(4*a3);
 b4=-1i*kx+1i*ky*delta3/(2*a3);


 a0=abs(r)*(cos(t)+cos(ts));

%  coeff=(k^2)*(a0^2)/4*sqrt(1/a3)*sqrt(1/a4)*exp(-(kz^2)*sstd^2)*exp(-0.25*(ky^2)/a3)*exp((b4)^2/(4*a4));

 coeff=(k^2)*(a0^2)/4*sqrt(1/a3)*sqrt(1/a4)*exp(-0.25*(ky^2)/a3)*exp((b4)^2/(4*a4));


end 