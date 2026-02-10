clear all, close all
% cd(1)=complex(4.66,-.29); %SMC=10%
% cd(2)=complex(11.3,-1.279); %SMC=25%
% cd(3)=complex(23.1,-3.04);    %SMC=40%
% f=13.57;  %GHz
f=1.5;  %GHz
lambda=30/f; %cm
Col=['b', 'r', 'g','c','k','m'];
SM=[5 20 40];
sz=[0.5 1.5 3];  %cm
theta_sp=0:90;
waven=2*pi/lambda;  %cm-1
beckspiz=2*waven*cosd(theta_sp);

for i=1:length(SM)
cd(i)=costantedielettrica(SM(i)/100,f);
refr_index=sqrt(cd(i));
[rH(:,i),rV(:,i)]=fresnel_Matlab(1,refr_index,theta_sp);
LR(:,i)=(abs((rH(:,i)+rV(:,i))./2)).^2;
RR(:,i)=(abs((rH(:,i)-rV(:,i))./2)).^2;
figure(1)
plot(theta_sp,abs(rH(:,i)).^2,'b','Linewidth',2)
hold on, plot(theta_sp,abs(rV(:,i)).^2,'r','Linewidth',2)
text(20,abs(rH(1,i))^2,num2str(SM(i)))

figure(22)
hold on, plot(theta_sp,LR(:,i),'r','Linewidth',2)
hold on, plot(theta_sp,RR(:,i),'b','Linewidth',2)
text(20,abs(LR(1,i)),num2str(SM(i)))

figure(2)
hold on, plot(theta_sp,db(LR(:,i),'power'),'r','Linewidth',2)
hold on, plot(theta_sp,db(RR(:,i),'power'),'b','Linewidth',2)
Leg{i}=[num2str(SM(i)) '%'];
end
figure(1)
box on
xlabel('Incidence Angle')
ylabel('Reflectivity')
legend('H','V','Location','northwest')

figure(2),box on
ylim([-30 0]);
xlabel('Incidence Angle')
ylabel('Reflectivity [dB]')
legend('LR','RR','Location','southwest')
figure(22),box on
xlabel('Incidence Angle')
ylabel('Reflectivity')
legend('LR','RR')
for isz=1:3
    factor=(beckspiz*sz(isz)).^2;
    factor=exp(-factor);
    Legsz{isz}=[num2str(sz(isz)) 'cm'];
for ism=1:3
Vrough(:,ism,isz)=abs(rV(:,ism)).*factor(:);
Hrough(:,ism,isz)=abs(rH(:,ism)).*factor(:);
LRrough(:,ism,isz)=LR(:,ism).*factor(:);
RRrough(:,ism,isz)=RR(:,ism).*factor(:);
end
end
figure

for isz=1:3

for ism=2:2
% hold on, h(ism)=plot(theta_sp,db(LRrough(:,ism,isz),'power'), Col(ism),'Linewidth',2);
% hold on, plot(theta_sp,db(RRrough(:,ism,isz),'power'),Col(ism),'Linewidth',2)
hold on, h(isz)=plot(theta_sp,db(Vrough(:,ism,isz),'power'), Col(isz),'Linewidth',2);
hold on, plot(theta_sp,db(Hrough(:,ism,isz),'power'),[Col(isz) '--'],'Linewidth',2)

end
% title([' sz = ' num2str(sz(isz)) ' cm'] )
title([' SM = ' num2str(SM(ism)) ' %'] )
hold on, %ylim([-40 0]);
legend(h, Legsz,'Location','southeast')
xlabel('Incidence Angle')
ylabel('Reflectivity [dB]')
end
% figure
% for isz=1:2:3
% ipl=1;
% 
% for ith=11:15:45
% if isz==1
% 
% hold on, plot(SM,db(Hrough(ith,:,isz),'power'),Col(ipl),'Linewidth',2)
% else
%     hold on, plot(SM,db(Hrough(ith,:,isz),'power'),[Col(ipl) '--'],'Linewidth',2)
% end 
% ipl=ipl+1;
% end
% end
% xlabel('Soil Moisture [%]')
% ylabel('Reflectivity H [dB]')
% 
% figure
% for isz=1:2:3
% ipl=1;
% 
% for ith=11:15:45
% if isz==1
% 
% hold on, plot(SM,db(Vrough(ith,:,isz),'power'),Col(ipl),'Linewidth',2)
% else
%     hold on, plot(SM,db(Vrough(ith,:,isz),'power'),[Col(ipl) '--'],'Linewidth',2)
% end 
% ipl=ipl+1;
% end
% end
% xlabel('Soil Moisture [%]')
% ylabel('Reflectivity V [dB]')
% figure
% for isz=1:2:3
% ipl=1;
% for ism=1:2:3
% if isz==1
%     hold on, plot(theta_sp,db(Hrough(:,ism,isz),'power'),Col(ipl),'Linewidth',2,'DisplayName','')
% else
%     hold on, plot(theta_sp,db(Hrough(:,ism,isz),'power'),[Col(ipl) '--'],'Linewidth',2)
% end 
% ipl=ipl+1;
% end
% end
% xlabel('Incidence Angle')
% ylabel('Reflectivity H [dB]')
% 
% figure
% for isz=1:2:3
% ipl=1;
% for ism=1:2:3
% 
% if isz==1
% 
% hold on, plot(theta_sp,db(Vrough(:,ism,isz),'power'),Col(ipl),'Linewidth',2)
% else
%     hold on, plot(theta_sp,db(Vrough(:,ism,isz),'power'),[Col(ipl) '--'],'Linewidth',2)
% end 
% ipl=ipl+1;
% end
% end
% xlabel('Incidence Angle')
% ylabel('Reflectivity V [dB]')
% 
% figure
% k=0;
% isz=1;   %Il rapporto RR/LR è indipendente da sz
% for ith=10:5:30
%     k=k+1;
% hold on, plot(SM,db(RRrough(ith,:,isz)./LRrough(ith,:,isz),'power'),'Linewidth',2)
% Leginc{k}=[num2str(ith) '°'];
% 
% end
% hold on
% xlabel('SM %')
% ylabel('RR/LR [dB]')
% ylim([-50 0]);
% grid on
% box on
% legend(Leginc)
% 
% figure(33)
% for isz=1:3
% for ism=1:3
% hold on, hold on, plot(db(LRrough(5:40,isz,ism),'power'),db(RRrough(5:40,isz,ism),'power'),['*' Col(ism)])
% end
% end
% box on
% xlabel('LR','Fontsize',15)
% ylabel('RR','Fontsize',15)
% % legend(Leg,'Location','northeast','Fontsize',15)
% figure(44)
% for isz=1:3
% for ism=1:3
% hold on, hold on, plot(db(LRrough(5:40,isz,ism),'power'),db(RRrough(5:40,isz,ism),'power'),['*' Col(isz)])
% end
% end
% box on
% xlabel('LR','Fontsize',15)
% ylabel('RR','Fontsize',15)
% 
% LRvec=[];
% RRvec=[];
% for isz=1:3
% for ism=1:3
% LRvec=[LRvec; LRrough(1:60,isz,ism)];
% RRvec=[RRvec; RRrough(1:60,isz,ism)];
% end
% end
% % polfitlin = fit(LRvec, RRvec, 'poly1');
% % figure(34)
% % plot(LRvec,RRvec,'*','Linewidth',2)
% % hold on, plot([0 0.4],[polfitlin.p2 polfitlin.p1*(0.4)+polfitlin.p2],'Linewidth',2,'Color','red')
% % text(0.2,0.05,['RR = ' num2str(polfitlin.p1) '*LR +' num2str(polfitlin.p2)],'Color','red','FontSize',10)
% % box on
% % ylim([0 0.1]);
% % xlabel('LR','Fontsize',15)
% % ylabel('RR','Fontsize',15)
% % % legend(Leg,'Location','northeast','Fontsize',15)
% % ft2=fittype('k*x');
% % pp2=fit(LRvec, RRvec, ft2);
% % hold on, plot([0 0.4],[0 pp2.k*(0.4)],'Linewidth',2,'Color','green')
% % text(0.07,0.07,['RR = ' num2str(pp2.k) '*LR'],'Color','green','FontSize',10)
% % figure
% % plot(db(LRvec,'power'),db(RRvec,'power'),'*','Linewidth',2)
% % hold on, plot([-35 0], [-35 0],'k'); 
% % ylim([-35 -5]); xlim([-35 -5])
% % text(-30,-10,['RR = ' num2str(db(pp2.k,'power')) '+LR'],'Color','red','FontSize',10)
% 
% SM=[5 10 20 30 40 50 60 70 80];
% LR=[];
% RR=[];
% for i=1:length(SM)
% cd(i)=costantedielettrica(SM(i)/100,f);
% refr_index=sqrt(cd(i));
% [rve(i),rvh(i)]=fresnel_Matlab(1,refr_index,20);
% LR(i)=(abs((rve(i)+rvh(i))./2)).^2;
% RR(i)=(abs((rve(i)-rvh(i))./2)).^2;
% end
% 
% figure(10), plot(SM,LR)
%  xlabel('SM %','Fontsize',15)
% ylabel('LR','Fontsize',15)
