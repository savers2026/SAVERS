function[egc]=costantedielettrica(gmoi,f)
% Questa funzione calcola la costante dielettrica del suolo in funzione
% della frequenza e della percentuale di umidità del suolo

%gmoi=0:0.01:0.4;

lambda=30/f;         %lunghezza d'onda in cm  
alpha=0.65;          
beta=1.1;            
rob=1.3;             % densità tipica della terra(g/cm^3) ???
ross=2.65;           % porosità del suolo (g/cm^3)
epsss=4.7;           % permittività del suolo solido
g1=1+(rob/ross)*(epsss^alpha-1);
f0=f/18.64; 
fghz=f*10^9;
if(f<=4)
   g2=4.9+74.1/(1+f0*i);
   g2b=0.107./(2.*pi.*fghz.*8.854.*10^-12).*(ross-rob)./(ross.*gmoi); 
   g2=(g2-g2b*i).^alpha;
else
   g2=(4.9+74.1/(1+f0*i))^alpha;
end
g3=g1+(gmoi).^beta.*(g2-1);
g4=(log(g3))/alpha;
egc=exp(g4);    %costante dielettrica del suolo (variabile immaginaria)

%plot(gmoi,real(egc),gmoi,-imag(egc))
return 
end