function y = sinc_dc(argsinc)

%Created by Davide Comite on Nov 2021

                    i0sinc        = argsinc == 0;                                                                                    
                    sincF         = sin(pi*argsinc)./(pi*argsinc);                                                      
                    sincF(i0sinc) = 1; 

                    y = sincF;

