% Sep 2021 for HSAVES - Davide Comite

function  y = rotation_clockwise_mat(p, alpha, tag)

%The routine applies rotation formulas over the plane xy assuming
%a clockwise rotation (i.e., downwards)

%check if p is column vector
if isrow(p), p = p.'; end
   

switch tag
    
    case 'rad'

    
          %Aclock    = [  cos(alpha)      -sin(alpha)    0;...
          %               sin(alpha)       cos(alpha)    0;...
          %                  0                 0         1    ];
                      
          %y = Aclock*p;
          

          y(:,:,1) =  p(:,:,1).*cos(alpha) - p(:,:,2).*sin(alpha) +     0;
          y(:,:,2) =  p(:,:,1).*sin(alpha) + p(:,:,2).*cos(alpha) +     0;
          y(:,:,3) =                                                 p(:,:,3);              
    
    
    case 'deg'
  
    
          %Aclock    = [  cosd(alpha)      -sind(alpha)    0;...
          %               sind(alpha)       cosd(alpha)    0;...
          %                  0                   0         1    ];    
        
          %y = Aclock*p;    
          
                  
          y(:,:,1) =  p(:,:,1).*cosd(alpha) - p(:,:,2).*sind(alpha) +     0;
          y(:,:,2) =  p(:,:,1).*sind(alpha) + p(:,:,2).*cosd(alpha) +     0;
          y(:,:,3) =                                                 p(:,:,3);            
        
end


end

