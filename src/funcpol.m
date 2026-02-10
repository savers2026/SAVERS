%%%%%%%%%%%%%%
function thetampl = funcpol(theta, vartheta,bias)
    thetampl = exp(-(theta).^2/vartheta)-bias; 
end
