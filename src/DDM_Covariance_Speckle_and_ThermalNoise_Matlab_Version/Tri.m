function y = Tri(x, tc)
    if abs(x) <= tc
        y = 1.0 - abs(x/tc);
    else
        y=0.0;
    end
end
