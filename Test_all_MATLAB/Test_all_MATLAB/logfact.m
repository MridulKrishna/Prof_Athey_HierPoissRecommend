function lf = logfact(n)
    lf = 0;
    
    for i = 1:n
        lf = lf + log(i);
    end
end