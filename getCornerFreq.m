function K0 = getCornerFreq(K,Pxx)
    Ki = K~=0;
    K = K(Ki);
    Pxx = Pxx(Ki);

    Kint = linspace(min(log10(K)),max(log10(K)),numel(K));
    Pint = interp1(log10(K),log10(Pxx),Kint);

    % Find change points
    ilt=findchangepts(Pint,'MaxNumChanges',1,'Statistic','linear');
    K0 = 10.^Kint(ilt);

end