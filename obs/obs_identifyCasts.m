function cast = obs_identifyCasts(EX)

XYT = [EX.t, EX.x, EX.y];
if sum(isnan(XYT))==3
    cast=nan;
    return
end

Unixyt = unique(XYT,'rows');
cast = nan(length(XYT),1);
C=1;
for k=1:size(Unixyt,1);
    cast(XYT==Unixyt(k))=C;
    C=C+1;
end
    
end
    
