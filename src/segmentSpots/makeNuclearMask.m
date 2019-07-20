function nuclearMask = makeNuclearMask(ellipsesFrame, dim)

radScale = 1.3;

nuclearMask = false(dim(1), dim (2)); 

for e = 1:length(ellipsesFrame)
    
    ceny = ellipsesFrame(e, 1);
    cenx = ellipsesFrame(e, 2);
    rad = ellipsesFrame(e,3)*radScale;
    if cenx - rad >= 0 & cenx + rad <= dim(1) & ceny - rad >= 0 & ceny + rad <= dim(2)
    for x = ceil(cenx - rad): ceil(cenx + rad)
        for y = ceil(ceny - rad): ceil(ceny + rad)
            if norm([x, y] - [cenx, ceny]) < rad 
                nuclearMask(x, y) = true;
            end
        end
    end
    
end




end