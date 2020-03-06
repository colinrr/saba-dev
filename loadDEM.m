function [A,X,Y,ginfo] = loadDEM(demfile,roi)

if nargin<2
    roi=[];
end

[dem, demR] = geotiffread(demfile);
ginfo = geotiffinfo(demfile);

if isempty(roi)
    A = double(dem); aR = demR;
else
    [A,aR] = CropGeotiff(dem,demR,roi);
end

[X,Y] = pixcenters(aR,size(A),'makegrid');

end