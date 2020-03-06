function [T,Ref,update_time] = maskPolyClip(PTfile,polyFile,saveflag)
% PTfile = path to plumeTrack output file
% polyFile = path to Polygon file (from preProcThermal)
% saveflag = true save PTfile_poly.mat

%% 
disp('===== Clipping masks with manual polygons =====')

load(PTfile)
load(polyFile)

% Assuming all images are the same size (bloody well should be)
imsz = size(T.Mask{1});

Idx = cellfun(@(x) str2num(x), T.Properties.RowNames);

for ii=1:length(Polys)
%     Polys(ii).mask = poly2mask(Polys(ii).P(:,1), Polys(ii).P(:,2), imsz(1), imsz(2) );
%     if ~Polys(ii).Type
%         Polys(ii).mask = ~Polys(ii).mask;
%     end

    for pidx = Polys(ii).Idx;
        if ismember(pidx,Idx)
            T.Mask{num2str(pidx)} = T.Mask{num2str(pidx)}.*Polys(ii).mask;
        end
    end
end



if saveflag
    [fdir,fname,fext] = fileparts(PTfile);
    fout = fullfile(fdir,[fname '_poly' fext])
end

end