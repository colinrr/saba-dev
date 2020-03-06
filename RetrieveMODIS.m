function [SS,SSraw] = RetrieveMODIS(hdfile,sds_fields,roi)
% SS = RetrieveRawMODIS(hdfile,sds_fields,roi)
% IN:   hdfile     = full path to MODIS hdf file
%       sds_fields = cell array, with each entry being a two element cell
%               containing {output field name, field name from HDF file}
%               eg. {'Temperature', 'Retrieved_Temperature_Profile'}
%       roi        = OPTIONAL 2 elementcell array containing lat/lon limits
%               of interest. Eg. { [lat1 lat2], [lon1 lon2]}
%           ---> Doesn't seem to work with hdfread right now...
% OUT:  SS = Data structure, where each entry is a structure containing data and
%            attributes for the named field
%
% C Rowell, September 2018

    if nargin<3
        roi = [];
    end
%     
%     if isempty(roi)
%         readfun = @(x) hdfread(x);
%     else
%         readfun = @(x) hdfread(x,'Box',roi);
%     end
    % meta_eos = hdfinfo(hdfile,'eos');
    
    fprintf('Retrieving data from file:\n\t%s\n',hdfile)
    meta     = hdfinfo(hdfile);
    
    % Lat/lon/pressure are kinda their own thing...
    geoi = strcmp({meta.Vgroup.Vgroup.Name},'Geolocation Fields');
    geo = meta.Vgroup.Vgroup(geoi).SDS;

    for a=1:length(geo)
        name   = geo(a).Name;
        s.data = hdfread(geo(a));
        attr   = geo(a).Attributes;
        for b = 1:length(attr)
            s.(genvarname(attr(b).Name)) = attr(b).Value;
        end
        SS.(name) = s;
        clear s
    end

    % Retrieve SDS fields along with meta data
    sdsi = strcmp({meta.Vgroup.Vgroup.Name},'Data Fields');
    pi   = strcmp({meta.Vgroup.Vgroup(sdsi).Vdata.Name},'Pressure_Level');
    sds = meta.Vgroup.Vgroup(sdsi).SDS;
    pinf = meta.Vgroup.Vgroup(sdsi).Vdata(pi);

    % Pressure
    SS.Pressure.data = hdfread(pinf);
    SS.Pressure.data = SS.Pressure.data{:};
    SS.Pressure.units = pinf.DataAttributes(2).Value;

    % Get rough time for now
%     [~,fn,~] = fileparts(hdfile);
%     fn = split(fn,'.');
    
    % The rest
    for jj = 1:length(sds_fields)
        name = sds_fields{jj}{1};
        field = sds_fields{jj}{2};
        s.data = hdfread(hdfile,field);
        di = strcmp({sds.Name},field);
        if ~any(di); error('Field name not found in data file.');end
        attr = sds(di).Attributes;
        for b = 1:length(attr)
            s.(genvarname(attr(b).Name)) = attr(b).Value;
        end
        SS.(name) = s;
        clear s
    end
    
    SSraw = SS;
    SS = CorrectRawMODIS(SS,roi);
end

function SS = CorrectRawMODIS(SS,roi)
% Applies correction to data retrieved from an hdf file, using the built in
% scale factors, offsets, and valid ranges from meta data. Also apply my
% own roi since the hdfread one isn't working atm.

sfields = fieldnames(SS);
corr_fields = {'valid_range','x_FillValue','scale_factor','add_offset'};

%% Apply data correction
disp('Applying data correction to retrieved data')
for ii = 1:length(sfields)
    ss = SS.(sfields{ii});
    field_check = isfield(ss,corr_fields);
    
%     if strcmp(sfields{ii},'Surface_Pressure')
%         figure
%         imagesc(SS.(sfields{ii}).data)
%     end

    if all(field_check)
        % Check for non-fill values in Non-valid ranges
        valid_check = ~or(ss.data(:)==ss.x_FillValue,and(ss.data(:)>ss.valid_range(1),ss.data(:)<ss.valid_range(2)));
        if any(valid_check)
            fprintf('Warning! Invalid values found in:\t%s\n',sfields{ii})
            ss.invalid = valid_check;
        end
        
        % Convert ot double and clear out fill values
        ss.data = double(ss.data);
        ss.valid_range = double(ss.valid_range);
        nn = ss.data==ss.x_FillValue;
        ss.data(nn) = NaN;
        ss.x_FillValue = NaN;
        
        % OR just clear out all non-valid values?
%         nn = and(ss.data(:)>=ss.valid_range(1),ss.data(:)<=ss.valid_range(2));
%         ss.data(~nn) = NaN;
        
        % Apply scale and offset
        ss.data = ss.scale_factor*(double(ss.data) - ss.add_offset);
        ss.valid_range = ss.scale_factor*(double(ss.valid_range) - ss.add_offset);
        
        ss.Correction_Applied = true;
        SS.(sfields{ii}) = ss;
    else
        fprintf('Skipping correction for: %s\n\t Missing correction fields:\n',sfields{ii})
        corr_fields(~field_check)
    end
    
    % Check elevations versus surface elevation (or pressures)
%     if isfield(SS,'Surface_Pressure')
end
    
%% Cut values below surface
% if isfield(SS,'Surface_Elevation')
%     disp('')
% elseif isfield(SS,'Surface_Pressure')
%     disp('')
% end
% for ii = 1:length(sfields)
% end
%% Apply ROI here
% Straight up boxing it for now
if ~isempty(roi)
    lats = roi{1};
    lons = roi{2};
    
    lati = and(SS.Latitude.data>=min(lats),SS.Latitude.data<=max(lats));
    loni = and(SS.Longitude.data>=min(lons),SS.Longitude.data<=max(lons));
    
    if ~any(and(lati,loni))
        figure
        plot(SS.Longitude.data(:),SS.Latitude.data(:),'ob')
        title('MODIS Pixels and ROI')
        hold on
        plot([lons lons [lons'; lons']],[[lats';lats'] lats lats],'k','LineWidth',2)
    end
    
    for ii = 1:length(sfields)
        % Quick size check
        ddims = size(SS.(sfields{ii}).data);
        gdims = size(lati);
        [A,B] = ismember(gdims,ddims);
        if all(A)
            xx = SS.(sfields{ii}).data;
            if ndims(xx)==2
                xx = xx(and(lati,loni));
            elseif ndims(xx)==3
                xx = xx(:,and(lati,loni));
            end
            SS.(sfields{ii}).data = xx;
        end
    end
end

% Visualize lat/lons
% figure
% scatter(SS.Longitude.data(:),SS.Latitude.data(:),10,[0.5 0.5 0.5],'.')
end
