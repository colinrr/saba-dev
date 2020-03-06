%% ================= PARAMETER INPUT ==================

project_24A
% project_24B
% project_25A
% project_25B


%% ================= DRIVER SWITCHES AND WORKFLOW ==================
%  ------ Data conversion and registration -------
% flag_tif2mat    = false; % Out of date, use asc2mat
flag_asc2mat     = false;  % IRB Ascii to .mat conversion
flag_pixelreg    = false;  % Register thermal images (correct shaking etc)

% ------ plumeTracker (mask generation), geometry and initial calcs ------
flag_mapPixels   = false;   % Generate mapping function to convert pixels to meters
flag_plumetrack  = false;   % Run Bombrun plume segmentation algorithm
flag_poly2mask   = false;   % Apply manual polygons to plume masks
flag_plumecalcs  = false;   % Basic H, v, A calcs for segmented plumes

% ------ pretty pictures --------
flag_image2dem   = false;   % Project thermal images onto the DEM
flag_riseDiag    = false;   % Plot rise diagram

% ----- Data Cube Workflow: re-grid images, spectral analysis, thermCube --------
flag_interpTherm = false;   % Interpolate thermal images/masks to regular x,z grids
flag_spectra1D   = false;   % Calc 1D (horizontal) spectra and histograms for a sliding window
flag_thermCube   = false;   % Get thermal data cube

% ----- Data Cube Workflow: velocimetry --------
flag_gradientVid = false;   % Make video of thermal gradient from thermal data cube
flag_scaledVid   = false;   % Make scaled thermal video for optic flow analysis
flag_velocimetry = true;   % Optical Flow analysis on thermal or gradient video
flag_thermCorr   = false;   % Cross-correlation analysis on thermal data cube
% flag_spectra1x1D = false; % Anisotropy spectral analysis - NOT YET IMPLEMENTED


%% ========================== DO THE THING ==========================
% if flag_tif2mat
%     irbTif2Mat(tifDir,matDir,tif_params,glob_spec)
% end

if flag_asc2mat
    irbAsc2Mat(ascDir,matDir,asc_params,glob_spec)
end

if flag_pixelreg
%     registered = plumePixelReg(matDir,reg_maskX,ref_idx,reg_idx,outputDir);
end

if flag_mapPixels
    [px2geo,geom] = mapPixels(obsLLE,refLLE,refPix,hfov,vfov,imsz);
end

if flag_plumetrack
    %                                           oDir,      ref, deb, fin, dt
    [content,ref] = mainTrackPlume(procDir,outputDir, ref, deb,fin,dN,fixpix,[plume_params,png,gif,video],delT); 
%     content
end

if flag_poly2mask
    [T,Ref,update_time] = maskPolyClip(paramf,polyFile,true);
end

if flag_plumecalcs
    [dat,geom] = plumeTrackCalcs(outputDir,obsLLE,refLLE,refPix,px2geo,geom,plotCalcs);
end

if flag_image2dem
    dem_frames = 2; %[1:size(dat,1)];

    lp=image2dem(dat,px2m,geom,demfile,dem_roi,matDir,dem_frames);
end

if flag_interpTherm
    interpThermal(matDir,interpDir,{mat_params,geomf},geomf,interpIdx,[],[],polyFile)
end

if flag_spectra1D
    [Tspectral,sParam,Kolm,oname] = plumeSpec1D(matDir,interpDir,{paramf},geomf,sParam,ofile);
%     load(fullfile(interpDir,oname))
%     plotPlumeHist(Tspectral,sParam,interpDir,histidx,histflags)
end  

if flag_thermCube
    [D,thermFile] = getThermCube(interpDir,mat_params,geomf,thermIdx,fixpix,ROI,cubeDir);
%     plotThermStats
end

if flag_gradientVid
    [vidParams,vidParFile] = thermCube2gradient(thermFile,vidPar,ovid,vidFlags);
end
if flag_scaledVid
    [vidParams,vidParFile] = thermCube2video(thermFile,vidPar,ovid);
end
if flag_velocimetry
    [V] = thermOpticFlow(velVid, vidParFile, FBparams, opticPlotFrames);
end

if flag_thermCorr
    load(thermFile)
    C = thermalCorrelation(D,fnames,taperN,Nwin,Nover,filt,defHighPass);    
end

if flag_riseDiag
%     D = plotRiseDiagram(matDir,fullfile(outputDir,'geometry.mat'),ventPix,profileType,Ridx);
    RiseandSpec
end

%%  Run a few vid scaling tests


% load(thermFile)
% 
% % idxt = [150:50:500];
% idxt = [850:-5:5];
% ROI  = [400 641 1 250];
% 
% hi_thresh = 98; % Percentile
% lo_thresh = 0.15; % Threshold (fraction) for percent difference in counts b/w non-masked and masked histograms
% 
% Tscale0 = [230 424.74]; % Initial scaling to normalize images universally
% 
% Scaling for image adjusting (equivalent to low/high in values in
% imadjust)

% Tminmax = [%240 400;
% %            0.15 0.98;
% %            0.15 0.99;
%            0.10 0.99;
% %            0.10 310;
%             ];
%        
% % lohi_in = [0 1];
% %            0.15 0.7];
% gamma   = [0.7 0.85 1.0];
% 
% hi_temps = zeros(size(idxt));
% lo_temps = hi_temps;
% hi_vals  = hi_temps;
% lo_vals  = lo_temps;
% 
% ny = 1; %size(lohi_in,1);
% nx = ceil(size(Tminmax,1)/ny);
% ng = numel(gamma);
% nh = 0;
% 
% for ii = 1:numel(idxt)
%     Frame = flipud(D.T(:,:,idxt(ii)));
%     Frame = Frame(ROI(1):ROI(2),ROI(3):ROI(4));
%     
%     Mask  = flipud(D.mask(:,:,idxt(ii)));
%     Mask  = Mask(ROI(1):ROI(2),ROI(3):ROI(4));
%     
%     FrMask = Frame.*Mask;% idxt = [50:50:800]; % 144 300 600];
% idxt = fliplr(idxt);
% 
% ROI  = [400 641 1 250];
% 
% hi_thresh = 'prctile'; % Or 'temperature'
% 
% Tminmax = [240 400;
%            240 350];
% lohi_in = [0 1];
% %            0.15 0.7];
% gamma   = [0.6 0.8 1.0];
% 
% nx = size(Tminmax,1);
% ny = size(lohi_in,1);
% ng = numel(gamma);
% 
% for ii = 1:numel(idxt)
%     Frame = flipud(D.T(:,:,idxt(ii)));
%     Frame = Frame(ROI(1):ROI(2),ROI(3):ROI(4));
%     
%     figure('position',[50 50 1800 1000])
%     for ti = 1:nx
%         Tmn = Tminmax(ti,1);
%         Tmx = Tminmax(ti,2);
% 
%         Fr1 = Frame;
%         % Remove min
%         Fr1 = Fr1 - Tmn;
%         Fr1(Fr1<0) = 0;    
% 
%         % Remove max and scale to 1
%         Fr1 = Fr1./(Tmx-Tmn);
%         Fr1(Fr1>1) = 1;
% 
%         [FrN,FrEdge]=histcounts(Fr1);
%         FrMid = FrEdge(1:end-1)+diff(FrEdge/2);    
% 
%         for li = 1:ny
%             
%             lh_in = lohi_in(li,:);
% 
%             % n = nx*(li-1)+ti
%             nh = nx*(2*li - 1) + ti;
% 
%             histax = tightSubplot(2*ny,nx,nh);
%             plot(histax,FrMid,FrN,'Color',[0.2 0.2 0.2],'LineWidth',3)
%             xlim(histax,[0 1])
%             grid(histax,'on')
%             hold(histax,'on')
%             set(histax,'YScale','log')
%             yl = ylim(histax);
%             plot([lh_in; lh_in], [yl' yl'],'--','Color',[0.6 0.6 0.6],'LineWidth',2 )
% 
%             for gi = 1:ng
%     %             Fr = Frame;
%                 Fr = Fr1;
% 
%                 gam = gamma(gi);
% 
%                 nim = nx*ng*2*(li - 1) + (ti-1)*ng + gi;
% 
%                 % i=1,j=1,g=3, n = 3
%                 % i=1,j=2,g=2  n = 5; 
%                 % i=2,j=2,g=2  n = 17
% 
% 
% 
% 
% 
% 
%                 % Apply imadjust
%                 Fr = imadjust(Fr,lh_in,[],gam);
% 
%                 tightSubplot(2*ny,nx*ng,nim);
%                 imagesc(Fr)
%                 colormap(gray)
%                 box off
%                 text(0.6,0.9,sprintf('g=%.2f',gam),'Units','normalized','Color','w','Fontsize',12)
%                 set(gca,'XTick',[],'YTick',[])
% 
%                 histogram(histax,Fr)
% 
%             end
% 
% 
%                 text(histax,0.3,0.9,sprintf('T=[%.0f %.0f], lohi=[%.2f %.2f]',Tmn,Tmx,lh_in(1),lh_in(2)),'Units','normalized','Fontsize',12)
%                 legend(histax,'Raw','lo in','hi in',sprintf('g %.2f',gamma(1)),sprintf('g %.2f',gamma(2)),sprintf('g %.2f',gamma(3)))
%         end
%     end
% 
% end

%     FrLin  = FrMask(FrMask>0);
%     
%     [NFr0,edges]=histcounts(Frame);
%     Frbins = edges(1:end-1) + diff(edges)/2;
%     
%     Nm0 = histcounts(FrLin,edges);
%     prcdiff = (NFr0-Nm0)./NFr0;
%     
%     % Lower bound temperature where full histogram and mask histogram
%     % differ by more than <lo_thresh>%.
% %     Tprcmin = Frbins(find(prcdiff>lo_thresh,1,'last')); 
% %     Tprcmax = prctile(FrLin,hi_thresh);
%     
% %     figure('position',[50 50 1800 1000])
%     
%     % Decision for temperature limits
%     for ti = 1:size(Tminmax,1)
% 
%         if and(Tminmax(ti,1)>=0,Tminmax(ti,1)<=1)
%             Tmn = Frbins(find(prcdiff>Tminmax(ti,1),1,'last')); 
% %             Tmn = Tprcmin;
%         else
%             Tmn = Tminmax(ti,1);
%         end
% 
%         if and(Tminmax(ti,2)>=0,Tminmax(ti,2)<=1)
%             Tmx = prctile(FrLin,Tminmax(ti,2)*100);
% %             Tprcmax = Tmx;
%         else
%             Tmx = Tminmax(ti,2);
%         end
%     
%         hi_temps(ii) = Tmx;
%         lo_temps(ii) = Tmn;
%         % ----- Initial Conversion to [0 1] scale ------
%         Tmn = (Tmn-Tscale0(1))/diff(Tscale0);
%         Tmx = (Tmx-Tscale0(1))/diff(Tscale0);
%         
%         Fr1 = Frame;
%         % Remove min
%         Fr1 = Fr1 - Tscale0(1);
%         Fr1(Fr1<0) = 0; 
% 
%         % Remove max and scale to 1
%         Fr1 = Fr1./diff(Tscale0);
%         Fr1(Fr1>1) = 1;
%         % ----------------------------------------------
% 
%         hi_vals(ii) = Tmx;
%         lo_vals(ii) = Tmn;
%         
%         % Normalized histograms
%         [FrN,FrEdge]=histcounts(Fr1);
%         FrMid = FrEdge(1:end-1)+diff(FrEdge/2);
% 
% %         for li = 1:ny
% 
% %             lh_in = lohi_in(li,:);
% 
%         % n = nx*(li-1)+ti
% %             nh = nx*(2*li - 1) + ti;
%         nh =  nx*ceil(ti/nx)+ti; %(nx+1)*(ti);
%         % ti = 1,    nh = 3
%         % ti = 2     nh = 4
%         % ti = 3     nh = 7
%         % ti = 4     nh = 8
% 
%         
%         
% %         histax = tightSubplot(2*ny,nx,nh);
% %         plot(histax,FrMid,FrN,'Color',[0.2 0.2 0.2],'LineWidth',3)
% %         xlim(histax,[0 1])
% %         grid(histax,'on')
% %         hold(histax,'on')
% %         set(histax,'YScale','log')
% %         yl = ylim(histax);
% %         plot([Tmn Tmx; Tmn Tmx], [yl' yl'],'--','Color',[0.6 0.6 0.6],'LineWidth',2 )
% 
%         for gi = 1:ng
%     %             Fr = Frame;
%             Fr = Fr1;
% 
%             gam = gamma(gi);
% 
% %             nim = nx*ng*2*(li - 1) + (ti-1)*ng + gi;
% %             nim = ng*(ti-1)+gi;
% %             nim = nx*ng*(ceil(ti/nx)-1) + gi;
%             nim = ng*(nx*ceil(ti/nx) + ti - nx -1) + gi;
%             % i=1,j=1,g=3, n = 3
%             % i=1,j=2,g=2  n = 5; 
%             % i=2,j=2,g=2  n = 17
% 
%             % Apply imadjust
%             Fr = imadjust(Fr,[Tmn Tmx],[],gam);
% 
% %             tightSubplot(2*ny,nx*ng,nim);
% %             imagesc(Fr)
% %             colormap(gray)
% %             caxis([0 1])
% %             box off
% %             text(0.6,0.9,sprintf('g=%.2f',gam),'Units','normalized','Color','w','Fontsize',12)
% %             set(gca,'XTick',[],'YTick',[])
% % 
% %             histogram(histax,Fr,FrEdge)
% 
%         end
% 
% %             text(histax,0.3,0.9,sprintf('T=[%.2f %.2f], lohi=[%.2f %.2f]',Tminmax(ti,1),Tminmax(ti,2),Tmn,Tmx),'Units','normalized','Fontsize',12)
% % %             text(histax,0.3,0.8,sprintf('Tprc=[%.0f %.0f]',Tprcmin,Tprcmax),'Units','normalized','Fontsize',12)            
% %             legend(histax,'Raw','lo in','hi in',sprintf('g %.2f',gamma(1)),sprintf('g %.2f',gamma(2)),sprintf('g %.2f',gamma(3)))
% %         end
%     end
% end
