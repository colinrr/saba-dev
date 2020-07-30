%July 14, 2019
%Find a given number of vertical or horizontal intensity profiles in square grayscale images; find mean profile from this set; build an image with this profile;
%subtract from initial image; return both the result and the inverted result
%
%INPUTS:
%inimage = grayscale image
%no_profiles = number of profiles in y (vertical) or x (horizontal)
%startingposition = pixel location for first profile
%vert_horiz = ascii 'vert' or 'horiz' to toggle vertical or horizontal; default is vertical; if
%horiz, image is rotated 90 degrees
%profiles
%polyorder = order of polynomial used in moving polynomial fit to intensity
%profile via polyfitp
%widthforpoly = width of the window for polyfitp
%run = run number


function [new_image,invert_new_image,smooth_mean_profile_y,mean_profile_y,Icx,Icy,image_to_remove,inimage] =...
    find_remove_profile(inimage,no_profiles,startingposition,vert_horiz,polyorder,widthforpoly,run)
clf
figure(run+8000),imshow(inimage) %raw image
SIZEimage=size(inimage);
inimage=imresize(inimage,[min(SIZEimage) min(SIZEimage)]);  %make the image square for calcs
%inimage=double(inimage);
image_to_remove=ones(min(SIZEimage) , min(SIZEimage));
sizeimage=size(inimage);
start=startingposition; %pixel location for first profile
yy= [1 sizeimage(2)];
interval = round((sizeimage(1))/no_profiles);
figure(run+9000),imshow(inimage),hold on
colindex=0;


if isequal(vert_horiz,'horiz')==1
    inimage=imrotate(inimage,90);
    [num2str(no_profiles),' Horizontal Intensity Profiles Used']
end
if isequal(vert_horiz,'vert')==1, inimage=imrotate(inimage,0);
    [num2str(no_profiles),' Vertical Intensity Profiles Used']
end

colindex=0;
for xx = start:interval:sizeimage
    colindex=colindex+1;
    [Icx(:,colindex),Icy(:,colindex),Ic(:,colindex),xi,yi]=improfile(inimage,[xx xx],yy);
    plot(Ic(:,colindex),Icy(:,colindex),'r'),hold on
    smooth_profile_gauss(:,colindex) = FilterGauss(Ic(:,colindex),4,3);
end

%mean_profile_y = mean(smooth_profile_gauss,2);  %Sum across columns
mean_profile_y = mean(Ic,2);
%stack_profile_y = stack_profile_y./colindex;
%stack_profile_y = mean(Ic,2);
%     %Plot profiles
%     line(Ic,Icy,'color','red'),hold on
plot(mean_profile_y,Icy,'y'),hold on

fit_profile_y = polyfitp(Icy(:,1),mean_profile_y,polyorder,widthforpoly,round(widthforpoly/2));
smooth_mean_profile_y=FilterMean(fit_profile_y,4);
smooth_mean_profile_y=smooth_mean_profile_y';
plot(fit_profile_y,Icy,'m--','LineWidth',2),hold on
plot(smooth_mean_profile_y,Icy,'g','LineWidth',4)
figure(run+5+9000),plot(smooth_mean_profile_y')  
[xgray,grayprofile]=ginput(1);%remove gray so that gradient min is at zero
smooth_mean_profile_y=smooth_mean_profile_y'-ones(length(smooth_mean_profile_y),1)*grayprofile;
image_to_remove=image_to_remove.*smooth_mean_profile_y;
figure(run+10+9000),imshow(mat2gray(image_to_remove))
%image_to_remove=image_to_remove.*-smooth_mean_profile_y;
%image_to_remove=image_to_remove.*mean_profile_y;
new_image=imsubtract(double(inimage),image_to_remove);
figure(run+10000),imshow(new_image)
new_image=mat2gray(new_image);
new_image=imgaussfilt(new_image,0.5);

figure(run+11000),imshow(new_image)
invert_new_image=(1-new_image);
figure(run+12000),imshow(invert_new_image)


