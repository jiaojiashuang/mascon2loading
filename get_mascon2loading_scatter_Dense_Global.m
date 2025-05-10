function [loading_meter_u_global,loading_meter_n_global,loading_meter_e_global] = get_mascon2loading_scatter_Dense_Global(GreensFunction,grid_EWH_m,gps_lon,gps_lat,dense_factor,window_radius) %
% This code was also used for the loading calculation in several of our previous works.
% e.g., 
% Jiao et al.,2022,JGR: Solid Earth. --Basin mass changes in Finland from GRACE: Validation and explanation. [https://doi.org/10.1029/2021JB023489]
% Pan et al.,2022,EPSL. --Transient hydrology-induced elastic deformation and land subsidence in Australia constrained by contemporary geodetic measurements.[https://doi.org/10.1016/j.epsl.2022.117556]


%%
% INPUT:
% GreensFunction --> Three column (angular distance[degree],radial GreensFunction[nomalization],Tangential GreensFunction[nomalization])
% grid_EWH_m --> a global regular gridded mass field (e.g., 180*360, 360*720, 720*1440), can be GRACE mascon solutions. [mass Unit is meter of EWH]
% calculated_point_lon --> One longitude column, corresponds to latitude
% calculated_point_lat --> One latitude column, corresponds to longitude
% dense_factor --> A constant for bilinear interpolation, when the size of the input grid_EWH_m is 360*720, a dense_factor of 50 is recommended
% window_radius --> Define the range of the near field, usually, 10 is enough


%
% OUTPUT:
% loading_meter_u_global --> [mm,nn], here mm stands for the number of time nodes,
%                             nn stands for the number of calculated points
% loading_meter_n_global --> [mm,nn]
%
% loading_meter_e_global --> [mm,nn]


% ********************************
% Written by 
% JIAO Jiashuang 7/6/2021
% SGG, Wuhan University, China
% jiashuang.jiao@whu.edu.cn
% jiaojiashuang@foxmail.com
% ********************************


%%
if nargin < 6, window_radius = 10; end % 
if nargin < 5, dense_factor = 50; end % 
%%
nrow = size(grid_EWH_m,1);

if nrow==2880
    type=0.0625;
elseif nrow==1440
    type=0.125;
elseif nrow==720
    type=0.25;
elseif nrow==360
    type=0.5;
elseif nrow==180
    type=1;
else
    error('Wrong input grid size');
end

%%

% 自适应近场矩形边界选取
lon_near_west = round(gps_lon)-window_radius-type/2;
lon_near_east = round(gps_lon)+window_radius+type/2;
lat_near_north = round(gps_lat)+window_radius+type/2;
lat_near_south = round(gps_lat)-window_radius-type/2;
% 考虑到窗口半径(默认10度)可能过大
if lat_near_south < -90+type/2; lat_near_south = -90+type/2;end % 不允许超出南极
if lat_near_north > 90-type/2; lat_near_north = 90-type/2;end% 不允许超出北极
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 以下多个 'if' 为了实现自适应边界调整而设计, 此时输入的 gps_lon和gps_lat不能是矢量, 只能是一个点, 否则 && 无法判别
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lon_near_west < 0 && lon_near_east < 360 % 左侧超出
    %
    [mask_temp1]=get_mask_grid_inside(type/2,lon_near_east,lat_near_south,lat_near_north,type); % left right down up
    [loading_meter_u_near(:,1),loading_meter_n_near(:,1),loading_meter_e_near(:,1)] = get_mascon2loading_scatter_Dense(GreensFunction,mask_temp1,grid_EWH_m,gps_lon,gps_lat,dense_factor);
    [mask_temp2]=get_mask_grid_inside(360 + lon_near_west,360 - type/2,lat_near_south,lat_near_north,type); % left right down up
    [loading_meter_u_near(:,2),loading_meter_n_near(:,2),loading_meter_e_near(:,2)] = get_mascon2loading_scatter_Dense(GreensFunction,mask_temp2,grid_EWH_m,gps_lon,gps_lat,dense_factor);
    
    
    [loading_meter_u_far,loading_meter_n_far,loading_meter_e_far] = get_mascon2loading_scatter(GreensFunction,1-mask_temp1-mask_temp2,grid_EWH_m,gps_lon,gps_lat);
    
    loading_meter_u_global = sum(loading_meter_u_near,2) + loading_meter_u_far;
    loading_meter_n_global = sum(loading_meter_n_near,2) + loading_meter_n_far;
    loading_meter_e_global = sum(loading_meter_e_near,2) + loading_meter_e_far;
    %
elseif lon_near_east > 360 && lon_near_west > 0 % 右侧超出
    %
    [mask_temp1]=get_mask_grid_inside(lon_near_west,360 - type/2,lat_near_south,lat_near_north,type); % left right down up
    [loading_meter_u_near(:,1),loading_meter_n_near(:,1),loading_meter_e_near(:,1)] = get_mascon2loading_scatter_Dense(GreensFunction,mask_temp1,grid_EWH_m,gps_lon,gps_lat,dense_factor);
    [mask_temp2]=get_mask_grid_inside(type/2,lon_near_east - 360,lat_near_south,lat_near_north,type); % left right down up
    [loading_meter_u_near(:,2),loading_meter_n_near(:,2),loading_meter_e_near(:,2)] = get_mascon2loading_scatter_Dense(GreensFunction,mask_temp2,grid_EWH_m,gps_lon,gps_lat,dense_factor);
    
    
    [loading_meter_u_far,loading_meter_n_far,loading_meter_e_far] = get_mascon2loading_scatter(GreensFunction,1-mask_temp1-mask_temp2,grid_EWH_m,gps_lon,gps_lat);
    
    loading_meter_u_global = sum(loading_meter_u_near,2) + loading_meter_u_far;
    loading_meter_n_global = sum(loading_meter_n_near,2) + loading_meter_n_far;
    loading_meter_e_global = sum(loading_meter_e_near,2) + loading_meter_e_far;
    %
elseif lon_near_east > 360 && lon_near_west < 0 % 左右都超出(此时窗口宽度大于360度,错误!)
    
    error('Wrong input window_radius');
    
else % 常规一般情况
    [mask_temp]=get_mask_grid_inside(lon_near_west,lon_near_east,lat_near_south,lat_near_north,type); % left right down up
    [loading_meter_u_near,loading_meter_n_near,loading_meter_e_near] = get_mascon2loading_scatter_Dense(GreensFunction,mask_temp,grid_EWH_m,gps_lon,gps_lat,dense_factor);
    [loading_meter_u_far,loading_meter_n_far,loading_meter_e_far] = get_mascon2loading_scatter(GreensFunction,1-mask_temp,grid_EWH_m,gps_lon,gps_lat);
    
    
    loading_meter_u_global = loading_meter_u_near + loading_meter_u_far;
    loading_meter_n_global = loading_meter_n_near + loading_meter_n_far;
    loading_meter_e_global = loading_meter_e_near + loading_meter_e_far;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [mask]=get_mask_grid_inside(lon_west,lon_east,lat_south,lat_north,type) % left right down up

% 
% JIAO Jiashuang 12/04/2019
% jiaojiashuang@foxmail.com

if nargin < 5, type = 1; end % default 1 degree interval

if type==1 % input value must end with '.5'
    mask=zeros(180,360); % for 1 degree interval
    for nlat=(89.5-lat_north+1)/1:(89.5-lat_south+1)/1 % lat,
        for nlon=(lon_west+0.5)/1:(lon_east+0.5)/1 % lon
            mask(nlat,nlon)=1;
        end
    end
    
elseif type==0.5 % input value must end with '.25'
    mask=zeros(360,720); % for 0.5 degree interval
    for nlat=(89.75-lat_north+0.5)/0.5:(89.75-lat_south+0.5)/0.5 % lat,
        for nlon=(lon_west+0.25)/0.5:(lon_east+0.25)/0.5 % lon
            mask(nlat,nlon)=1;
        end
    end
    
    
elseif type==0.25 % input value must end with '.125'
    
    mask=zeros(720,1440); % for 0.25 degree interval
    % lon1 = 0.125 ;   % left
    % lon2 = 59.875;   % right
    % lat1 =  89.875;   % up
    % lat2 =  0.125;  % down
    
    for nlat=(89.875-lat_north+0.25)/0.25:(89.875-lat_south+0.25)/0.25 % lat,
        for nlon=(lon_west+0.125)/0.25:(lon_east+0.125)/0.25 % lon
            mask(nlat,nlon)=1;
        end
    end
    
elseif type==0.125 % input value must end with '.0625'
    
    mask=zeros(1440,2880); % for 0.125 degree interval
    for nlat=(89.9375-lat_north+0.125)/0.125:(89.9375-lat_south+0.125)/0.125 % lat,
        for nlon=(lon_west+0.0625)/0.125:(lon_east+0.0625)/0.125 % lon
            mask(nlat,nlon)=1;
        end
    end
    
elseif type==0.0625 % input value must end with '.0625'
    
    mask=zeros(2880,5760); % for 0.0625 degree interval
    for nlat=(89.96875-lat_north+0.0625)/0.0625:(89.96875-lat_south+0.0625)/0.0625 % lat,
        for nlon=(lon_west+0.03125)/0.0625:(lon_east+0.03125)/0.0625 % lon
            mask(nlat,nlon)=1;
        end
    end
    
else
    error('Wrong input type: 1, 0.5 or 0.25');
end
end
