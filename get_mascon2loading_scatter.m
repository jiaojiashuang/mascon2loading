function [loading_meter_u,loading_meter_n,loading_meter_e] = get_mascon2loading_scatter(GreensFunction,mask_mascon_area,grace_mascon,calculated_point_lon,calculated_point_lat)
% This is a sub-function of 'get_mascon2loading_scatter_Dense_Global'
%
%%
% Code for calculating 3-D loading displacement caused by surface mass changes in idealized spherical reference frame;
% Input Load mass['meter' in the unit of EWH] must be global mass grids,such as GRACE MASCON solutions
%
% Here, calculated displacement points[GPS points] are discrete points.
%
% INPUT:
% GreensFunction --> Three column (angular distance[degree],radial GreensFunction[nomalization],Tangential GreensFunction[nomalization])
% mask_mascon_area --> Usually, we don't need to use all global mass grids
% calculated_point_lon --> One longitude column, corresponds to latitude
% calculated_point_lat --> One latitude column, corresponds to longitude
%
%
% OUTPUT:
% loading_meter_u --> [mm,nn], here mm stands for the number of time nodes,
%                             nn stands for the number of calculated points
% loading_meter_n --> [mm,nn]
%
% loading_meter_e --> [mm,nn]
%
%
% ********************************
% JIAO Jiashuang 28/1/2021
% jiaojiashuang@foxmail.com
%
% Last modified 13/4/2021
% JIAO Jiashuang
% SGG, Wuhan University, China
% jiashuang.jiao@whu.edu.cn
% jiaojiashuang@foxmail.com
% ********************************
%%
% % a  =  6378136.3; % semi-major axis of ellipsoid [m]
a = 6371000;% [m], refers to Farrell1972_Deformation_Earth_Surface_Loads
% e2 = 0.0; % Eccentricity
ro = 1000; % water density: kg/m^3
%
%
%%
mask_mascon_area(mask_mascon_area >0 )=1; % 确保 mask 值只含有 0,1

%%
[nrow,ncol] = size(mask_mascon_area);

if (nrow==2880 && ncol==5760)
    lon=0.03125:0.0625:359.96875;
    lat=89.96875:-0.0625:-89.96875;
    type=0.0625;
elseif (nrow==1440 && ncol==2880)
    lon=0.0625:0.125:359.9375;
    lat=89.9375:-0.125:-89.9375;
    type=0.125;
elseif (nrow==720 && ncol==1440) % for altimetry data 0.25
    lon = 0.125:0.25:359.875;
    lat = 89.875:-0.25:-89.875;
    type=0.25;
elseif (nrow==360 && ncol==720)
    lon=0.25:0.5:359.75;
    lat=89.75:-0.5:-89.75;
    type=0.5;
elseif (nrow==180 && ncol==360)
    lon=0.5:359.5;
    lat=89.5:-1:-89.5;
    type=1;
else
    error('Wrong input grid size');
end

%%
[lon_global,lat_global]=meshgrid(lon,lat);
%%
area_grid = cos(lat_global*pi/180).*(type*pi/180)*(type*pi/180)* a* a/10^6;% unit is km^2;
%%
[mass_point_lat,mass_point_lon,zero_index] = get_loading_parameter_lat_lon(mask_mascon_area,nrow,ncol,lon_global,lat_global);
%%
num_gps = length(calculated_point_lon);

jjj = size(grace_mascon,3);

% pre-set for speed
loading_meter_u = zeros(jjj,num_gps);
loading_meter_n = zeros(jjj,num_gps);
loading_meter_e = zeros(jjj,num_gps);

for iii = 1:num_gps
    
    tic
    [G_radial,G_tangent_n,G_tangent_e] = get_loading_parameter_grn(calculated_point_lat(iii),calculated_point_lon(iii), mass_point_lat, mass_point_lon, a, GreensFunction);
    
    for www = 1:jjj%
        
        tic
        
        grace_mascon_temp = grace_mascon(:,:,www);
        [mass_kg] = get_loading_parameter_mass(grace_mascon_temp,area_grid,nrow,ncol,zero_index,ro);
        [disp_u,disp_n,disp_e] = get_loading_disp(G_radial,G_tangent_n,G_tangent_e,mass_kg);
        
        loading_meter_u(www,iii) = disp_u;
        loading_meter_n(www,iii) = disp_n;
        loading_meter_e(www,iii) = disp_e;
        
        toc
    end
    
end


end

function [disp_u,disp_n,disp_e] = get_loading_disp(G_radial,G_tangent_n,G_tangent_e,mass_kg)


% sum of individual segments on vertical deformations
%D_radial = G_radial.*grace_xyz_area.*grace_xyz_mass.*1e6;
D_radial_vector = G_radial.*mass_kg; %
%         nany = isnan(D_radial_vector);
%         D_radial_vector(nany) = 0;
disp_u = sum(sum(D_radial_vector));


% not sure below
N_vector = G_tangent_n.*mass_kg;% positive is north
%         nany = isnan(N_vector);
%         N_vector(nany) = 0;

E_vector = G_tangent_e.*mass_kg;% positive is east
%         nany = isnan(E_vector);
%         E_vector(nany) = 0;

disp_n = sum(sum(N_vector));
disp_e = sum(sum(E_vector));

end

function [mass_point_lat,mass_point_lon,zero_index] = get_loading_parameter_lat_lon(mask_mascon_area,nrow,ncol,global_ewh_grid_lon,global_ewh_grid_lat)

%[nrow,ncol] = size(mask_mascon_area);
%global_mass_kg_grid = msk_mar.*1.*1e6.*rho_water/1000;% unit: mmWE -> kg; Here the area of each MAR grid is 1 km^2
msk_mar_vector = reshape(mask_mascon_area,[nrow*ncol,1]);
%   nan_yes = isnan(msk_mar_vector);
%msk_mar_vector(nan_yes) = [];

[zero_index,~] = find(msk_mar_vector(:,1) == 0);
% msk_mar_vector(m,:)=[];

mass_point_lon = reshape(global_ewh_grid_lon,[nrow*ncol,1]);
mass_point_lon(zero_index,:)=[];
mass_point_lat = reshape(global_ewh_grid_lat,[nrow*ncol,1]);
mass_point_lat(zero_index,:)=[];

end

function [mass_kg] = get_loading_parameter_mass(global_ewh_grid,area_grid,nrow,ncol,zero_index,rho_water)

%[nrow,ncol] = size(global_ewh_grid);
%area_grid = get_mask2area_grid(global_ewh_grid); %unit: km^2
global_mass_kg_grid = global_ewh_grid.*area_grid.*1e6.*rho_water;%unit: kg
%global_mass_kg_grid = global_ewh_grid.*1.*1e6.*rho_water/1000;% unit: mmWE -> kg; Here the area of each MAR grid is 1 km^2
mass_kg = reshape(global_mass_kg_grid,[nrow*ncol,1]);
mass_kg(zero_index,:)=[];

%     nan_yes = isnan(mass_kg);
%     mass_kg(nan_yes) = [];
%     mass_point_lon = reshape(global_ewh_grid_lon,[nrow*ncol,1]);
%     mass_point_lon(nan_yes) = [];
%     mass_point_lat = reshape(global_ewh_grid_lat,[nrow*ncol,1]);
%     mass_point_lat(nan_yes) = [];

end

function [G_radial,G_tangent_n,G_tangent_e] = get_loading_parameter_grn(calculated_point_lat,calculated_point_lon, mass_point_lat, mass_point_lon,a,grn)
% [ s12, alpha12, alpha21 ] = inverseVincenty(calculated_point_lat, calculated_point_lon, mass_point_lat, mass_point_lon, a, e2 );
azimuth = get_azimuthRad_inside(mass_point_lat, mass_point_lon, calculated_point_lat, calculated_point_lon);
% azimuth = alpha12 - pi 
%
%
%%

% NOTE: 'sind' or 'cosd' here, not 'sin' or 'cos'
dist = acos(sind(calculated_point_lat).*sind(mass_point_lat)+cosd(calculated_point_lat).*cosd(mass_point_lat).*cosd(calculated_point_lon-mass_point_lon));
%
%由于matlab存在舍入误差，有时候上式acos中得到的值会略大于1或者小于-1,那样就出来复数，
dist=real(dist);
%
% loading calculation, depending on spherical distance
%grn = grn1;
theta = (grn(:,1))*pi/180;
% theta = deg2rad(s);

L1_radial(:,1) = grn(:,2);
L2_tangent(:,1) = grn(:,3);

%
%de-normalization of loading
L1_radial = L1_radial ./ (a.*theta.*(10^12));
L2_tangent = L2_tangent ./ (a.*theta.*(10^12));

%
% interpolation of Green's function coefficients to match mean spherical distance
% of spherical segment - 3rd degree spline
% G_radial = interp1(theta,L1_radial,dist,'spline');
% G_tangent = interp1(theta,L2_tangent,dist,'spline');
G_radial = interp1(theta,L1_radial,dist,'linear','extrap');%2022-4-11发现并更换插值方法 %2022-9-19此处加'extrap'可实现线性外推,避免nan
G_tangent = interp1(theta,L2_tangent,dist,'linear','extrap');

% sum of individual segments on deformations
% in direction of meridian and first vertical
G_tangent_n = G_tangent.*(cos(azimuth));
G_tangent_e = G_tangent.*(sin(azimuth));


end

function [ azR ] = get_azimuthRad_inside(mass_lat1,mass_lon1,gps_lat2,gps_lon2)
% 
% Compute the azimut[radian] from 'start' to 'end' points at the location of 'start point'[mass_lat1,mass_lon1] 
% start point: mass_lat,mass_lon
% end point: gps_lat,gps_lon
%
% The value of 'azR'[azR12] is same as 'alpha12' derived from Function
% 'inverseVincenty' bellow:
% [ ~, alpha12, ~] = inverseVincenty(mass_point_lat, mass_point_lon, calculated_point_lat, calculated_point_lon, a, e2 );
%
% The equations refered to the link bellow：
% https://dothinking.github.io/2017-03-08-%E7%90%83%E9%9D%A2%E8%B7%9D%E7%A6%BB%E4%B8%8E%
% E6%96%B9%E4%BD%8D%E8%A7%92%E5%85%AC%E5%BC%8F%E7%9A%84%E6%8E%A8%E5%AF%BC%EF%BC%9A%E8%A7%A3%E4%B8%89%E8%A7%92%E5%BD%A2%E6%B3%95/
%
% azR12:
azR = atan2(cosd(gps_lat2).*sind(gps_lon2-mass_lon1),(cosd(mass_lat1).*sind(gps_lat2)-sind(mass_lat1).*cosd(gps_lat2).*cosd(gps_lon2-mass_lon1)));
% azR21:
% azR = atan2(-cosd(mass_lat1).*sind(gps_lon2-mass_lon1),(cosd(gps_lat2).*sind(mass_lat1)-sind(gps_lat2).*cosd(mass_lat1).*cosd(gps_lon2-mass_lon1)));


% la = lonU2 - lonU1
% from U1 to U2:
% alpha12 = atan2(cos(U2).*sin(la),(cos(U1).*sin(U2)-sin(U1).*cos(U2).*cos(la)));
% from U2 to U1:
% alpha21 = atan2(cos(U1).*sin(la),(-sin(U1).*cos(U2)+cos(U1).*sin(U2).*cos(la)));
end
%

