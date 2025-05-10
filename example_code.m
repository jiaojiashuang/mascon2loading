% ------- example code for 'get_mascon2loading_scatter_Dense_Global'--------%
%
% written by,
% Jiashuang Jiao
% 2023-10-25
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.load the necessary data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Green's Function
load LoadExampleData.mat GreensFunction_prem_CF
% load longitude and Latitude coordinates of GNSS stations
load LoadExampleData.mat gps_lon_lat
% load EWH mass anomaly field, Unit:m
load LoadExampleData.mat ewh_jpl_m_200501_200512
% load time nodes [not necessary for calculation]
load LoadExampleData.mat time_jpl_200501_200512

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.Calculate the 3-D load deformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GreensFunction=GreensFunction_prem_CF;
grid_EWH_m=ewh_jpl_m_200501_200512;
gps_lon=gps_lon_lat(:,1);
gps_lat=gps_lon_lat(:,2);
%
number_gps = length(gps_lon);%the number of calculated points
number_month = size(grid_EWH_m,3);%the number of time nodes
%
for iii = 1:number_gps

    [loading_meter_u_global_temp,loading_meter_n_global_temp,loading_meter_e_global_temp] = get_mascon2loading_scatter_Dense_Global(GreensFunction,grid_EWH_m,gps_lon(iii,1),gps_lat(iii,1),50,10); %
    loading_meter_u_global(:,iii) = loading_meter_u_global_temp;% [number_month X number_gps]
    loading_meter_n_global(:,iii) = loading_meter_n_global_temp;% [number_month X number_gps]
    loading_meter_e_global(:,iii) = loading_meter_e_global_temp;% [number_month X number_gps]
    %
    fprintf('%d %s\n',iii,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('%d %d %d %s\n',iii,iii,iii,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('%d %d %d %d %d %s\n',iii,iii,iii,iii,iii,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('%d %d %d %d %d %d %d %s\n',iii,iii,iii,iii,iii,iii,iii,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('%d %d %d %d %d %d %d %d %d %s\n',iii,iii,iii,iii,iii,iii,iii,iii,iii,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('%d %d %d %d %d %d %d %d %d %d %d %s\n',iii,iii,iii,iii,iii,iii,iii,iii,iii,iii,iii,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('%d %d %d %d %d %d %d %d %d %d %d %d %d %s\n',iii,iii,iii,iii,iii,iii,iii,iii,iii,iii,iii,iii,iii,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
end


