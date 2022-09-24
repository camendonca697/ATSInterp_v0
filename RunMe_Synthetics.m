%
% THIN SHEET AUTOMATIC INTERPRETATION USING REGULARIZED SECOND DERIVATIVES   
%
% SIMULATION WITH 2D PRISMATIC MODELS
%
clear
clc
help RunMe_Synthetics
% 1) profile parameters (distance in km)
zo=0;   % sensor height above the ground
xi=0.0; % profile initial position
xf=100; % profile final position
dd=0.35; % sampling interval
x0=[xi:dd:xf]';
nx=length(x0);
% 2) prism parameters
zT=4;   % reference to depth to the top 
HH=24;  % prisms depth to the bottom 
b=2*dd; % reference prism width
zmax=24; % maximum depth for model representation
%
% 3) Geomagnetic and magnetization parameters
Tinc=30;    % inclination of geomagnetic field
MincI=Tinc; % induced magnetization inclination
MincR=90;   % remanent magnetization inclination
%
% 4) testing models (see in aux01_TestingModels.m)
sgm=2; % noise level to corrupt data (nT)
for nMOD=1:10
aux01_TestingModels
    % A) model response
    [H Z ftt]=multiprism(x0,pp);
    TT=sqrt(H.*H+Z.*Z);
    TTr=TT+randn(nx,1)*sgm;
    % step 1) x-derivative from NON-CORRUPTED DATA
    gz2=zeros(nx,1);
    aR2 = diff(TT(:),2)/(dd*dd);
    gz2(:,1)=[0;aR2;0];
    % step 2) x-derivative from CORRUPTED DATA (non-regularized)
    gz2N=zeros(nx,1); % padding zeros to fit vector length
    aR2 = diff(TTr(:),2)/(dd*dd);
    gz2N(:,1)=[0;aR2;0];
    % step 3) x-derivative from CORRUPTED DATA (regularized)
    aR1 = rdiff(TTr(:),dd,'tikhonov');
    gz2R= rdiff(aR1,dd,'tikhonov');
    aux03_ZaPlotting
    pause(1)
end

