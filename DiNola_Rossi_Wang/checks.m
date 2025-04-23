clear
clc
close all

% StationaryDist_FHorz_Case1_e
% StationaryDist_FHorz_Case1_Iteration_e_raw
%  0.999999998835758

load jone_cpu.mat
load jone_gpu.mat

% Reshape jone_gpu(a,z,e) ==> jone_gpu(a,z1,z2,e)
n_a  = size(jone_cpu,1);
n_e  = size(jone_cpu,4);
n_z1 = 7;
n_z2 = 2;
jone_gpu = reshape(jone_gpu,[n_a,n_z1,n_z2,n_e]);

err = max(abs(jone_gpu-jone_cpu),[],"all")