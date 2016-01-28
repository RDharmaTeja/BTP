function [] = postprocess(proc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
X = load(strcat('process_xcoord_',num2str(proc),'.txt'));
Y = load(strcat('process_ycoord_',num2str(proc),'.txt'));
P = load(strcat('pressure_',num2str(proc),'.txt'));
mesh(X,Y,P);
end

