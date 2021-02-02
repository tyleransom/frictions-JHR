%**************************************************************************
% Program outline:
% 1. Estimate wage parameters via OLS
% 2. Estimate flexbile mlogit and generate CCP's
% 3. Using CCP's and E[ln wage], estimate structural flow utility parameters
%**************************************************************************
%==========================================================================
% Import data
%==========================================================================
clear all;
clc;
addpath [REDACTED]

pathstring = pwd;
pathstring = lower(pathstring);
pathstring = pathstring(end-11:end);

% global Beta
nloc   = str2num(pathstring(1:2));
money  = pathstring(7:10); % enter 'wage' or 'earn'
time   = 'annual'; % enter 'annual' or 'trimes'
Beta     =.9;
% Beta     = 0;
sample = upper(pathstring(end-1:end)); % enter 'HS' or 'BA'

modelFitFrict

% mkdir wageBA
% mkdir wageHS
% mkdir earnBA
% mkdir earnHS
% 
% 03loc/earnBA
