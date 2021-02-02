%**************************************************************************
% Program outline:
% 1. Compute correlations between offers, layoffs, wages, and amenities
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

backOutCorrelationsFrict
