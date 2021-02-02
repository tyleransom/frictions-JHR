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
nloc   = str2double(pathstring(1:2));
money  = pathstring(7:10); % enter 'wage' or 'earn'
time   = 'annual'; % enter 'annual' or 'trimes'
Beta     =.9;
% Beta     = 0;
sample = upper(pathstring(end-1:end)); % enter 'HS' or 'BA'

locname{1,1}='Atl';locname{2,1}='Aus';locname{3,1}='Bal';locname{4,1}='Bos';locname{5,1}='Chi';locname{6,1}='Cin';locname{7,1}='Cle';locname{8,1}='Col';locname{9,1}='Dal';locname{10,1}='Den';locname{11,1}='Det';locname{12,1}='Hou';locname{13,1}='Ind';locname{14,1}='Kan';locname{15,1}='Kno';locname{16,1}='Los';locname{17,1}='Mia';locname{18,1}='Mil';locname{19,1}='Min';locname{20,1}='NYC';locname{21,1}='Phi';locname{22,1}='Phx';locname{23,1}='Pit';locname{24,1}='PDX';locname{25,1}='Pro';locname{26,1}='Ric';locname{27,1}='Riv';locname{28,1}='Sac';locname{29,1}='SDG';locname{30,1}='SFO';locname{31,1}='Sea';locname{32,1}='STL';locname{33,1}='Tam';locname{34,1}='Vir';locname{35,1}='WDC';locname{36,1}='Nes';locname{37,1}='Nem';locname{38,1}='Mas';locname{39,1}='Mam';locname{40,1}='Ens';locname{41,1}='Enm';locname{42,1}='Wns';locname{43,1}='Wnm';locname{44,1}='Sas';locname{45,1}='Sam';locname{46,1}='Ess';locname{47,1}='Esm';locname{48,1}='Wss';locname{49,1}='Wsm';locname{50,1}='Mts';locname{51,1}='Mtm';locname{52,1}='Pas';locname{53,1}='Pam';locname{54,1}='Ala';locname{55,1}='Haw';
stabb{1,1}='AL';stabb{2,1}='AK';stabb{4,1}='AZ';stabb{5,1}='AR';stabb{6,1}='CA';stabb{8,1}='CO';stabb{9,1}='CT';stabb{10,1}='DE';stabb{11,1}='DC';stabb{12,1}='FL';stabb{13,1}='GA';stabb{15,1}='HI';stabb{16,1}='ID';stabb{17,1}='IL';stabb{18,1}='IN';stabb{19,1}='IA';stabb{20,1}='KS';stabb{21,1}='KY';stabb{22,1}='LA';stabb{23,1}='ME';stabb{24,1}='MD';stabb{25,1}='MA';stabb{26,1}='MI';stabb{27,1}='MN';stabb{28,1}='MS';stabb{29,1}='MO';stabb{30,1}='MT';stabb{31,1}='NE';stabb{32,1}='NV';stabb{33,1}='NH';stabb{34,1}='NJ';stabb{35,1}='NM';stabb{36,1}='NY';stabb{37,1}='NC';stabb{38,1}='ND';stabb{39,1}='OH';stabb{40,1}='OK';stabb{41,1}='OR';stabb{42,1}='PA';stabb{44,1}='RI';stabb{45,1}='SC';stabb{46,1}='SD';stabb{47,1}='TN';stabb{48,1}='TX';stabb{49,1}='UT';stabb{50,1}='VT';stabb{51,1}='VA';stabb{53,1}='WA';stabb{54,1}='WV';stabb{55,1}='WI';stabb{56,1}='WY';
birthStateDivXwalk(1,:) = [6 15 46 47 47 47 47 47 47 47 47];birthStateDivXwalk(2,:) = [16 24 27 28 29 30 31 52 53 54 55];birthStateDivXwalk(4,:) = [10 22 50 51 51 51 51 51 51 51 51];birthStateDivXwalk(5,:) = [2 9 12 48 49 49 49 49 49 49 49];birthStateDivXwalk(6,:) = [16 24 27 28 29 30 31 52 53 54 55];birthStateDivXwalk(8,:) = [10 22 50 51 51 51 51 51 51 51 51];birthStateDivXwalk(9,:) = [4 20 25 36 37 37 37 37 37 37 37];birthStateDivXwalk(10,:) = [1 3 17 21 26 33 34 35 44 45 45];birthStateDivXwalk(11,:) = [1 3 17 21 26 33 34 35 44 45 45];birthStateDivXwalk(12,:) = [1 3 17 21 26 33 34 35 44 45 45];birthStateDivXwalk(13,:) = [1 3 17 21 26 33 34 35 44 45 45];birthStateDivXwalk(15,:) = [16 24 27 28 29 30 31 52 53 54 55];birthStateDivXwalk(16,:) = [10 22 50 51 51 51 51 51 51 51 51];birthStateDivXwalk(17,:) = [5 6 7 8 11 13 18 19 32 40 41];birthStateDivXwalk(18,:) = [5 6 7 8 11 13 18 19 32 40 41];birthStateDivXwalk(19,:) = [14 19 32 42 43 43 43 43 43 43 43];birthStateDivXwalk(20,:) = [14 19 32 42 43 43 43 43 43 43 43];birthStateDivXwalk(21,:) = [6 15 46 47 47 47 47 47 47 47 47];birthStateDivXwalk(22,:) = [2 9 12 48 49 49 49 49 49 49 49];birthStateDivXwalk(23,:) = [4 20 25 36 37 37 37 37 37 37 37];birthStateDivXwalk(24,:) = [1 3 17 21 26 33 34 35 44 45 45];birthStateDivXwalk(25,:) = [4 20 25 36 37 37 37 37 37 37 37];birthStateDivXwalk(26,:) = [5 6 7 8 11 13 18 19 32 40 41];birthStateDivXwalk(27,:) = [14 19 32 42 43 43 43 43 43 43 43];birthStateDivXwalk(28,:) = [6 15 46 47 47 47 47 47 47 47 47];birthStateDivXwalk(29,:) = [14 19 32 42 43 43 43 43 43 43 43];birthStateDivXwalk(30,:) = [10 22 50 51 51 51 51 51 51 51 51];birthStateDivXwalk(31,:) = [14 19 32 42 43 43 43 43 43 43 43];birthStateDivXwalk(32,:) = [10 22 50 51 51 51 51 51 51 51 51];birthStateDivXwalk(33,:) = [4 20 25 36 37 37 37 37 37 37 37];birthStateDivXwalk(34,:) = [20 21 23 38 39 39 39 39 39 39 39];birthStateDivXwalk(35,:) = [10 22 50 51 51 51 51 51 51 51 51];birthStateDivXwalk(36,:) = [20 21 23 38 39 39 39 39 39 39 39];birthStateDivXwalk(37,:) = [1 3 17 21 26 33 34 35 44 45 45];birthStateDivXwalk(38,:) = [14 19 32 42 43 43 43 43 43 43 43];birthStateDivXwalk(39,:) = [5 6 7 8 11 13 18 19 32 40 41];birthStateDivXwalk(40,:) = [2 9 12 48 49 49 49 49 49 49 49];birthStateDivXwalk(41,:) = [16 24 27 28 29 30 31 52 53 54 55];birthStateDivXwalk(42,:) = [20 21 23 38 39 39 39 39 39 39 39];birthStateDivXwalk(44,:) = [4 20 25 36 37 37 37 37 37 37 37];birthStateDivXwalk(45,:) = [1 3 17 21 26 33 34 35 44 45 45];birthStateDivXwalk(46,:) = [14 19 32 42 43 43 43 43 43 43 43];birthStateDivXwalk(47,:) = [6 15 46 47 47 47 47 47 47 47 47];birthStateDivXwalk(48,:) = [2 9 12 48 49 49 49 49 49 49 49];birthStateDivXwalk(49,:) = [10 22 50 51 51 51 51 51 51 51 51];birthStateDivXwalk(50,:) = [4 20 25 36 37 37 37 37 37 37 37];birthStateDivXwalk(51,:) = [1 3 17 21 26 33 34 35 44 45 45];birthStateDivXwalk(53,:) = [16 24 27 28 29 30 31 52 53 54 55];birthStateDivXwalk(54,:) = [1 3 17 21 26 33 34 35 44 45 45];birthStateDivXwalk(55,:) = [5 6 7 8 11 13 18 19 32 40 41];birthStateDivXwalk(56,:) = [10 22 50 51 51 51 51 51 51 51 51];
birthStateLocXwalk(1,:) = [46 47 47 47 47 47 47];birthStateLocXwalk(2,:) = [54 54 54 54 54 54 54];birthStateLocXwalk(4,:) = [22 50 51 51 51 51 51];birthStateLocXwalk(5,:) = [48 49 49 49 49 49 49];birthStateLocXwalk(6,:) = [16 27 28 29 30 52 53];birthStateLocXwalk(8,:) = [10 50 51 51 51 51 51];birthStateLocXwalk(9,:) = [20 36 37 37 37 37 37];birthStateLocXwalk(10,:) = [21 44 45 45 45 45 45];birthStateLocXwalk(11,:) = [35 44 45 45 45 45 45];birthStateLocXwalk(12,:) = [17 33 44 45 45 45 45];birthStateLocXwalk(13,:) = [1 44 45 45 45 45 45];birthStateLocXwalk(15,:) = [55 55 55 55 55 55 55];birthStateLocXwalk(16,:) = [50 51 51 51 51 51 51];birthStateLocXwalk(17,:) = [5 32 40 41 41 41 41];birthStateLocXwalk(18,:) = [5 6 13 40 41 41 41];birthStateLocXwalk(19,:) = [42 43 43 43 43 43 43];birthStateLocXwalk(20,:) = [14 42 43 43 43 43 43];birthStateLocXwalk(21,:) = [6 46 47 47 47 47 47];birthStateLocXwalk(22,:) = [48 49 49 49 49 49 49];birthStateLocXwalk(23,:) = [36 37 37 37 37 37 37];birthStateLocXwalk(24,:) = [3 35 44 45 45 45 45];birthStateLocXwalk(25,:) = [4 25 36 37 37 37 37];birthStateLocXwalk(26,:) = [11 40 41 41 41 41 41];birthStateLocXwalk(27,:) = [19 42 43 43 43 43 43];birthStateLocXwalk(28,:) = [46 47 47 47 47 47 47];birthStateLocXwalk(29,:) = [14 32 42 43 43 43 43];birthStateLocXwalk(30,:) = [50 51 51 51 51 51 51];birthStateLocXwalk(31,:) = [42 43 43 43 43 43 43];birthStateLocXwalk(32,:) = [50 51 51 51 51 51 51];birthStateLocXwalk(33,:) = [4 36 37 37 37 37 37];birthStateLocXwalk(34,:) = [20 21 38 39 39 39 39];birthStateLocXwalk(35,:) = [50 51 51 51 51 51 51];birthStateLocXwalk(36,:) = [20 38 39 39 39 39 39];birthStateLocXwalk(37,:) = [44 45 45 45 45 45 45];birthStateLocXwalk(38,:) = [42 43 43 43 43 43 43];birthStateLocXwalk(39,:) = [6 7 8 40 41 41 41];birthStateLocXwalk(40,:) = [48 49 49 49 49 49 49];birthStateLocXwalk(41,:) = [24 52 53 53 53 53 53];birthStateLocXwalk(42,:) = [20 21 23 38 39 39 39];birthStateLocXwalk(44,:) = [25 36 37 37 37 37 37];birthStateLocXwalk(45,:) = [44 45 45 45 45 45 45];birthStateLocXwalk(46,:) = [42 43 43 43 43 43 43];birthStateLocXwalk(47,:) = [15 46 47 47 47 47 47];birthStateLocXwalk(48,:) = [2 9 12 48 49 49 49];birthStateLocXwalk(49,:) = [50 51 51 51 51 51 51];birthStateLocXwalk(50,:) = [36 37 37 37 37 37 37];birthStateLocXwalk(51,:) = [26 34 35 44 45 45 45];birthStateLocXwalk(53,:) = [24 31 52 53 53 53 53];birthStateLocXwalk(54,:) = [44 45 45 45 45 45 45];birthStateLocXwalk(55,:) = [5 18 19 40 41 41 41];birthStateLocXwalk(56,:) = [50 51 51 51 51 51 51];
birthStateDivXwalk = [birthStateDivXwalk birthStateDivXwalk+nloc];
birthStateLocXwalk = [birthStateLocXwalk birthStateLocXwalk+nloc];
filestring = mfilename
if isempty(filestring)
   filestring=getenv('MFNAME') 
end

typer        = str2double(filestring(25:25))
yr           = str2double(filestring(14:17))
locstring    = filestring(11:13)
ager         = str2double(filestring(9:10))
bstatestring = filestring(22:23)
diarystub    = filestring(1:25);

LYtemp                   = find(strcmpi(locname,locstring)==1)
birthStateFIPS           = find(strcmpi(stabb,bstatestring)==1);
birthLocVec              = birthStateLocXwalk(birthStateFIPS,:);
birthDivVec              = birthStateDivXwalk(birthStateFIPS,:);
birthLoci                = zeros(1,2*nloc);
birthLoci(1,birthLocVec) = 1;
birthDivi                = zeros(1,2*nloc);
birthDivi(1,birthDivVec) = 1;

susumeLocs = [6 10 29];

cflFrictS

% mkdir wageBA
% mkdir wageHS
% mkdir earnBA
% mkdir earnHS
% 
% 03loc/earnBA
