%==========================================================================
% Import data
%==========================================================================
clear all;
clc;
addpath [REDACTED]

load(['Results25hmm2007BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMuhmm07t1,seHMuhmm07t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25hmm2011BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMuhmm11t1,seHMuhmm11t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25hmm2007BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMuhmm07t0,seHMuhmm07t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25hmm2011BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMuhmm11t0,seHMuhmm11t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));

load(['Results25lmm2007BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMulmm07t1,seHMulmm07t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25lmm2011BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMulmm11t1,seHMulmm11t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25lmm2007BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMulmm07t0,seHMulmm07t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25lmm2011BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMulmm11t0,seHMulmm11t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));

load(['Results25mhm2007BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMumhm07t1,seHMumhm07t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25mhm2011BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMumhm11t1,seHMumhm11t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25mhm2007BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMumhm07t0,seHMumhm07t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25mhm2011BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMumhm11t0,seHMumhm11t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));

load(['Results25mlm2007BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMumlm07t1,seHMumlm07t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25mlm2011BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMumlm11t1,seHMumlm11t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25mlm2007BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMumlm07t0,seHMumlm07t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25mlm2011BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMumlm11t0,seHMumlm11t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));

load(['Results25mmh2007BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMummh07t1,seHMummh07t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25mmh2011BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMummh11t1,seHMummh11t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25mmh2007BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMummh07t0,seHMummh07t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25mmh2011BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMummh11t0,seHMummh11t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));

load(['Results25mml2007BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMumml07t1,seHMumml07t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25mml2011BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMumml11t1,seHMumml11t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25mml2007BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMumml07t0,seHMumml07t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25mml2011BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMumml11t0,seHMumml11t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));

load(['Results25hhh2007BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMuhhh07t1,seHMuhhh07t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25hhh2011BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMuhhh11t1,seHMuhhh11t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25hhh2007BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMuhhh07t0,seHMuhhh07t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25hhh2011BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMuhhh11t0,seHMuhhh11t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));

load(['Results25lll2007BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMulll07t1,seHMulll07t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25lll2011BornILt1'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMulll11t1,seHMulll11t1] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25lll2007BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMulll07t0,seHMulll07t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
load(['Results25lll2011BornILt0'],'HeatMap*','Covar*','nloc','LYtemp','yr');
[bHMulll11t0,seHMulll11t0] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));

dlmwrite(['DirectionalRegressionsOutput.csv'],[bHMuhmm07t1 seHMuhmm07t1 bHMulmm07t1 seHMulmm07t1 bHMuhmm11t1 seHMuhmm11t1 bHMulmm11t1 seHMulmm11t1 bHMuhmm07t0 seHMuhmm07t0 bHMulmm07t0 seHMulmm07t0 bHMuhmm11t0 seHMuhmm11t0 bHMulmm11t0 seHMulmm11t0]);
dlmwrite(['DirectionalRegressionsOutput.csv'],[nan(1,16)],'-append');
dlmwrite(['DirectionalRegressionsOutput.csv'],[bHMumhm07t1 seHMumhm07t1 bHMumlm07t1 seHMumlm07t1 bHMumhm11t1 seHMumhm11t1 bHMumlm11t1 seHMumlm11t1 bHMumhm07t0 seHMumhm07t0 bHMumlm07t0 seHMumlm07t0 bHMumhm11t0 seHMumhm11t0 bHMumlm11t0 seHMumlm11t0],'-append');
dlmwrite(['DirectionalRegressionsOutput.csv'],[nan(1,16)],'-append');
dlmwrite(['DirectionalRegressionsOutput.csv'],[bHMummh07t1 seHMummh07t1 bHMumml07t1 seHMumml07t1 bHMummh11t1 seHMummh11t1 bHMumml11t1 seHMumml11t1 bHMummh07t0 seHMummh07t0 bHMumml07t0 seHMumml07t0 bHMummh11t0 seHMummh11t0 bHMumml11t0 seHMumml11t0],'-append');
dlmwrite(['DirectionalRegressionsOutput.csv'],[nan(1,16)],'-append');
dlmwrite(['DirectionalRegressionsOutput.csv'],[bHMuhhh07t1 seHMuhhh07t1 bHMulll07t1 seHMulll07t1 bHMuhhh11t1 seHMuhhh11t1 bHMulll11t1 seHMulll11t1 bHMuhhh07t0 seHMuhhh07t0 bHMulll07t0 seHMulll07t0 bHMuhhh11t0 seHMuhhh11t0 bHMulll11t0 seHMulll11t0],'-append');
dlmwrite(['DirectionalRegressionsOutput.csv'],[nan(1,16)],'-append');


