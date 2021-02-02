if Beta==0
%    delete(['backOutCorrelationsFrict_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
    diary (['backOutCorrelationsFrict_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
else                               
%    delete(['backOutCorrelationsFrict_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
    diary (['backOutCorrelationsFrict_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
end

if Beta==0
    load(['strucFrictHet_',money,'_',num2str(nloc),'loc',sample,'Beta0.mat'],'bstruc');
else
    load(['strucFrictHet_',money,'_',num2str(nloc),'loc',sample,'.mat'],'bstruc');
end

load [REDACTED]sippCombinedNHWmaleAnnual.mat lnpop55
load LLCresults
nodefl = load('[REDACTED]LLCresults.mat','wHat','bwage');
nodefl.bwageCons = nodefl.bwage(1)-.363;
nodefl.wHat = nodefl.wHat+nodefl.bwageCons;
bwageCons = 7.20736301359329;
wHat = wHat+bwageCons;

regdum = zeros(55,9);
regdum(1,5) = 1;
regdum(2,7) = 1;
regdum(3,5) = 1;
regdum(4,1) = 1;
regdum(5,3) = 1;
regdum(6,[3 6]) = 1;
regdum(7,3) = 1;
regdum(8,3) = 1;
regdum(9,7) = 1;
regdum(10,8) = 1;
regdum(11,3) = 1;
regdum(12,7) = 1;
regdum(13,3) = 1;
regdum(14,4) = 1;
regdum(15,6) = 1;
regdum(16,9) = 1;
regdum(17,5) = 1;
regdum(18,3) = 1;
regdum(19,[3 4]) = 1;
regdum(20,[1 2]) = 1;
regdum(21,[2 5]) = 1;
regdum(22,8) = 1;
regdum(23,2) = 1;
regdum(24,9) = 1;
regdum(25,1) = 1;
regdum(26,5) = 1;
regdum(27,9) = 1;
regdum(28,9) = 1;
regdum(29,9) = 1;
regdum(30,9) = 1;
regdum(31,9) = 1;
regdum(32,[3 4]) = 1;
regdum(33,5) = 1;
regdum(34,5) = 1;
regdum(35,5) = 1;
regdum(36,1) = 1;
regdum(37,1) = 1;
regdum(38,2) = 1;
regdum(39,2) = 1;
regdum(40,3) = 1;
regdum(41,3) = 1;
regdum(42,4) = 1;
regdum(43,4) = 1;
regdum(44,5) = 1;
regdum(45,5) = 1;
regdum(46,6) = 1;
regdum(47,6) = 1;
regdum(48,7) = 1;
regdum(49,7) = 1;
regdum(50,8) = 1;
regdum(51,8) = 1;
regdum(52,9) = 1;
regdum(53,9) = 1;
regdum(54,9) = 1;
regdum(55,9) = 1;


%==========================================================================
% Regress lambda, delta on locational characteristics (population and
% Census division dummies)
%==========================================================================
for k=1:9
   eval(['regdum',num2str(k),'=regdum(:,k)*ones(1,10);']);
end
alphal   = reshape(bstruc(1:35)*ones(1,10),350,1);
wagel    = reshape(wHat(1:35,:),350,1);
ndwagel  = reshape(nodefl.wHat(1:35,:),350,1);
deltal   = reshape(deltabar(1:35,:),350,1);
lambdal  = reshape(lambdabar(1:35,:),350,1);
wdriftl  = reshape(bARw(1:35)*ones(1,10),350,1);
wuncerl  = reshape(bARw(nloc+2:nloc+36)*ones(1,10),350,1);
URdriftl = reshape((rho_hat_urate(1,1:35)')*ones(1,10),350,1);
URrhol   = reshape((rho_hat_urate(2,1:35)')*ones(1,10),350,1);
URuncerl = reshape((sig_hat_urate(1,1:35)')*ones(1,10),350,1);
lnpopl   = reshape(lnpop55(1:35,:),350,1);
regdum1l = reshape(regdum1(1:35,:),350,1);
regdum2l = reshape(regdum2(1:35,:),350,1);
regdum3l = reshape(regdum3(1:35,:),350,1);
regdum4l = reshape(regdum4(1:35,:),350,1);
regdum5l = reshape(regdum5(1:35,:),350,1);
regdum6l = reshape(regdum6(1:35,:),350,1);
regdum7l = reshape(regdum7(1:35,:),350,1);
regdum8l = reshape(regdum8(1:35,:),350,1);
regdum9l = reshape(regdum9(1:35,:),350,1);
[ba,sea] = lscov([ones(35,1) lnpopl(1:35) regdum1l(1:35) regdum2l(1:35) regdum3l(1:35) regdum4l(1:35) regdum5l(1:35) regdum6l(1:35) regdum7l(1:35) regdum8l(1:35)],alphal(1:35));
[~,~,~,~,sa] = regress(alphal(1:35),[ones(35,1) lnpopl(1:35) regdum1l(1:35) regdum2l(1:35) regdum3l(1:35) regdum4l(1:35) regdum5l(1:35) regdum6l(1:35) regdum7l(1:35) regdum8l(1:35)]);
[bw,sew] = lscov([ones(350,1) lnpopl regdum1l regdum2l regdum3l regdum4l regdum5l regdum6l regdum7l regdum8l],wagel  );
[~,~,~,~,sw] = regress(wagel  ,[ones(350,1) lnpopl regdum1l regdum2l regdum3l regdum4l regdum5l regdum6l regdum7l regdum8l]);
[bndw,sendw] = lscov([ones(350,1) lnpopl regdum1l regdum2l regdum3l regdum4l regdum5l regdum6l regdum7l regdum8l],ndwagel);
[~,~,~,~,sndw] = regress(ndwagel  ,[ones(350,1) lnpopl regdum1l regdum2l regdum3l regdum4l regdum5l regdum6l regdum7l regdum8l]);
[bd,sed] = lscov([ones(350,1) lnpopl regdum1l regdum2l regdum3l regdum4l regdum5l regdum6l regdum7l regdum8l],deltal );
[~,~,~,~,sd] = regress(deltal ,[ones(350,1) lnpopl regdum1l regdum2l regdum3l regdum4l regdum5l regdum6l regdum7l regdum8l]);
[bl,sel] = lscov([ones(350,1) lnpopl regdum1l regdum2l regdum3l regdum4l regdum5l regdum6l regdum7l regdum8l],lambdal);
[~,~,~,~,sl] = regress(lambdal,[ones(350,1) lnpopl regdum1l regdum2l regdum3l regdum4l regdum5l regdum6l regdum7l regdum8l]);
[bwd,sewd] = lscov([ones(35,1) lnpopl(1:35) regdum1l(1:35) regdum2l(1:35) regdum3l(1:35) regdum4l(1:35) regdum5l(1:35) regdum6l(1:35) regdum7l(1:35) regdum8l(1:35)],wdriftl(1:35));
[~,~,~,~,swd] = regress(wdriftl(1:35),[ones(35,1) lnpopl(1:35) regdum1l(1:35) regdum2l(1:35) regdum3l(1:35) regdum4l(1:35) regdum5l(1:35) regdum6l(1:35) regdum7l(1:35) regdum8l(1:35)]);
[bwu,sewu] = lscov([ones(35,1) lnpopl(1:35) regdum1l(1:35) regdum2l(1:35) regdum3l(1:35) regdum4l(1:35) regdum5l(1:35) regdum6l(1:35) regdum7l(1:35) regdum8l(1:35)],wuncerl(1:35));
[~,~,~,~,swu] = regress(wuncerl(1:35),[ones(35,1) lnpopl(1:35) regdum1l(1:35) regdum2l(1:35) regdum3l(1:35) regdum4l(1:35) regdum5l(1:35) regdum6l(1:35) regdum7l(1:35) regdum8l(1:35)]);
[bURd,seURd] = lscov([ones(35,1) lnpopl(1:35) regdum1l(1:35) regdum2l(1:35) regdum3l(1:35) regdum4l(1:35) regdum5l(1:35) regdum6l(1:35) regdum7l(1:35) regdum8l(1:35)],URdriftl(1:35));
[~,~,~,~,sURd] = regress(URdriftl(1:35),[ones(35,1) lnpopl(1:35) regdum1l(1:35) regdum2l(1:35) regdum3l(1:35) regdum4l(1:35) regdum5l(1:35) regdum6l(1:35) regdum7l(1:35) regdum8l(1:35)]);
[bURr,seURr] = lscov([ones(35,1) lnpopl(1:35) regdum1l(1:35) regdum2l(1:35) regdum3l(1:35) regdum4l(1:35) regdum5l(1:35) regdum6l(1:35) regdum7l(1:35) regdum8l(1:35)],URrhol(1:35));
[~,~,~,~,sURr] = regress(URrhol(1:35),[ones(35,1) lnpopl(1:35) regdum1l(1:35) regdum2l(1:35) regdum3l(1:35) regdum4l(1:35) regdum5l(1:35) regdum6l(1:35) regdum7l(1:35) regdum8l(1:35)]);
[bURu,seURu] = lscov([ones(35,1) lnpopl(1:35) regdum1l(1:35) regdum2l(1:35) regdum3l(1:35) regdum4l(1:35) regdum5l(1:35) regdum6l(1:35) regdum7l(1:35) regdum8l(1:35)],URuncerl(1:35));
[~,~,~,~,sURu] = regress(URuncerl(1:35),[ones(35,1) lnpopl(1:35) regdum1l(1:35) regdum2l(1:35) regdum3l(1:35) regdum4l(1:35) regdum5l(1:35) regdum6l(1:35) regdum7l(1:35) regdum8l(1:35)]);


[ba sea ba./sea]
sa(1)
[bndw sendw bndw./sendw]
sndw(1)
[bw sew bw./sew]
sw(1)
[bd sed bd./sed]
sd(1)
[bl sel bl./sel]
sl(1)
[bwd sewd bwd./sewd]
swd(1)
[bwu sewu bwu./sewu]
swu(1)
[bURd seURd bURd./seURd]
sURd(1)
[bURr seURr bURr./seURr]
sURr(1)
[bURu seURu bURu./seURu]
sURu(1)
togethertable = cat(2,[ba sea ba./sea],[bndw sendw bndw./sendw],[bw sew bw./sew],[bd sed bd./sed],[bl sel bl./sel],[bwd sewd bwd./sewd],[bwu sewu bwu./sewu],[bURd seURd bURd./seURd],[bURr seURr bURr./seURr],[bURu seURu bURu./seURu]);
togetherR2 = cat(2,sa(1),sndw(1),sw(1),sd(1),sl(1),swd(1),swu(1),sURd(1),sURr(1),sURu(1));
dlmwrite('CorrelationsTable.csv',[togethertable;nan(1,size(togethertable,2))]);
dlmwrite('CorrelationsTable.csv',togetherR2,'-append');


%==========================================================================
% Correlation matrices of interest
%==========================================================================
for t=1:10
%     disp(['Year ',num2str(2003+t),':']);
%     [CorrLLM(:,:,t)    ,PLLM(:,:,t)    ]=corr([wHat(:,t) deltabar(:,t) lambdabar(:,t)]);
    [RankCorrLLM(:,:,t),RankPLLM(:,:,t)]=corr([bstruc(1:55) wHat(:,t) deltabar(:,t) lambdabar(:,t) wHat(:,t).*lambdabar(:,t) wHat(:,t).*(1-deltabar(:,t))],'type','Spearman');
end

% [squeeze(RankCorrLLM(1,2,:)) squeeze(RankCorrLLM(1,3,:)) squeeze(RankCorrLLM(1,4,:)) squeeze(RankCorrLLM(2,3,:)) squeeze(RankCorrLLM(2,4,:)) squeeze(RankCorrLLM(3,4,:))]
% [squeeze(RankPLLM(1,2,:)) squeeze(RankPLLM(1,3,:)) squeeze(RankPLLM(1,4,:)) squeeze(RankPLLM(2,3,:)) squeeze(RankPLLM(2,4,:)) squeeze(RankPLLM(3,4,:))]
fprintf('\n %4s %13s %13s %13s %13s %13s %13s %13s %13s\n','Year','Alph,wage','Alph,delta','Alph,lambda','wage,delta','wage,lambda','delta,lambda','w,Ewlam','w,Ewdel');
for t=1:10
    fprintf('%4d %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f \n',2003+t,RankCorrLLM(1,2,t),RankCorrLLM(1,3,t),RankCorrLLM(1,4,t),RankCorrLLM(2,3,t),RankCorrLLM(2,4,t),RankCorrLLM(3,4,t),RankCorrLLM(2,5,t),RankCorrLLM(2,6,t));
    fprintf('%4d %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f \n',2003+t,RankPLLM(1,2,t),RankPLLM(1,3,t),RankPLLM(1,4,t),RankPLLM(2,3,t),RankPLLM(2,4,t),RankPLLM(3,4,t),RankPLLM(2,5,t),RankPLLM(2,6,t));
end

% plot(squeeze(RankCorrLLM(1,2,:)))
% plot(squeeze(RankCorrLLM(1,3,:)))
% plot(squeeze(RankCorrLLM(1,4,:)))
% plot(squeeze(RankCorrLLM(2,3,:)))
% plot(squeeze(RankCorrLLM(2,4,:)))
% plot(squeeze(RankCorrLLM(3,4,:)))

%==========================================================================
% Correlation matrices of interest --- only proper places (no Census divisions)
%==========================================================================
for t=1:10
%     disp(['Year ',num2str(2003+t),':']);
%     [CorrLLM(:,:,t)    ,PLLM(:,:,t)    ]=corr([wHat(:,t) deltabar(:,t) lambdabar(:,t)]);
    [RankCorrLLM(:,:,t),RankPLLM(:,:,t)]=corr([bstruc([1:35 54 55]) wHat([1:35 54 55],t) deltabar([1:35 54 55],t) lambdabar([1:35 54 55],t) wHat([1:35 54 55],t).*lambdabar([1:35 54 55],t) wHat([1:35 54 55],t).*(1-deltabar([1:35 54 55],t))],'type','Spearman');
end

% [squeeze(RankCorrLLM(1,2,:)) squeeze(RankCorrLLM(1,3,:)) squeeze(RankCorrLLM(1,4,:)) squeeze(RankCorrLLM(2,3,:)) squeeze(RankCorrLLM(2,4,:)) squeeze(RankCorrLLM(3,4,:))]
% [squeeze(RankPLLM(1,2,:)) squeeze(RankPLLM(1,3,:)) squeeze(RankPLLM(1,4,:)) squeeze(RankPLLM(2,3,:)) squeeze(RankPLLM(2,4,:)) squeeze(RankPLLM(3,4,:))]
fprintf('\n %4s %13s %13s %13s %13s %13s %13s %13s %13s\n','Year','Alph,wage','Alph,delta','Alph,lambda','wage,delta','wage,lambda','delta,lambda','w,Ewlam','w,Ewdel');
for t=1:10
    fprintf('%4d %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f \n',2003+t,RankCorrLLM(1,2,t),RankCorrLLM(1,3,t),RankCorrLLM(1,4,t),RankCorrLLM(2,3,t),RankCorrLLM(2,4,t),RankCorrLLM(3,4,t),RankCorrLLM(2,5,t),RankCorrLLM(2,6,t));
    fprintf('%4d %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f \n',2003+t,RankPLLM(1,2,t),RankPLLM(1,3,t),RankPLLM(1,4,t),RankPLLM(2,3,t),RankPLLM(2,4,t),RankPLLM(3,4,t),RankPLLM(2,5,t),RankPLLM(2,6,t));
end

%==========================================================================
% Correlation matrices of interest --- only cities (no Census divisions, no AK/HI)
%==========================================================================
for t=1:10
%     disp(['Year ',num2str(2003+t),':']);
    [CorrLLM(:,:,t),    PLLM(:,:,t)    ]=corr([bstruc([1:35]) wHat([1:35],t) deltabar([1:35],t) lambdabar([1:35],t) wHat([1:35],t).*lambdabar([1:35],t) wHat([1:35],t).*(1-deltabar([1:35],t))]);
    [RankCorrLLM(:,:,t),RankPLLM(:,:,t)]=corr([bstruc([1:35]) wHat([1:35],t) deltabar([1:35],t) lambdabar([1:35],t) wHat([1:35],t).*lambdabar([1:35],t) wHat([1:35],t).*(1-deltabar([1:35],t))],'type','Spearman');
end

% [squeeze(RankCorrLLM(1,2,:)) squeeze(RankCorrLLM(1,3,:)) squeeze(RankCorrLLM(1,4,:)) squeeze(RankCorrLLM(2,3,:)) squeeze(RankCorrLLM(2,4,:)) squeeze(RankCorrLLM(3,4,:))]
% [squeeze(RankPLLM(1,2,:)) squeeze(RankPLLM(1,3,:)) squeeze(RankPLLM(1,4,:)) squeeze(RankPLLM(2,3,:)) squeeze(RankPLLM(2,4,:)) squeeze(RankPLLM(3,4,:))]
fprintf('\n %4s %13s %13s %13s %13s %13s %13s %13s %13s\n','Year','Alph,wage','Alph,delta','Alph,lambda','wage,delta','wage,lambda','delta,lambda','w,Ewlam','w,Ewdel');
for t=1:10
    fprintf('%4d %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f \n',2003+t,RankCorrLLM(1,2,t),RankCorrLLM(1,3,t),RankCorrLLM(1,4,t),RankCorrLLM(2,3,t),RankCorrLLM(2,4,t),RankCorrLLM(3,4,t),RankCorrLLM(2,5,t),RankCorrLLM(2,6,t));
    fprintf('%4d %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f \n',2003+t,RankPLLM(1,2,t),RankPLLM(1,3,t),RankPLLM(1,4,t),RankPLLM(2,3,t),RankPLLM(2,4,t),RankPLLM(3,4,t),RankPLLM(2,5,t),RankPLLM(2,6,t));
end

aminotaur = reshape(bstruc(1:35)*ones(1,10),350,1);
wHater = reshape(wHat(1:35,:),350,1);
ndwHater = reshape(nodefl.wHat(1:35,:),350,1);
lambdabarer = reshape(lambdabar(1:35,:),350,1);
deltabarer = reshape(deltabar(1:35,:),350,1);

[RankCorrLLM,RankPLLM]=corr([aminotaur wHater(:) deltabarer(:) lambdabarer(:) wHater(:).*lambdabarer(:) wHater(:).*(1-deltabarer(:))],'type','Spearman');
fprintf('%13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f \n',RankCorrLLM(1,2),RankCorrLLM(1,3),RankCorrLLM(1,4),RankCorrLLM(2,3),RankCorrLLM(2,4),RankCorrLLM(3,4),RankCorrLLM(2,5),RankCorrLLM(2,6));
fprintf('%13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f \n',RankPLLM(1,2),RankPLLM(1,3),RankPLLM(1,4),RankPLLM(2,3),RankPLLM(2,4),RankPLLM(3,4),RankPLLM(2,5),RankPLLM(2,6));
[RankCorrLLM,RankPLLM]=corr([aminotaur ndwHater(:) deltabarer(:) lambdabarer(:) ndwHater(:).*lambdabarer(:) ndwHater(:).*(1-deltabarer(:))],'type','Spearman');
fprintf('%13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f \n',RankCorrLLM(1,2),RankCorrLLM(1,3),RankCorrLLM(1,4),RankCorrLLM(2,3),RankCorrLLM(2,4),RankCorrLLM(3,4),RankCorrLLM(2,5),RankCorrLLM(2,6));
fprintf('%13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f \n',RankPLLM(1,2),RankPLLM(1,3),RankPLLM(1,4),RankPLLM(2,3),RankPLLM(2,4),RankPLLM(3,4),RankPLLM(2,5),RankPLLM(2,6));


return

for t=1:10,temper(t,1)=RankCorrLLM(2,5,t);end
plot(temper)

for t=1:10,temper2(t,1)=RankCorrLLM(2,6,t);end
plot(temper2)

for j=1:nloc
   [temp,ptemp]=corr([wHat(j,:);deltabar(j,:)]','type','Spearman');
   RankCorrTime(j,1)=temp(1,2);
   RankPTime(j,1)   =ptemp(1,2);
end
summarize(RankCorrTime);
[RankCorrTime RankPTime RankPTime<=.1]
RankCorrTime(RankPTime<=.1)

for j=1:nloc
   [temp,ptemp]=corr([wHat(j,:);deltabar(j,:)]');
   RankCorrTime(j,1)=temp(1,2);
   RankPTime(j,1)   =ptemp(1,2);
end
summarize(RankCorrTime);
[RankCorrTime RankPTime RankPTime<=.1]
RankCorrTime(RankPTime<=.1)
