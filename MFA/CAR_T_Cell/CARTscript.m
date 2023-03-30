%% CAR-T cell flux and confidence interval estimation script

xmlfile='CARTcell_ccm';
xlsname='EGFRt.xlsx'; %or Rituximab.xlsx
load CAR_T.mat
load EGFRt.mat %or Rituximab.mat


% % account for low variance to avoid overfitting
minvar=0.004^2; % max(0.2*avgvar,1.6e-05);
var=diag(covar);
var(var<minvar)=minvar;
var=diag(var);

%%
input.GLC_IN__C=isotopomervector([1 1 0 0 0 0],1);
input.GLC_IN__C2=isotopomervector([1 1 1 1 1 1],1);
input.GLC_IN__C3=isotopomervector([0 0 0 0 0 0],1);

input.Gln_IN__C=isotopomervector([0 0 0 0 0],1);
input.Gln_IN__C2=isotopomervector([0 0 0 0 0],1);
input.Gln_IN__C3=isotopomervector([0 0 0 0 0 ; 1 1 1 1 1], [0.43 0.57]);

input.Glu_IN__C=isotopomervector([0 0 0 0 0],1);
input.Glu_IN__C2=isotopomervector([0 0 0 0 0],1);
input.Glu_IN__C3=isotopomervector([0 0 0 0 0],1);

input.AC_IN__C=isotopomervector([0 0],1);
input.AC_IN__C2=isotopomervector([0 0],1);
input.AC_IN__C3=isotopomervector([0 0],1);

input.OAA_IN__C=isotopomervector([0 0 0 0],1);
input.OAA_IN__C2=isotopomervector([0 0 0 0],1);
input.OAA_IN__C3=isotopomervector([0 0 0 0],1);

input.CO2_IN__C=isotopomervector(0,1);
input.CO2_IN__C2=isotopomervector(0,1);
input.CO2_IN__C3=isotopomervector(0,1);

free_net=rand(size(model.kernel_net,2),1);
free_xch=rand(size(model.kernel_xch,2),1);

[fmea,fmeaStr]=xlsreadfmea(xlsname,model,free_net,free_xch);
[ineq,ineqStr]=xlsreadineq(xlsname,model,free_net,free_xch);  

simulate=str2func(xmlfile);

%% evaluate the best fit fluxes and labeling pattern
[net_opt,xch_opt,info,net_,xch_,info_,fval_,I]=pargloptflux(simulate,model,free_net,free_xch,ineq,[],input,mea,fmea,var,12);

[tflux,tiso,tscore]=writetableplus(simulate,model,net_opt,xch_opt,input,mea,fmea,var);
sim=simulate(net_opt,xch_opt,input);

%% estimate confidence interval
[lb,ub,hs,net_nopt,xch_nopt,nscore,nflag,I]=parconfest3_1(simulate,model,net_opt,xch_opt,ineq,input,mea,fmea,var);
