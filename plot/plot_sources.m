% load true sources if available (synthetic test)
src_true = load('../sources_true.dat');

% load initial sources
src_init_ = dlmread('../bsources.dat');
src_init = src_init_(2:end,:);
nsrc = src_init_(1,1);

% load mean and variance
chains=[5];
ave=zeros(src_init_(1,1),src_init_(1,2));
var=zeros(src_init_(1,1),src_init_(1,2));
ave_ = ave; var_ =var;
ntotal = 0;
for i = chains
    ave_tmp = dlmread(['../Results/resume/run_time_mean_sources_' num2str(i) '.dat']);
    var_tmp = dlmread(['../Results/resume/run_time_var_sources_' num2str(i) '.dat']);
    ave_(:,1:3) = ave_tmp(2:2:end,1:3); ave_(:,4) = ave_tmp(3:2:end,1);
    var_(:,1:3) = var_tmp(2:2:end,1:3); var_(:,4) = var_tmp(3:2:end,1);
    num = dlmread(['../Results/resume/run_time_info_' num2str(i) '.dat']);
    nthin = num(3,1); 
    ave = ave + ave_.*nthin;
    var = var + var_.*nthin;
    ntotal = ntotal + nthin;
end
ave = ave./ntotal;
var = var./ntotal;
std = real((var - ave.^2).^0.5);

%%
err_init = src_init - src_true;
err = src_init + ave - src_true;
errL = -std/2;
errU = std/2;
figure;
errorbar((1:nsrc),err(:,4),errL(:,4),errU(:,4),'ro');
hold on; plot((1:nsrc),err_init(:,4),'b*');

figure;
errorbar((1:nsrc),err(:,1),errL(:,1),errU(:,1),'ro');
hold on; plot((1:nsrc),err_init(:,1),'b*');
figure;
errorbar((1:nsrc),err(:,2),errL(:,2),errU(:,2),'ro');
hold on; plot((1:nsrc),err_init(:,2),'b*');
figure;
errorbar((1:nsrc),err(:,3),errL(:,3),errU(:,3),'ro');
hold on; plot((1:nsrc),err_init(:,3),'b*');