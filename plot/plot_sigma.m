% plot the data noise level of markov chains
%
%

% load sigmas
sigmas = cell(1,n);
chains = [1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 18 19 20];
%chains = [ 15 11 10];
%%
for i = chains
    sigmas{i} = load(['../' num2str(i) '/sigma_' num2str(i) '.dat']);
%     misfit{i} = misfit{i}(:,2);
    %figure;
    %x=(1:length(misfit{i}));
    %semilogy(x,misfit{i}(:,1));
%     misfit{i} = load(['misfit.out' num2str(i)]);
end

%% plot
x=(1:length(sigmas{1}));
figure;
for i= chains
    len = length(sigmas{i});
    x=(1:len/2);
    semilogy(x,sigmas{i}(1:2:len,3));
%     plot(log(x),log(misfit{i}(:,1)));
    hold on;
end
% loglog(x,misfit{1}(1:len,1),x,misfit{2}(1:len,1),x,misfit{3}(1:len,1),...
%     x,misfit{4}(1:len,1),x,misfit{5}(1:len,1),x,misfit{6}(1:len,1),...
%     x,misfit{8}(1:len,1));
% loglog(x,misfit{1}(1:len,1),x,misfit{2}(1:len,1),x,misfit{3}(1:len,1),...
%     x,misfit{4}(1:len,1),x,misfit{5}(1:len,1),x,misfit{6}(1:len,1),...
%     x,misfit{7}(1:len,1),x,misfit{8}(1:len,1));
set(gca,'fontsize',16);
% axis([0 len 0 0.02]);
xlabel('Samples');
ylabel('misfit');
% print('misfits.png','-dpng','-r300');
