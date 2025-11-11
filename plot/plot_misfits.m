% plot the number of cells of markov chains
%
%

% load number of cells
n = 16;
misfit = cell(1,n);
%chains = [1:1:15,17 18 19 20];
chains = [1:8];
%chains = [ 15 11 10];
%%
x=(1:length(misfit{1}));
for i = chains
    misfit{i} = load(['../' num2str(i) '/Results/likelihood_' num2str(i) '.dat']);
%     misfit{i} = misfit{i}(:,2);
    %figure;
    %x=(1:length(misfit{i}));
    %semilogy(x,misfit{i}(:,1));
%     misfit{i} = load(['misfit.out' num2str(i)]);
end

%% plot
x=(1:length(misfit{1}));
len = length(misfit{1});
figure;

% hold on;
for i= chains
%     figure;
    x=(1:length(misfit{i}));
    semilogy(x,misfit{i}(:,2));
%       plot(x,misfit{i}(:,2));
    hold on;
end
%axis([0 len 0 2]);
set(gca,'fontsize',14);
% title('misfits');
xlabel('Samples');
ylabel('Misfits');
figure;
%%
% hold on;
for i= chains
    x=(1:length(misfit{i}));
    semilogy(x,misfit{i}(:,3));
%     plot(log(x),log(misfit{i}(:,1)));
     hold on;
end
% title('likelihood');
% loglog(x,misfit{1}(1:len,1),x,misfit{2}(1:len,1),x,misfit{3}(1:len,1),...
%     x,misfit{4}(1:len,1),x,misfit{5}(1:len,1),x,misfit{6}(1:len,1),...
%     x,misfit{8}(1:len,1));
% loglog(x,misfit{1}(1:len,1),x,misfit{2}(1:len,1),x,misfit{3}(1:len,1),...
%     x,misfit{4}(1:len,1),x,misfit{5}(1:len,1),x,misfit{6}(1:len,1),...
%     x,misfit{7}(1:len,1),x,misfit{8}(1:len,1));
% set(gca,'fontsize',16);
% axis([0 len 0 0.02]);
set(gca,'fontsize',14);
xlabel('Samples');
ylabel('Unweighted');
% print('misfits.png','-dpng','-r300');

%%
figure;
% hold on;
for i= chains
    x=(1:length(misfit{i}));
    sigma = (misfit{i}(:,2)./misfit{i}(:,1)./2).^0.5;
    %semilogy(x,misfit{i}(:,2));
    plot(x,sigma);
    hold on;
end
set(gca,'fontsize',14);
% title('sigmas');
xlabel('Samples');
ylabel('sigmas');
