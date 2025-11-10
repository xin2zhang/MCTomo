% plot the number of cells of markov chains
%
%

% load number of cells
n = 8;
ncell = cell(1,n);
%chains = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20];
%chains = [15 11 10];
chains = [1 2 3 4 5 6 7 8];
for i = chains
    ncell{i} = load(['../' num2str(i) '/Results/ncells_' num2str(i) '.dat']);
end
%% xlike
x=(1:length(ncell{1}));
xlike = [1000:1000:1200000];
figure;hold on;
for i= chains
%     figure;
    nlike_chain = x.*0.95;
    plot(nlike_chain,ncell{i});
    hold on;
end
% load time
time=load('time.dat');
time = time;
%% plot
x=(1:length(ncell{1}));
len = length(ncell{1});
figure;hold on;
for i= chains
%     figure;
    x=(1:length(ncell{i}))*time(i)/length(ncell{i});
    plot(x,ncell{i}(:,1));
    hold on;
end

% plot(x,ncell{1}(1:len,1),x,ncell{2}(1:len,1),x,ncell{3}(1:len,1),...
%     x,ncell{4}(1:len,1),x,ncell{5}(1:len,1),x,ncell{6}(1:len,1),...
%     x,ncell{7}(1:len,1),x,ncell{8}(1:len,1),'lineWidth',0.5);
set(gca,'fontsize',14);
box on;
% axis([0 1.3e6 0 110]);
xlabel('Time(s)','fontsize',14);
ylabel('Number of cells','fontsize',14);
%print('ncells.png','-dpng','-r300');
figure;hold on;
for i= chains
%     figure;
    x=(1:length(ncell{i}));
    plot(x,ncell{i}(:,1));
    hold on;
end
xlabel('Samples');
ylabel('Number of cells');
