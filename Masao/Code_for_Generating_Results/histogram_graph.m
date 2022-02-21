clear
close all

est = csvread('betasigma_record.csv');
est = est(1,:);
bh = zeros(1,25);
bh(1) = est(1);
bh(2) = est(2);
bh(3) = est(3);
bh(4) = est(7);
bh(5) = est(8);
bh(6) = est(13);
bh(7) = est(15);
bh(8) = est(4);
bh(9) = est(12);
bh(10) = est(9);
bh(11) = est(10).*est(9);
bh(12) = (est(10).*est(9)./10).^(est(12)+1);
bh(13) = est(11).*est(9);
bh(14) = (est(11).*est(9)./10).^(est(12)+1);
bh(15) = est(5);
bh(16) = (est(5)./5).^(est(12)+1);
bh(17) = est(6);
bh(18) = est(16);
bh(19) = est(17);
bh(20) = est(18);
bh(21) = est(19);
bh(22) = est(14);
bh(23) = est(20);
bh(24) = est(20).*est(21);
bh(25) = 1-est(22);

tr = bh;
boot = csvread('bootstrap_estimates.csv');
boot = unique(boot,'rows');
boot = boot(:,2:end);
distpara = csvread('bootstrap_betasigma_.csv');
distpara = unique(distpara,'rows');
distpara = distpara(:,2:7);
%%
b_boot = zeros(size(boot,1),25);
for i = 1:size(boot,1)
    distpara0 = distpara(i,:);
    xx = boot(i,:);
    est = [distpara0...
        [xx(1) xx(2)/(1+xx(1)) xx(3) ...
        xx(4)*10*xx(7)/10/9.5^(-xx(6)-1) xx(5)*10*xx(7)/10/8^(-xx(6)-1) ...
        xx(6:11).*[1 0.1 1 0.1 0.01 0.1 ] 0 0] xx(12:13) xx(14)];
    bh(1) = est(1);
    bh(2) = est(2);
    bh(3) = est(3);
    bh(4) = est(7);
    bh(5) = est(8);
    bh(6) = est(13);
    bh(7) = est(15);
    bh(8) = est(4);
    bh(9) = est(12);
    bh(10) = est(9);
    bh(11) = est(10).*est(9);
    bh(12) = (est(10).*est(9)./10).^(est(12)+1);
    bh(13) = est(11).*est(9);
    bh(14) = (est(11).*est(9)./10).^(est(12)+1);
    bh(15) = est(5);
    bh(16) = (est(5)./5).^(est(12)+1);
    bh(17) = est(6);
    bh(18) = est(16);
    bh(19) = est(17);
    bh(20) = est(18);
    bh(21) = est(19);
    bh(22) = est(14);
    bh(23) = est(20);
    bh(24) = est(20).*est(21);
    bh(25) = 1-est(22);
    b_boot(i,:) = bh;
end

%% calculating quantile
for i = 1:size(b_boot,2);
    q5(i) = quantile(b_boot(:,i),0.05);
    q95(i) = quantile(b_boot(:,i),0.95);
end
result = [tr' q5' q95'];
csvwrite('confidence.csv',result)
%% histogram
for i = 1:size(b_boot,2);
    subplot(5,5,i)
    histogram(b_boot(:,i))
    y1=get(gca,'ylim');
    hold on
    plot([tr(i) tr(i)],y1)
    hold off
    title_lb = ['$b(' num2str(i) ')$'];
    title(title_lb,'interpreter','latex')
end
export_fig('hist_boot.pdf','-nofontswap','-transparent')
