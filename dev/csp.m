clear, clc

mycolors = 'rgbcmyk'; %= {[1,0,0],[0,1,1],[0,0,1],[1,1,0]};

%% load data
%load ~/bgc/data/dae/dae437/dae437fullv.mat
%load ~/bgc/data/dae/dae438/dae438fullv.mat
%load ~/bgc/data/dae/dae439/dae439fullv.mat
load /bgc/data/dae/dae440/dae440fullv.mat

% V{1} = res.V(1, 1:1000000);
% V{2} = res.V(2, 1:1000000);
% V{3} = res.V(3, 1:1000000);
% V{4} = res.V(4, 1:1000000);
% for c = 1:1
%     [e, t] = Spikes(V{c}, -2);
%     spk{c} = {e , t};
% end

%%
Vs = res.V(:,1:1000000);

%%
[e, t] = SpikesMulti(Vs, 2, -2);
    

%%

figure(111), 
for i = 1:size(e,1), 
    plot(squeeze(e(i,:,:))'), 
    pause(0.15), 
end

%%

figure(123), clf
for i = 1:size(e,2), 
    subplot(4,1,i)
    plot(squeeze(e(:,i,:))', mycolors(i));
end

%%

figure(124), clf
for i = 1:size(e,2), 
    subplot(4,1,i)
    plot(mean(squeeze(e(:,i,:))), mycolors(i));
end


%%
for i  = 1: size(e,1)
    b = [];
    for j = 1: size(e,2)
        b = [b, squeeze(e(i,j,:))']; 
    end
    ee(i, :) = squeeze(b);
end

[pc, zscores, pcvars] = princomp(ee');
%%
figure(210), scatter(std(ee(:,:)'), pc(:,1), 'filled')
figure(211), scatter(std(ee(:,1:51)'), pc(:,1), 'filled')
figure(212), scatter(std(ee(:,52:102)'), pc(:,1), 'filled')
figure(213), scatter(std(ee(:,103:153)'), pc(:,1), 'filled')
figure(214), scatter(std(ee(:,154:204)'), pc(:,1), 'filled')


%%
for a = 1:4
    for b = 1: 4
        for i = 1: size(e,1), 
            corrs(a, b, i) = corr(squeeze(e(i, a, :)), squeeze(e(i, b, :))); 
        end
    end
end

figure, imagesc(mean(corrs,3));

% 
% 
% % generate filters
% cd ~/Desktop/ica/
% LPf = LPFliterGen;
% HPf = HPFliterGen;
% 
% % load filter coefficientes
% %load ~/Desktop/ica/filers.mat
% 
% % lowpass filter
% V_lp = filter(, 1, res.V);
% 
% 
% % highpass filter
% %V_hp = filter(NumHP, 1, res.V);
% V_lp = filter(NumLP, 1, res.V);
% 
% %V_hp_lp = filter(NumLP, 1, V_hp);
% V_lp_hp = filter(NumLP, 1, V_lp);
% 
% figure(13), clf
% plot(res.V(10055000:10070000))
% %hold on, plot(V_hp_lp(10055000:10070000), 'm')
% hold on, plot(V_lp_hp(10055000:10070000), 'm')
% 
% 
% % ica | pca


