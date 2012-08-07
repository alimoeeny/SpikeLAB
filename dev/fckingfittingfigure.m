load /Users/ali/Desktop/matlab.735066.5651.mat

datadelta = PopPSTH(1,:) - PopPSTH(2,:);


PsychPerformanceAtEnd = 75;
lag = 0;

x = [500 + lag 2000];
y = [0 100-PsychPerformanceAtEnd];

%% ============= linear fit

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Upper = [Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel( 'x' );
ylabel( 'y' );
grid on


p1 = fitresult.p1;
p2 = fitresult.p2;
%Pflip(x) = p1*x + p2;


%% 
Trange = [500:2000];
h = figure(1357); clf, hold on,
ft = [];
leastsqrt = [];
for fri = [1:1:50]
    for t = Trange
        ft(fri, t) =  fri - (2 * fri * ((p1 *t + p2) / 100));
    end
    leastsqrt(fri) = sqrt(sum(abs(datadelta([1000:2000]) - ft(fri,[1000:2000]))));
    plot(Trange, ft(fri, Trange), 'Color', [1/fri 0 fri/50]);
end

figure (1380); clf, hold on
plot(leastsqrt);

%%

PlotRange = [1:2000];

figure(1999), clf, hold on,
plot(PopPSTH(1,PlotRange));
plot(PopPSTH(2,PlotRange));
plot(datadelta(PlotRange), 'LineWidth',3, 'Color', [1 0 0]);
plot([501:2000], ft(find(leastsqrt==min(leastsqrt)),[501:2000]), 'LineWidth',3, 'Color', [0 0 1]);
for i = 1:size(ft, 1)
    plot([501:2000],ft(i,[501:2000]));
end
% plot(ft(1,PlotRange));
% plot(ft(end,PlotRange));

reflinexy(500,100);
reflinexy(1000,100);