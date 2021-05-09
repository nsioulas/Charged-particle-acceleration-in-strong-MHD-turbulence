over1=0.001:0.01:1;
%over1=0.9;%:tM(10):tM(600);
J=0;
for over=over1
J=J+1;
% Fig. 2d
fileFormat = 'eps';    % Empty for no output
figFlag = 1;        % Figure label
nbins = 45;
binType = 'log';
lErrBars = true;

i_omit = 24;    % Omit a histogram point; set to zero to include all

tshot = over;
fileDir = 'F2o3DR10t30';
fileTag = 'F2o_t20';
%FP = load('WFP_F2o_t20.dat');
nbins=45;
binType='log';
% Fit limits in eV
fitBodyType = 'maxw';
Wf_bodyLim1 = 1.e0;
Wf_bodyLim2 = 1.e8;

Wf_tailLim3 = 1.e3;    % Tail fit
Wf_tailLim4 = 5.e4;


Wf_tailLim1 = 1.e3;    % Tail fit
Wf_tailLim2 = 5.e5;



xlims = [1.e-0 1*1.e8];
ylims = [1.e-14 5.e0];
xticks = [1.e0  1.e2  1.e4  1.e6  1.e8 1.e9];
yticks = [1.e-15 1.e-12 1.e-9 1.e-6 1.e-3 1.e0];

if exist('setHomeDir', 'file') == 2
    setHomeDir;
else
    homeDir = pwd;
end
parentDir = strcat(homeDir, '/', 'cmp4R');
%%
if exist('setFontSize', 'file') == 2
    setFontSize;
else
  
lsize = 16; % Label fontsize
    nsize = 16; % Axis fontsize

end

lexport = true;
if isempty(fileFormat), lexport = false; end
maxw = fittype('exp(-x./a) ./ sqrt(a.*pi.*x)', 'coefficients', 'a', 'independent', 'x');

if isempty(fileTag)
    fileTag = fileDir;
end
dir = strcat(parentDir, '/', fileDir);


WkinetA = Wkinet_eV;
WkinetI = W0_eV;



tM =tM;
[~,tidx] = min(abs(tM-tshot));
WfeV = Wkinet_eV(:,tidx);
WieV = W0_eV;

for j=tidx
kpM(j) =nnz(Wkinet_eV(:,j));
end

lossCoeff = kpM(:,tidx)/length(W0_eV);
%lossCoeff = 1;

[xWi,~,~,hWi,p1Wi,p2Wi] = histeb(WieV, nbins, binType);
[xWf,NWf,~,hWf,p1Wf,p2Wf] = histeb(WfeV, nbins, binType);
hWf = lossCoeff * hWf;
p1Wf = lossCoeff * p1Wf;
p2Wf = lossCoeff * p2Wf;

% Fit initial distribution
idxiReg = floor(length(xWi)/2) : length(xWi);
maxwfitInit = fit(xWi(idxiReg),hWi(idxiReg),maxw,'StartPoint',1)
aInit = coeffvalues(maxwfitInit);

% Fit final distribution
% Tail
%   Limits (indices)
[~,idxftail1] = min(abs(xWf-Wf_tailLim1));
[~,idxftail2] = min(abs(xWf-Wf_tailLim2));

[~,idxftail3] = min(abs(xWf-Wf_tailLim3));
[~,idxftail4] = min(abs(xWf-Wf_tailLim4));

%   Finetunning
idxftail1 = idxftail1;
idxftail2 = idxftail2;

idxftail3 = idxftail3;
idxftail4 = idxftail4;

idxfRegtail1 = idxftail1 : idxftail2;   % Region
idxfRegtail2 = idxftail3 : idxftail4;


plawfitTail = polyfit(log(xWf(idxfRegtail1)),log(hWf(idxfRegtail1)),1)
zfTail = plawfitTail(1);
bfTail = plawfitTail(2);

plawfitTail2 = polyfit(log(xWf(idxfRegtail2)),log(hWf(idxfRegtail2)),1)
zfTail2 = plawfitTail2(1);
bfTail2 = plawfitTail2(2);



% Body
%   Limits (indices)
[~,idxfbody1] = min(abs(xWf-Wf_bodyLim1));
[~,idxfbody2] = min(abs(xWf-Wf_bodyLim2));

idxRegbody1 = idxfbody1 : idxfbody2;  % Region

maxwFitCorrCoeff = 0.01;
if strcmp(fitBodyType, 'maxw')
    maxwfitBody = fit(xWf(idxRegbody1)/maxwFitCorrCoeff,hWf(idxRegbody1)*maxwFitCorrCoeff,maxw,'StartPoint',10000)
    aBody = coeffvalues(maxwfitBody)* maxwFitCorrCoeff;
elseif strcmp(fitBodyType, 'plaw')
    plawfitBody = polyfit(log(xWf(idxRegbody1)),log(hWf(idxRegbody1)),1)
    zBody = plawfitBody(1);
    bBody = plawfitBody(2);
end

%%
figure(1); clf; hold on; box on
set(gca, 'xscale', binType, 'yscale', binType)

x = xWf;
y = exp(-x./aInit) ./ sqrt(aInit.*pi.*x);
ph_i = plot(x,y,'co-.');
if (i_omit < 1 || i_omit > nbins)
    hWf_ = hWf;
    p1Wf_ = p1Wf;
    p2Wf_ = p2Wf;
else
    hWf_ = [1*hWf(1:i_omit-1); 0.5*(hWf(i_omit-1)+hWf(i_omit+1)); 1*hWf(i_omit+1:end)];
    p1Wf_ = [1*p1Wf(1:i_omit-1); 0.5*(p1Wf(i_omit-1)+p1Wf(i_omit+1)); 1*p1Wf(i_omit+1:end)];
    p2Wf_ = [1*p2Wf(1:i_omit-1); 0.5*(p2Wf(i_omit-1)+p2Wf(i_omit+1)); 1*p2Wf(i_omit+1:end)];
end
if lErrBars
    ph_eb = ploteb(xWf,hWf_, [], {p1Wf_,p2Wf_}, 'black', binType, 'hhy', 0.5);
    ph_f = ph_eb(1);
else
    ph_f = plot(xWf,hWf_,'bs-');
end
%{
plot(xWf(idxfbody1),hWf(idxfbody1),'bx')
plot(xWf(idxfbody2),hWf(idxfbody2),'bx')
plot(xWf(idxftail1),hWf(idxftail1),'bo')
plot(xWf(idxftail2),hWf(idxftail2),'bo')
%}
x = xWf(idxftail1-4:end);
y = exp(bfTail) * x.^(zfTail);
ph_ft = plot(x, y, '--', 'LineWidth',2, 'Color',[0 0.75 1]);

x2 = xWf(idxftail3-2:end-7);
y2 = exp(bfTail2) * x2.^(zfTail2);
%ph_ft2 = plot(x2, y2, '--', 'LineWidth',2, 'Color',[0.75 0.5 0.94]);




if strcmp(fitBodyType, 'maxw')
    x = xWf(idxfbody1:end);
    y = lossCoeff* exp(-x./aBody) ./ sqrt(aBody.*pi.*x);
    fitBodyLgnd = ['Maxwellian fit, kT = ', num2str(aBody/1000,4), ' keV'];
end
if strcmp(fitBodyType, 'plaw')
    x = xWf(idxRegbody1);
    y = exp(bBody) * x.^(zBody);
    fitBodyLgnd = ['Power-law fit, z = ' num2str(-zfTail,2)];
end
if ~isempty(fitBodyType)
    ph_fb = plot(x, y, 'ro-.','LineWidth',.3);
end

%ph_FP = plot(FP(:,1), FP(:,2), 'k');
hold on
axis([xlims ylims])
if ~isempty(xticks)
    set(gca,'XTick',xticks)
end
if ~isempty(yticks)
    set(gca,'YTick',yticks)
end
ylabel('p(W)', 'Fontsize', lsize)
xlabel('W_{kin} [eV]', 'Fontsize', lsize)

if isempty(fitBodyType)
    legend([ph_i, ph_f,ph_ft2, ph_ft], 'Location', 'Best',...
           {'Initial distribution',...
            ['Distribution after t = ',num2str(tM(tidx),2), ' sec'],...
            ['Power-law fit, k = ' num2str(-zfTail2,2)],...
            ['Power-law fit, k = ' num2str(-zfTail,2)],...
            'Fokker-Planck solution'...
           },'FontSize',10)
else
    legend([ph_i, ph_f, ph_fb, ph_ft], 'Location', 'Best',...
           {'Initial distribution',...
            ['Distribution after t = ',num2str(tM(tidx)), ' sec'],...
            fitBodyLgnd,....
            ['Power-law fit, z = ' num2str(-zfTail,2)],...
            'Fokker-Planck solution'...
           },'FontSize',10)
end

set(gca, 'Fontsize', nsize)

if lexport
    % Save fig
    set(gcf,'windowstyle','normal');
    fname = ['nWfitFP_' fileTag '.' fileFormat];
    hgexport(gcf, fname, hgexport('factorystyle'), 'Format', fileFormat);
    %set(gcf,'windowstyle','docked');
end
indexa1(J)=-zfTail;
%tm1((J))=over1;
end
%{hold off
figure (20)
plot(over1,indexa1)
%}
  set(gca, 'Fontsize', nsize)
set(gcf,'paperpositionmode','auto');
    set(gcf,'windowstyle','normal');
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
%set(gca,'fontweight','bold')


opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';

if imodel == 1

    temp=['png_files_upper\f8a.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f8a.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f8a.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f8a.png']; saveas(gca,temp);
    temp1=['new_eps\f8a.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f8a.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f8a P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f8a P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f8a P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
end

clear WkinetA WkinetI