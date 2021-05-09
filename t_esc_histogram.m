 lsize = 16; % Label fontsize
    nsize = 16; % Axis fontsize

start=20;    rd1start=1;      end2=30; %number of bins
ending=end2;     rd1end=end2;  %nP=length(rmsd);
tacc=2;


t_esc1=t_esc(t_esc>0);
t_esc1_=t_esc1;
nP=length(t_esc1);
%split2=linspace(log10(min(t_esc1)),log10(max(t_esc1)),end2);
%split2=logspace(0,log10(max(t_esc1)),end2);
split2=logspace(log10(min(t_esc1)),log10(max(t_esc1)),end2);
t_esc_half=zeros(1,end2);
p1=zeros(1,end2);






for jN=1
 
  t_esc1=t_esc1_;
   t_esc_half(jN)=(1/2)*(split2(jN));
for i=1:nP
        if t_esc1(i)>=split2(jN)

        
t_esc1(i)=0;

       end
end

    
        p1(1,jN)=nnz(t_esc1(:));
 
     
    
   % dw1=dw1_;
 %t_esc1_=t_esc1;
end

t_esc1=t_esc1_;
 
for jN=2:end2
    disp(['percentage: ' num2str(jN/end2)])
   t_esc1=t_esc1_;
    t_esc_half(jN)=(1/2)*(split2(jN)-split2(jN-1))+split2(jN-1);
for i=1:nP
       
        if t_esc1(i)<=split2(jN-1) || t_esc1(i)>=split2(jN)

       
t_esc1(i)=0;

        end
       
end
       
p1(1,jN)=nnz(t_esc1(:));
        
    t_esc1=t_esc1_;
end

total=sum(p1);
relative=zeros(1,end2);
for jN=1:end2
     if jN==1
         relative(1,jN)=p1(1,jN)/(total*(split2(jN)));
         %relative(1,jN)=p1(1,jN)/total;
     else
   relative(1,jN)=p1(1,jN)/(total*(split2(jN)-split2(jN-1)));
    %relative(1,jN)=p1(1,jN)/total;
     end
end



figure(2)
fitResults1 = polyfit(log10(t_esc_half(start:ending)),log10(relative(start:ending)),1);
pol=polyval(fitResults1,log10(t_esc_half(start:ending)));
a=fitResults1(1);
b=fitResults1(2);
polyfit_str = ['<(?r)^2> ~ ? ^{'  (num2str(a)) ' } '] % + ' num2str(b)]  polyfit_str will be : y = 4*x + 2


fit1=t_esc_half(start:ending).^a*10^(b);


loglog((t_esc_half(rd1start:rd1end)),(relative(rd1start:rd1end)),'blue-','LineWidth',2)
hold on
loglog((t_esc_half(start:(length(pol)+start-1))),fit1,'cyan--','LineWidth',2)
legend({ strcat('Distribution of t_{esc} '),strcat(' Power-law fit, k:',num2str(a))},'Location','SouthWest','FontSize',12)%,'fit for [0-20 sec]'})
xlabel('t_{esc}  [sec]', 'Fontsize', lsize)
ylabel('p(t_{esc})', 'Fontsize', lsize)

xlim([10^(-4),10^3])
ylim([10^(-7),10^0])

set(gca, 'Fontsize', nsize)
set(gcf,'paperpositionmode','auto');
    set(gcf,'windowstyle','normal');
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','bold')

opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';

hold off



if imodel == 1

    temp=['png_files_upper\f5b.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f5b.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f5b.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f5b.png']; saveas(gca,temp);
    temp1=['new_eps\f5b.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f5b.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f5b P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f5b P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f5b P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
 end
