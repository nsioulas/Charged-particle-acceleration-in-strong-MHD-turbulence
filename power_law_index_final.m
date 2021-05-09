

lsize = 16; % Label fontsize
    nsize = 16; % Axis fontsize
    


J=0;

start=27;    rd1start=1;      end2=45; %number of bins
ending=end2-1;     rd1end=end2; % nP=length(rmsd);



a12=[600  720 750  790 800 830 850  ];%1*(ntM/100);
%a12=1:99:800;
%a12=700:10:860;
%a=zeros(1,length(a12));
%a12=310:5:400;
index=zeros(1,length(a12));
a1=zeros(1,length(a12));


%% Gia Maxwell-boltzman katanomh
%W0_eV = W0.*erg2eV;

W0_eV1=W0_eV;
 
  %rd1start=1;      end2=45; %number of bins
  %rd1end=end2;  nP=length(rmsd);
tacc=2;

split2=logspace(0,log10(max(W0_eV)),end2);

W_esc_half1=zeros(1,end2);
p2=zeros(1,end2);


for jN=1
 
  W0_eV=W0_eV1;
   W_esc_half1(jN)=(1/2)*(split2(jN));
for i=1:nP
        if W0_eV(1,i)>=split2(jN)

           W0_eV(1,i)=0;

       end
end

        p2(1,jN)=nnz(W0_eV(:));

end

W0_eV=W0_eV1;
 
for jN=2:end2
    disp(['percentage: ' num2str(jN/end2)])
   %w_tacc=w_tacc1;
    W_esc_half1(jN)=(1/2)*(split2(jN)-split2(jN-1))+split2(jN-1);
for i=1:nP
       
        if W0_eV(1,i)<=split2(jN-1) || W0_eV(1,i)>=split2(jN)

        W0_eV(1,i)=0;

        end
       
end
       
    p2(1,jN)=nnz(W0_eV(:));   
    W0_eV=W0_eV1;
end

total=sum(p2);
relative1=zeros(1,end2);
for jN=1:end2
     if jN==1
         relative1(1,jN)=p2(1,jN)/(total*(split2(jN)));
     elseif jN>=end2-1
            relative1(1,jN)=p2(1,jN)/(total*(split2(jN)-split2(jN-1))); 
     else
    relative1(1,jN)=p2(1,jN)/(total*(split2(jN)-split2(jN-1)));
     end
end








for i=a12

 J=J+1;
 
tacc=i;
 %if i>=200 && i<=250
   % ending=ending-1; 
 %end
keV=0;   % 1: x axis in keV 
         % 0: x axis in eV

index(i)=i;


%split2=linspace((min(Wkinet_eV(:,tacc))),(max(Wkinet_eV(:,tacc))),end2);
split2=logspace(0,log10(max(Wkinet_eV(:,tacc))),end2);
%split2=logspace(log10(min(Wkinet_eV(:,tacc))),log10(max(Wkinet_eV(:,tacc))),end2);
W_esc_half=zeros(1,end2);
p1=zeros(1,end2);


w_tacc=zeros(nP,1);
w_tacc1=zeros(nP,1);
for i=1:nP
    
w_tacc(i,1)= Wkinet_eV(i,tacc);
w_tacc1(i,1)= Wkinet_eV(i,tacc);
end

for jN=1
 
  w_tacc=w_tacc1;
   W_esc_half(jN)=(1/2)*(split2(jN));
for ip=1:nP
        if w_tacc(ip,1)>=split2(jN)

        
w_tacc(ip,1)=0;

       end
end

        p1(1,jN)=nnz(w_tacc(:));
end

w_tacc=w_tacc1;
 
for jN=2:end2
    disp(['percentage: ' num2str(jN/end2)])
   w_tacc=w_tacc1;
    W_esc_half(jN)=(1/2)*(split2(jN)-split2(jN-1))+split2(jN-1);
for ip=1:nP
       
        if w_tacc(ip,1)<=split2(jN-1) || w_tacc(ip,1)>=split2(jN)

w_tacc(ip,1)=0;

        end
       
end
       
p1(1,jN)=nnz(w_tacc(:));
        
    w_tacc=w_tacc1;
end

total=sum(p1);
relative=zeros(1,end2);
for jN=1:end2
     if jN==1
         relative(1,jN)=p1(1,jN)/(total*(split2(jN)));
     else
    relative(1,jN)=p1(1,jN)/(total*(split2(jN)-split2(jN-1)));
    
     end
end



if keV
    W_esc_half=W_esc_half/10^3;
   % W_esc_half1=W_esc_half1/10^3;
end


figure(1)
fitResults1 = polyfit(log10(W_esc_half(start:ending)),log10(relative(start:ending)),1);
pol=polyval(fitResults1,log10(W_esc_half(start:ending)));
a=fitResults1(1);
b=fitResults1(2);
polyfit_str = ['<(?r)^2> ~ ? ^{'  (num2str(a)) ' } '] % + ' num2str(b)]  polyfit_str will be : y = 4*x + 2
a1(J)=a;
fit1=W_esc_half(start:ending).^a*10^(b);

   %loglog((W_esc_half1(rd1start:rd1end)),(relative1(rd1start:rd1end)),'red-','LineWidth',2)

   if J==1
        loglog((W_esc_half1(rd1start:rd1end)),(relative1(rd1start:rd1end)),'red-','LineWidth',2)
        
   end
   hold on
   loglog((W_esc_half(rd1start:rd1end)),(relative(rd1start:rd1end)),'LineWidth',2)
   
   %loglog((W_esc_half(start:(length(pol)+start-1))),fit1,'cyan--','LineWidth',2);
   

end
 hold off
for u=1:length(a12)
       str = {strcat('kT=100eV')} ; % at the end of first loop, z being loop output;
  for i=1:length(a12)
      str = [str , strcat('t=' , num2str(tM(a12(i)),1),' s')]%, '[sec]'),strcat(' Power-law fit, k:',num2str(a1(i)))]; % after 2nd loop;
 % str = [str , strcat('Distr.. after t=' , num2str(tM(a12(i))), '[sec]')]; % after 2nd loop;
  end
     
end


legend(str(:),'Location', 'SouthWest','FontSize',12);
xlabel('w_{kinet}  [eV]','FontSize',lsize)
ylabel('p(w_{kinet})','FontSize',lsize)

xlim([10^(0),10^9])
ylim([1.e-13,10^(0)])
XTick = [ 10^(0)   10^(2)  10^(4)  10^(6)  10^(8) 10^(9)];
set(gca,'xtick',XTick)

YTick = [10^(-12) 10^(-10) 10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(0)];
set(gca,'ytick',YTick)


  set(gca, 'Fontsize', nsize)
set(gcf,'paperpositionmode','auto');
    set(gcf,'windowstyle','normal');
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')


opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';

if imodel == 1

    temp=['png_files_upper\f17.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f17.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f17.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f17.png']; saveas(gca,temp);
    temp1=['new_eps\f17.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f17.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f17 P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f17 P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f17 P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
 end


%{
figure (2)
index12=a1<0;
plot(tM(a12),abs(a1(index12))-0.15,'b-','LineWidth',2)
set(gca,'fontsize',17)
set(gcf,'paperpositionmode','auto');
box on
xlim([0 2.75])
  
  xlabel('Time [sec]','FontSize',lsize)
ylabel('Power-law index k','FontSize',lsize)

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

if upper_step_limit
    
    
     if imodel == 1
    temp=['png_files_upper\f8b.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f8b.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f8b.fig']; saveas(gca,temp2); 

     else
         temp=['png_files_upper\index_of_Wkinet-Time_noEnerg, Lmin=',num2str(Dr_min),', MLim=',num2str(Dr_max),', tfin=',num2str(tfin),', nP=',num2str(nP),'.png']; saveas(gca,temp);
    temp1=['eps_files_upper\index_of_Wkinet-Time_noEnerg, Lmin=',num2str(Dr_min), ', MLim=',num2str(Dr_max),', tfin=',num2str(tfin),', nP=',num2str(nP),'.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\index_of_Wkinet-Time_noEnerg, Lmin=',num2str(Dr_min),', MLim=',num2str(Dr_max),', tfin=',num2str(tfin),', nP=',num2str(nP),'.fig']; saveas(gca,temp2); 

         
     end
   
     if imodel == 1
 temp3=['png_files\index_of_Wkinet-Time, Lmin=',num2str(Dr_min),', tfin=',num2str(tfin),', nP=',num2str(nP),'.png']; saveas(gca,temp3);
    temp4=['eps_files\index_of_Wkinet-Time, Lmin=',num2str(Dr_min),', tfin=',num2str(tfin),', nP=',num2str(nP),'.eps']; saveas(gca,temp4); 
    temp5=['fig_files\index_of_Wkinet-Time, Lmin=',num2str(Dr_min),', tfin=',num2str(tfin),', nP=',num2str(nP),'.fig']; saveas(gca,temp5); 
    
     else
          temp3=['png_files\index_of_Wkinet-Time_noEnerg, Lmin=',num2str(Dr_min),', tfin=',num2str(tfin),', nP=',num2str(nP),'.png']; saveas(gca,temp3);
    temp4=['eps_files\index_of_Wkinet-Time_noEnerg, Lmin=',num2str(Dr_min),', tfin=',num2str(tfin),', nP=',num2str(nP),'.eps']; saveas(gca,temp4); 
    temp5=['fig_files\index_of_Wkinet-Time_noEnerg, Lmin=',num2str(Dr_min),', tfin=',num2str(tfin),', nP=',num2str(nP),'.fig']; saveas(gca,temp5); 
    
         
     end
end

%}