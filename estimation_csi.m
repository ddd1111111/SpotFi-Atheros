clear all;
close all;
correct=0;
%%
% 參數設定
case_num=1;
derad = pi/180;         % 弧度/角度
N = 3;                 % 陣元個數
vN = 56;                % 虛擬陣元個數
M = 1;                  % 信源數目
K = 56;                % 快拍數
fc = 5.22e9;            % 傳輸頻率
c = 3e8;                % 光速
dd = 2.85e-2;              % 陣元間距
d=0:dd:ceil((N/2-1))*dd;        % 陣元間距陣列
fdd = 312.5e3;          % 子載波頻率間距
fd = 0:fdd:ceil(56/2-1)*fdd;  % 子載波頻率間距陣列
MUSIC_iter=1
SubCarrInd = [-28:-1 1:28];
% theta = [-30 0 60];	% 待估計角度
% snr = 10;             % 信噪比
for k=1:168
    csi_num=k;
    load("C:\Users\user\Desktop\20200605\save_data\recv_csi_data_5G_3M_2.mat")
    %     load("C:\Users\user\Desktop\20200605\save_data\recv_csi_data_phase_offset_12.mat")
    %     load("C:\Users\user\Desktop\20200707_output\zxc-0\point01.mat")
    %% sanitization
    for ind = 1:MUSIC_iter
        sample_csi_trace = save_csi_data{1,csi_num}(1,1:168);
        csi_plot = reshape(sample_csi_trace, K, N);
        [PhsSlope, PhsCons] = removePhsSlope(csi_plot,N,SubCarrInd,K);
        ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,N) - PhsCons*ones(K,N) ));
        csi_plot = csi_plot.*ToMult;
        relChannel_noSlope = reshape(csi_plot, K, N, M);
        sample_csi_trace_sanitized = relChannel_noSlope(:);
        sample_csi_trace_sanitized_m(ind,:) = sample_csi_trace_sanitized;
    end
    %%for read record
    csi_0=save_csi_data{1,csi_num}(1,1:56);
    csi_1=save_csi_data{1,csi_num}(1,57:112);
    csi_2=save_csi_data{1,csi_num}(1,113:168);
    %%for 42point.mat
    %     if file{csi_num, 1}.csi~=0
    %     csi_0=reshape(file{csi_num, 1}.csi(1,1,:),1,[]);
    %     csi_1=reshape(file{csi_num, 1}.csi(2,1,:),1,[]);
    %     csi_2=reshape(file{csi_num, 1}.csi(3,1,:),1,[]);
    %     else
    %         continue
    %     end
    %%for sani
    %     csi_0=sample_csi_trace_sanitized_m(1:56);
    %     csi_1=sample_csi_trace_sanitized_m(57:112);
    %     csi_2=sample_csi_trace_sanitized_m(113:168);
    %{
figure(1)
hold on
plot(1:56,phase(csi_0),'r')
plot(1:56,phase(csi_1),'g')
plot(1:56,phase(csi_2),'b')
pause(5)
close all;
    %}
    
    %% remove phase shift
    for shift=1
        phase_csi_0=phase(csi_0);
        phase_csi_0=phase_csi_0-phase(csi_1(1)).*ones(1,56)+1.1;%校正
        [csi_0_x,csi_0_y]=pol2cart(phase_csi_0,1);
        csi_0=[csi_0_x+i*csi_0_y];
        
        phase_csi_2=phase(csi_2);
        phase_csi_2=phase_csi_2-phase(csi_1(1)).*ones(1,56)-0.25;%校正
        [csi_2_x,csi_2_y]=pol2cart(phase_csi_2,1);
        csi_2=[csi_2_x+i*csi_2_y];
        
        phase_csi_1=phase(csi_1);
        phase_csi_1=phase_csi_1-phase(csi_1(1)).*ones(1,56);
        [csi_1_x,csi_1_y]=pol2cart(phase_csi_1,1);
        csi_1=[csi_1_x+i*csi_1_y];
    end
    
    figure(2)
    hold on
    plot(1:56,phase(csi_0),'r')
    plot(1:56,phase(csi_1),'g')
    plot(1:56,phase(csi_2),'b')
    %% remove pi
    if case_num==0
        if abs(phase(csi_0(1)))>abs(phase(csi_2(1)))+pi/4
            csi_0=-csi_0;
        end
        if abs(phase(csi_2(1)))>abs(phase(csi_0(1)))+pi/4
            csi_2=-csi_2;
        end
    elseif case_num==-1
        if abs(phase(csi_0(1)))>abs(phase(csi_2(1)))+pi/4
            csi_2=-csi_2;
        end
        if abs(phase(csi_2(1)))>abs(phase(csi_0(1)))+pi/4
            csi_0=-csi_0;
        end
    end
    
    figure(3)
    hold on
    plot(1:56,phase(csi_0),'r')
    plot(1:56,phase(csi_1),'g')
    plot(1:56,phase(csi_2),'b')
    
    
    csv_tmp=[csi_0 csi_1 csi_2];
    
    for p =1:4
        
        if p==1
            X1(1,1:56)=csv_tmp(1:56);
            X1(2,1:56)=csv_tmp(57:112);
            X1(3,1:56)=csv_tmp(113:168);
        end
        if p==2
            X1(1,1:56)=-csv_tmp(1:56);
            X1(2,1:56)=csv_tmp(57:112);
            X1(3,1:56)=csv_tmp(113:168);
        end
        if p==3
            X1(1,1:56)=csv_tmp(1:56);
            X1(2,1:56)=-csv_tmp(57:112);
            X1(3,1:56)=csv_tmp(113:168);
        end
        if p==4
            X1(1,1:56)=csv_tmp(1:56);
            X1(2,1:56)=csv_tmp(57:112);
            X1(3,1:56)=-csv_tmp(113:168);
        end
        
        %%
        if (phase(X1(3,1))>phase(X1(2,1)))&&(phase(X1(2,1))>phase(X1(1,1)))||...
                (phase(X1(1,1))>phase(X1(2,1)))&&(phase(X1(2,1))>phase(X1(3,1)))
            
            c1= X1(1,1:29);
            c2= X1(2,1:29);
            c3= X1(3,1:29);
            for j=2:28
                c1 = [c1;X1(1,j:28+j)];
                c2 = [c2;X1(2,j:28+j)];
                c3 = [c3;X1(3,j:28+j)];
            end
            Xe = [[c1;c2].';[c2;c3].'].';
            
            Rxx=Xe*Xe'/K;           % 計算協方差矩陣
            [EV,D]=eig(Rxx);        % 特征值分解
            EVA=diag(D)';           % 取特徵值
            [EVA,I]=sort(EVA);      % 特徵值排列
            EV=fliplr(EV(:,I));     % 依特徵值大小排列特徵向量
            En=EV(:,M+1:vN);         % 雜訊子空間
            
            % 遍歷每個角度，計算空間譜
            for iang = 1:361
                for j = -100:1:100
                    ang(iang)=(iang-181)/2;
                    tof(j+101)=j;
                    phim=derad*ang(iang);
                    a=(exp(-1i*2*pi*d*fc*sin(phim)/c).'*exp (-1i*2*pi*fd*j*1e-9)).';
                    a= reshape(a,1,numel(a)).';
                    Pmusic(iang,j+101)=abs(a'*a)/(a'*En*En'*a);
                end
            end
            Pmusic=abs(Pmusic);         % 取絕對值
            Pmusic=10*log10(Pmusic);    % 轉dB
            
            %     作圖
            %         figure(2);
            %         plot(ang,Pmusic,'Linewidth',0.5);
            %         xlabel('入射角/(degree)');
            %         ylabel('空間譜/(dB)');
            %         set(gca, 'XTick',[-90:30:90]);
            %         grid on;
            %         hold on;
            %
            
            %     figure(i);
            %     [ptof, pang] = meshgrid(tof,ang) ;
            %     mesh(pang, ptof, Pmusic);
            
            [x y]=find(Pmusic==max(max(Pmusic)));
            AoA(k,p)=(x-180)/2
            AoA_check(k)= AoA(k,p);
            if k==1
                AoA_output(k)=AoA(k,p);
            end
            if k>1
                if k>2
                    if abs(AoA(k-1,1)-AoA_output(k-2))>5&&abs(AoA(k-1,2)-AoA_output(k-2))>5&&...
                            abs(AoA(k-1,3)-AoA_output(k-2))>5&&abs(AoA(k-1,4)-AoA_output(k-2))>5
                        AoA_output(k-1)=AoA_output(k-2);
                    end
                end
                if abs(AoA_check(k)- AoA_output(k-1))<5
                    AoA_output(k)=AoA_check(k)
                    ToF=y-101
                    rssi=save_rssi_data{1,(k)};
                    rssi_amp(k)=(70-mean(rssi))/27;  %P0=70;Gamma=2.7
                    dist(k)=10^rssi_amp(k);
                    figure(4);hold on;
                    axis([0, 6, -3, 3]);
                    plot(0,0,'+');
                    plot(cos(AoA_output(k)*pi/180)*dist(k),sin(AoA_output(k)*pi/180)*dist(k),'o')
                    error_dist(k)=((cos(AoA_output(k)*pi/180)*dist(k)-2.7)^2+(sin(AoA_output(k)*pi/180)*dist(k)-1.2)^2)^(1/2);
                    %     abc(k)=phase(csi_0(1));
                    continue
                end
            end
        end
    end
end
