clear; close all;
%%%%%%%% MUSIC for Uniform Linear Array%%%%%%%%
% 參數設定
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
sub_carrier_number=1;
% theta = [-30 0 60];	% 待估計角度
% snr = 10;             % 信噪比

% 構建信號模型
% A=exp(-1i*2*pi*d.'*sin(theta*derad));
% S=randn(M,K); X=A*S;
% X1=awgn(X,snr,'measured');

% MUSIC 演算法
% load('60.mat');     % 導入數據
% point_num=1;
% point_num=string(point_num);
% load('C:\Users\user\Desktop\資料\20210412\one_pack\point'+point_num+'_B.mat')
for for_5pt=1:19
    point_num=string(for_5pt);
    load('C:\Users\user\Desktop\資料\20210605\clear0\zxc-1\'+point_num+'.mat');
    % % load("C:\Users\user\Desktop\程式碼\test_AoA.mat")
    %
    % csv_tmp=save_csi(1,1:168);
    
    for calculate_number=1:180
        csv_tmp=[reshape(file{calculate_number*10,1}.csi(1,1,1:56),1,56) reshape(file{calculate_number*10,1}.csi(2,1,1:56),1,56) reshape(file{calculate_number*10,1}.csi(3,1,1:56),1,56)];
        
        phase_csi_1=phase(csv_tmp(57:112));
        abs_csi_1=abs(csv_tmp(57:112));
        phase_csi_0=phase(csv_tmp(1:56));
        abs_csi_0=abs(csv_tmp(1:56));
        %     phase_csi_0=phase_csi_0-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56))-0.00184675654972853*[1:56]-(0.688313191228861);
        phase_csi_0=phase_csi_0-0.00184675654972853*[1:56]-(0.688313191228861);
        
        [csi_0_x,csi_0_y]=pol2cart(phase_csi_0,1);%pol2cart:Polar to Cartesian Coordinates
        csi_0=[csi_0_x+i*csi_0_y];
        
        phase_csi_2=phase(csv_tmp(113:168));
        abs_csi_2=abs(csv_tmp(113:168));
        %     phase_csi_2=phase_csi_2-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56))+0.000208698745420363*[1:56]-(-0.285393881253289)-pi/2;
        phase_csi_2=phase_csi_2+0.000208698745420363*[1:56]-(-0.285393881253289)-pi/2;
        
        [csi_2_x,csi_2_y]=pol2cart(phase_csi_2,1);
        csi_2=[csi_2_x+i*csi_2_y];
        
        phase_csi_1=phase(csv_tmp(57:112));
        %     phase_csi_1=phase_csi_1-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56));
        phase_csi_1=phase_csi_1;
        [csi_1_x,csi_1_y]=pol2cart(phase_csi_1,1);
        csi_1=[csi_1_x+i*csi_1_y];
        %
        csv_tmp=[csi_0 csi_1 csi_2];
        % csv_tmp=[0.997273367837507 - 0.0737958657516527i,0.995419008687697 - 0.0956085620810322i,0.999207134247950 + 0.0398133503739587i,0.998395817083084 + 0.0566197176873941i,0.998119154128017 + 0.0613037858763313i,0.999850250078509 - 0.0173054158558167i,0.999999888723352 - 0.000471755534024664i,0.997704454557446 + 0.0677186928124599i,0.999130377235163 + 0.0416951949980071i,0.999891333215858 + 0.0147418370569938i,0.998097226759095 - 0.0616597594530207i,0.980837876178302 - 0.194825718666806i,0.988744951788655 - 0.149610896369381i,0.990793953514265 - 0.135378512621365i,0.947311304157903 - 0.320314365919877i,0.913768908234919 - 0.406234393353351i,0.876356053475711 - 0.481663853259176i,0.892773534424942 - 0.450505733848523i,0.890464896606738 - 0.455051939794955i,0.834005519963431 - 0.551756098897445i,0.715393546221157 - 0.698721742917105i,0.710246239700586 - 0.703953321599648i,0.669719800167155 - 0.742613889759723i,0.692739611606100 - 0.721187791432876i,0.753868122327844 - 0.657025763678938i,0.666623035375490 - 0.745395015214596i,0.618337752289934 - 0.785912478646975i,0.684626257203210 - 0.728894291339920i,0.697062388849860 - 0.717010478341096i,0.625698337674618 - 0.780065119224811i,0.622736550299963 - 0.782431587373939i,0.707220471650081 - 0.706993072440627i,0.606528009629157 - 0.795062119293388i,0.574885297964242 - 0.818234009427966i,0.565720468220308 - 0.824597084542867i,0.583923316792966 - 0.811808819923448i,0.653631601332798 - 0.756812876303728i,0.619366866270813 - 0.785101703581054i,0.649851612507234 - 0.760061103939510i,0.581732809231433 - 0.813379947296284i,0.561581964366215 - 0.827421112432226i,0.557096710047593 - 0.830447623666989i,0.593917080411823 - 0.804526259108487i,0.594442606017490 - 0.804138040482563i,0.529617420737334 - 0.848236634230999i,0.628108390384776 - 0.778125857383139i,0.626336838413892 - 0.779552541427254i,0.665146069723750 - 0.746713268886424i,0.636010353116189 - 0.771680523746077i,0.567133566926845 - 0.823625835719615i,0.603261882093246 - 0.797543165987469i,0.524284157926845 - 0.851543376315581i,0.617705769152822 - 0.786409297220805i,0.634379285743297 - 0.773021941357310i,0.522534497711706 - 0.852618143544444i,0.394525387536636 - 0.918885041008432i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,1.00000000000000 + 0.00000000000000i,0.741148272382810 + 0.671341372435795i,0.769675486134805 + 0.638435310774045i,0.674304570090268 + 0.738453347717633i,0.637353747422217 + 0.770571346889343i,0.648279970230963 + 0.761402048984203i,0.707721845257336 + 0.706491181647408i,0.667549479363279 + 0.744565438763991i,0.494448628283043 + 0.869206853395104i,0.506298776891931 + 0.862358132400765i,0.447199733470503 + 0.894434121880371i,0.558320220121461 + 0.829625537097023i,0.729160007855217 + 0.684343249359984i,0.549173501979866 + 0.835708361046586i,0.399629458464006 + 0.916676767420101i,0.631651957090773 + 0.775252091324749i,0.682779946061502 + 0.730624079302245i,0.709905181060760 + 0.704297262456052i,0.551691398103018 + 0.834048320697990i,0.459318916117480 + 0.888271429967588i,0.570625529009187 + 0.821210390608269i,0.797177754630141 + 0.603744670802854i,0.780057815496857 + 0.625707443205106i,0.805770822608184 + 0.592227474399264i,0.757725020327961 + 0.652573975552957i,0.556601372378438 + 0.830779701404915i,0.753673383233246 + 0.657249139524543i,0.840241802500974 + 0.542211871255061i,0.691594919893771 + 0.722285585331127i,0.738062997295981 + 0.674731807477959i,0.823917664124338 + 0.566709522369171i,0.853269247466572 + 0.521470604471461i,0.658988973601577 + 0.752152599325124i,0.862747984301294 + 0.505634171693384i,0.879409590772229 + 0.476065932049145i,0.886103720702353 + 0.463486996751200i,0.842198659579143 + 0.539167337477980i,0.696189782140900 + 0.717857776472893i,0.758122602321694 + 0.652112045471469i,0.675043014496313 + 0.737778373618887i,0.803616287075585 + 0.595147765808501i,0.823579836512508 + 0.567200363972054i,0.835714190821991 + 0.549164630378492i,0.787237803726753 + 0.616649528000694i,0.775679502002334 + 0.631127015879855i,0.895120492374454 + 0.445824297376573i,0.784652626181401 + 0.619935687169751i,0.810653484512542 + 0.585526197575885i,0.781292441266995 + 0.624165139381445i,0.874641787484776 + 0.484769784109361i,0.952805402242283 + 0.303581727806405i,0.959126129659173 + 0.282978916891374i,0.995626482910054 + 0.0934232654543603i,0.977164521137070 + 0.212484113831979i,0.982298819489578 + 0.187320658843017i,0.999999932376926 - 0.000367758267488303i,0.981225020559265 - 0.192866427945537i];
        
        for p = 1:4
            %     X1 = reshape(original_csi{i}.csi(:,1,1:56),3,56);
            %     c1= X1(1,1:29);
            %     c2= X1(2,1:29);
            %     c3= X1(3,1:29);
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
            
            figure(1);
            hold on
            plot(phase(X1(1,:)),'r','Linewidth',2)
            plot(phase(X1(2,:)),'g','Linewidth',2)
            plot(phase(X1(3,:)),'b','Linewidth',2)
            xlim([1 56]);ylim([-pi pi]);
            
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
                for j = -100:100
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
            figure(2);
            [ptof, pang] = meshgrid(tof,ang) ;
            mesh(pang, ptof, Pmusic);
            xlabel('AoA');ylabel('ToF');zlabel('Pmusic(dB)')
            
            [x y]=find(Pmusic==max(max(Pmusic)));
            AoA(calculate_number,p)=(x-181)/2
            ToF=y-101
            z(calculate_number,p)=max(max(Pmusic));
            
        end
    end
    cd("C:\Users\user\Desktop\資料\20210605\output\zxc-1")
    save(point_num+'_find_pmusic.mat','AoA','z');
end
