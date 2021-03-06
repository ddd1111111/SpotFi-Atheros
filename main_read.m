clear;
close all;
% sample CSI trace is a 90x1 vector where first 30 elements correspond to subcarriers for first rx antenna, second 30 correspond to CSI from next 30 subcarriers and so on.
% replace sample_csi_trace with CSI from Intel 5300 converted to 90x1 vector
calculate_number=0;
MUSIC_iter = 1;
% format short g;
sanitization=0;
sample_csi_traceTmp_r(1,168)=0;
particle_filter=0;
sub_carrier_number=27;

plot_map(210);

% figure(1);hold on;
% xlim([-50 50])
% ylim([-90 90])
% xlabel('ToF (ns)')
% ylabel('AoA (\circ)')
%
% figure(2);hold on;
% ground_truth_AoA=-70;
for for_5pt=7
    point_num=string(for_5pt);
    %zxc-0
    %     ground_truth=[34.9920201985587,30.9637565320735,26.5650511770780,21.8014094863518,16.6992442339936,11.3099324740202,5.71059313749966,0,-5.71059313749963,-11.3099324740202,-12.5288077091515,-6.34019174590991,0,6.34019174590992,12.5288077091515,18.4349488229220,23.9624889745782,29.0546040990771,33.6900675259798,37.8749836510982,41.1859251657097,45,49.3987053549955,54.4623222080256,60.2551187030578,56.3099324740202,50.1944289077348,45,40.6012946450045,36.8698976458440,74.0546040990771,71.5650511770780,68.1985905136482,63.4349488229220,56.3099324740202,45,26.5650511770780,0,-26.5650511770780,-45,-18.4349488229220,-14.0362434679265,-9.46232220802561,-8.13010235415598,-7.12501634890180,7.12501634890181,8.13010235415598,9.46232220802561,11.3099324740202,14.0362434679265,18.4349488229220,33.6900675259798,26.5650511770780,21.8014094863518,18.4349488229220,15.9453959009229,14.0362434679265];
    %zxc-1
    ground_truth=[-59.0362434679265,-53.1301023541560,-45,-33.6900675259798,-18.4349488229220,0,18.4349488229220,33.6900675259798,45,53.1301023541560,45,36.8698976458440,26.5650511770780,14.0362434679265,0,-14.0362434679265,-26.5650511770780,-36.8698976458440,-45,-51.3401917459099,-45,-39.8055710922652,-35.5376777919744,-32.0053832080835,-29.0546040990771,-23.9624889745782,-26.5650511770780,-29.7448812969422,-33.6900675259798,-38.6598082540901,-24.4439547804165,-19.9831065219000,-15.2551187030578,-10.3048464687660,-5.19442890773481,0,5.19442890773481,10.3048464687660,15.2551187030578,19.9831065219000,16.6992442339936,18.4349488229220,23.1985905136482,26.5650511770780,30.9637565320735,11.3099324740202,9.46232220802562,8.13010235415598,7.12501634890180,6.34019174590991,5.71059313749964,0,0,0,0,0,0];
    %user
    %     ground_truth=[5.71059313749966,0,-5.71059313749963,-11.3099324740202,-16.6992442339936,-21.8014094863518,-26.5650511770780,-30.9637565320735,-34.9920201985587,-38.6598082540901,-41.6335393365702,-37.8749836510982,-33.6900675259798,-29.0546040990771,-23.9624889745782,-18.4349488229220,-12.5288077091515,-6.34019174590991,0,6.34019174590992,7.12501634890181,8.13010235415598,9.46232220802561,11.3099324740202,14.0362434679265,0,0,0,0,0,26.5650511770780,0,-26.5650511770780,-45,-56.3099324740202,-63.4349488229220,-68.1985905136482,-71.5650511770780,-74.0546040990772,-75.9637565320735,-66.8014094863518,-60.2551187030578,-49.3987053549955,-45,-41.1859251657097,-32.0053832080835,-35.5376777919744,-39.8055710922652,-45,-51.3401917459099,-59.0362434679265,-53.1301023541560,-45,-38.6598082540901,-33.6900675259798,-29.7448812969422,-26.5650511770780];
    %     ground_truth=[-63.4 -36.8 -23.19 -18.43 0 18.43 23.96 35.53 45 54.56];
    % ground_truth=[-60 -36 0 24 45];
    ground_truth_AoA=ground_truth(str2num(point_num));
    if ground_truth(for_5pt)>0
        init=1;
    else
        init=-1;
    end
    for calculate_number=1:180
        
        
        load('C:\Users\user\Desktop\????\20210605\clear0\zxc-1\'+point_num+'.mat');
        %     if file{calculate_number, 1}.csi==0
        %         csv_tmp(1:56)=file{calculate_number-1, 1}.csi(1,1,1:56);csv_tmp(57:112)=file{calculate_number-1, 1}.csi(2,1,1:56);csv_tmp(113:168)=file{calculate_number-1, 1}.csi(3,1,1:56);
        %     else
        %         csv_tmp(1:56)=file{calculate_number, 1}.csi(1,1,1:56);csv_tmp(57:112)=file{calculate_number, 1}.csi(2,1,1:56);csv_tmp(113:168)=file{calculate_number, 1}.csi(3,1,1:56);
        %     end
        %         csv_tmp=save_csi(calculate_number,1:168);
        
        csv_tmp=[reshape(file{calculate_number*10,1}.csi(1,1,1:56),1,56) reshape(file{calculate_number*10,1}.csi(2,1,1:56),1,56) reshape(file{calculate_number*10,1}.csi(3,1,1:56),1,56)];
        rssi=file{calculate_number*10,1}.rssi;
        
        %         csv_tmp=csi_0(calculate_number,1:168);
        %         rssi=rssi_0(calculate_number,1:4);
        %         rssi=save_rssi(calculate_number,1:4);
        
        
        phase_csi_1=phase(csv_tmp(57:112));
        phase_csi_0=phase(csv_tmp(1:56));
        %zxc-1
        phase_csi_0=phase_csi_0-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56))-0.00184675654972853*[1:56]-(0.688313191228861);%-0.0103*[1:56]+0.2192;
        %         phase_csi_0=phase_csi_0-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56));%-0.0103*[1:56]+0.2192;
        %     phase_csi_0=phase_csi_0-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)+0.00184675654972853*[1:56]+(0.688313191228861);
        %         phase_csi_0=phase_csi_0-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56));
        
        [csi_0_x,csi_0_y]=pol2cart(phase_csi_0,1);%pol2cart:Polar to Cartesian Coordinates
        csi_0=[csi_0_x+i*csi_0_y];
        
        phase_csi_2=phase(csv_tmp(113:168));
        %zxc-1
        phase_csi_2=phase_csi_2-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56))+0.000208698745420363*[1:56]-(-0.285393881253289)-pi/2;%-0.3-0.0019*[1:56]+0.0246;
        %         phase_csi_2=phase_csi_2-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56))-pi/2;%-0.3-0.0019*[1:56]+0.0246;
        %     phase_csi_2=phase_csi_2-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-0.000208698745420363*[1:56]+(-0.285393881253289-pi/2);
        %         phase_csi_2=phase_csi_2-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56));
        
        [csi_2_x,csi_2_y]=pol2cart(phase_csi_2,1);
        csi_2=[csi_2_x+i*csi_2_y];
        
        phase_csi_1=phase(csv_tmp(57:112));
        phase_csi_1=phase_csi_1-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56));
        [csi_1_x,csi_1_y]=pol2cart(phase_csi_1,1);
        csi_1=[csi_1_x+i*csi_1_y];
        
        figure(2);
        subplot(1,2,1)
        xlim([1 56]);ylim([-2*pi 2*pi]);hold on;
        plot(1:56,phase(csi_0),'r')
        plot(1:56,phase(csi_1),'g')
        plot(1:56,phase(csi_2),'b')
        
        
        %     csv_tmp=[csi_0 csi_1 csi_2];
        
        %     plot(1:56,angle(csi_0),'r')
        %     plot(1:56,angle(csi_1),'g')
        %     plot(1:56,angle(csi_2),'b')
        %     tmp(calculate_number)=csi_0(1);
        %
        %%create and find up/down
        
        antenna1(1,1:56)=csi_0;
        antenna1(2,1:56)=-csi_0;
        antenna3(1,1:56)=csi_2;
        antenna3(2,1:56)=-csi_2;
        
        [up,up_num1(calculate_number)]=min(abs([mean(phase(antenna1(1,1:56)))-pi mean(phase(antenna1(2,1:56)))-pi]));
        [down,down_num1(calculate_number)]=min(abs([mean(phase(antenna1(1,1:56)))+pi mean(phase(antenna1(2,1:56)))+pi]));
        [up,up_num3(calculate_number)]=min(abs([mean(phase(antenna3(1,1:56)))-pi mean(phase(antenna3(2,1:56)))-pi]));
        [down,down_num3(calculate_number)]=min(abs([mean(phase(antenna3(1,1:56)))+pi mean(phase(antenna3(2,1:56)))+pi]));
        
        if init==1
            csv_tmp=[antenna1(up_num1(calculate_number),1:56) csi_1 antenna3(down_num3(calculate_number),1:56)];
        else
            csv_tmp=[antenna1(down_num1(calculate_number),1:56) csi_1 antenna3(up_num3(calculate_number),1:56)];
        end
        if (calculate_number==1)
            init_tmp(1)=cos(phase(csv_tmp(1)));
            init_tmp(2)=cos(phase(csv_tmp(113)));
        end
        if  abs(init_tmp(1)-cos(phase(csv_tmp(1))))>0.8
            csv_tmp(1:56)=-csv_tmp(1:56);
        end
        if  abs(init_tmp(2)-cos(phase(csv_tmp(113))))>0.8
            csv_tmp(113:168)=-csv_tmp(113:168);
        end
        
        %
        subplot(1,2,2)
        xlim([1 56]);ylim([-2*pi 2*pi]);hold on;
        plot(1:56,phase(csv_tmp(1:56)),'r')
        plot(1:56,phase(csv_tmp(57:112)),'g')
        plot(1:56,phase(csv_tmp(113:168)),'b')
        
        %     csv_tmp=[csv_tmp(113:168) csv_tmp(57:112) csv_tmp(1:56)];
        
        
        %     plot(1:56,angle(csv_tmp(1:56)),'r')
        %     plot(1:56,angle(csv_tmp(57:112)),'g')
        %     plot(1:56,angle(csv_tmp(113:168)),'b')
        %%
        
        %         csv_tmp=[csi_0 csi_1 csi_2];
        
        for num=1:MUSIC_iter
            sample_csi_traceTmp_r(num,:)=csv_tmp(1,((num*168)-167):(num*168));
        end
        
        cd("C:\Users\user\Desktop\?{???X\20210527\SpotFi-Atheros_0")
        
        %%         target_antenna = 1;
        rssi_amp(calculate_number)=(65-mean(rssi(1)))/30;  %P0=70;Gamma=2.7
        dist(calculate_number)=10^rssi_amp(calculate_number);
        %     AoA(calculate_number,1:4)=[180 180 180 180];
        for p=1:4
            
            
            if p==1
                sample_csi_traceTmp_r=[csv_tmp(1:56) csv_tmp(57:112) csv_tmp(113:168)];
            end
            if p==2
                sample_csi_traceTmp_r=[-csv_tmp(1:56) csv_tmp(57:112) csv_tmp(113:168)];
            end
            if p==3
                sample_csi_traceTmp_r=[-csv_tmp(1:56) csv_tmp(57:112) -csv_tmp(113:168)];
            end
            if p==4
                sample_csi_traceTmp_r=[csv_tmp(1:56) csv_tmp(57:112) -csv_tmp(113:168)];
            end
            
            
            fc =5220*10^6;
            M = 3;    % number of rx antennas
            %fs = 40e6; % channel bandwidth
            fs = 20e6;
            c = 3e8;  % speed of light
            % d = 6.2e-2;  % distance between adjacent antennas in the linear antenna array
            % d = 6.2e-2;
            d = 2.85e-2;
            % dTx = 2.6e-2;
            %SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; % WiFi subcarrier indices at which CSI is available
            SubCarrInd = [-28:-1 1:28]; % WiFi subcarrier indices at which CSI is available
            N = length(SubCarrInd); % number of subcarriers
            % subCarrSize = 128;  % total number fo
            fgap = 312.5e3; % frequency gap in Hz between successive subcarriers in WiFi
            lambda = c/fc;  % wavelength
            T = 1; % number of transmitter antennas
            
            % MUSIC algorithm requires estimating MUSIC spectrum in a grid. paramRange captures parameters for this grid
            % For the following example, MUSIC spectrum is caluclated for 101 ToF (Time of flight) values spaced equally between -25 ns and 25 ns. MUSIC spectrum is calculated for for 101 AoA (Angle of Arrival) values between -90 and 90 degrees.
            paramRange = struct;
            paramRange.GridPts = [101 361 1]; % number of grid points in the format [number of grid points for ToF (Time of flight), number of grid points for angle of arrival (AoA), 1]
            paramRange.delayRange = [-50 50]*1e-9; % lowest and highest values to consider for ToF grid.
            paramRange.angleRange = 90*[-1 1]; % lowest and values to consider for AoA grid.
            do_second_iter = 0;
            % paramRange.seconditerGridPts = [1 51 21 21];
            paramRange.K = floor(M/2)+1; % parameter related to smoothing.
            paramRange.L = floor(N/2); % parameter related to smoothing.
            paramRange.T = 1;
            paramRange.deltaRange = [0 0];
            
            maxRapIters = Inf;
            useNoise = 0;
            paramRange.generateAtot = 0; %2 to 0
            
            sample_csi_trace_m = sample_csi_traceTmp_r;
            for ind = 1:MUSIC_iter
                if sanitization==0
                    sample_csi_trace_sanitized_m(ind,:)=sample_csi_trace_m(ind,:);
                elseif sanitization==1
                    sample_csi_trace = sample_csi_trace_m(ind,:);
                    % ToF sanitization code (Algorithm 1 in SpotFi paper)
                    csi_plot = reshape(sample_csi_trace, N, M);
                    [PhsSlope, PhsCons] = removePhsSlope(csi_plot,M,SubCarrInd,N);
                    ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,M) - PhsCons*ones(N,M) ));
                    csi_plot = csi_plot.*ToMult;
                    relChannel_noSlope = reshape(csi_plot, N, M, T);
                    sample_csi_trace_sanitized = relChannel_noSlope(:);
                    sample_csi_trace_sanitized_m(ind,:) = sample_csi_trace_sanitized;
                end
            end
            
            % MUSIC algorithm for estimating angle of arrival
            % aoaEstimateMatrix is (nComps x 5) matrix where nComps is the number of paths in the environment. First column is ToF in ns and second column is AoA in degrees as defined in SpotFi paper
            aoaEstimateMatrix = backscatterEstimationMusic(sample_csi_trace_sanitized_m, M, N, c, fc,...
                T, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, MUSIC_iter, ones(2));
            
            
            [min_ToF,min_ToF_num]=min(aoaEstimateMatrix(:,1));
            
            AoA(calculate_number,p)=aoaEstimateMatrix(min_ToF_num,2);
            %             AoA_check(calculate_number)= AoA(calculate_number,p);
            
            %                             AoA_output(calculate_number)=AoA(calculate_number,p)
            
            
            if calculate_number==1
                %                 if (init==1&&AoA(calculate_number,p)>=0)||(init==-1&&AoA(calculate_number,p)<=0)
                if init==1&&AoA(calculate_number,p)>=0
                    AoA_output(calculate_number)=AoA(calculate_number,p);
                elseif init==-1&&AoA(calculate_number,p)<=0
                    AoA_output(calculate_number)=AoA(calculate_number,p);
                elseif init==1&&AoA(calculate_number,p)<=0
                    AoA_output(calculate_number)=asin((sin(AoA(calculate_number,p)*pi/180)*pi+pi)/pi)*180/pi;
                elseif init==-1&&AoA(calculate_number,p)>=0
                    AoA_output(calculate_number)=asin((sin(AoA(calculate_number,p)*pi/180)*pi-pi)/pi)*180/pi;
                end
                
                %%particle filter par
                if particle_filter==1
                    % point is the number of particles to be used
                    point = 400;
                    % X0, Y0 are the initial positions.
                    X0 = AoA_output(calculate_number);
                    
                    % GTData is the ground truth data of size [N,2].
                    for loop_num=1:1
                        GTData(loop_num,1) = 1;
                    end
                    % dnoise_x, dnoise_y are the standard deviations of the dynamic noise.
                    dnoise_x = 3;
                    
                    % mu_x and mu_y are vectors of means for the measurement noise.
                    mu_x = 0;
                    
                    % mnoise_x, mnoise_y are vectors of standard deviations for the measurement noise.
                    mnoise_x = 3;
                    
                    Xm = repmat([X0;0],1,point);  % Array of state variable matrices at time t
                    Xm1 = repmat([X0;0],1,point); % Array of state variable matrices at time t-1
                    Wm = ones(2, point)/point;  % Array of weights at time t
                    Wm1 = ones(2, point)/point; % Array of weights at time t
                    Index = zeros(1,point); % Index, used for resampling
                    xout(1,1) = X0;
                end
            end
            
            
            
        end
        if calculate_number>1
            [no_used,min_p]=min(abs(AoA_output(calculate_number-1)-AoA(calculate_number,:)));
            AoA_output(calculate_number)=AoA(calculate_number,min_p)
            %         [no_used,min_p]=min(abs(ground_truth_AoA-AoA(calculate_number,:)))
            %         AoA_output(calculate_number)=AoA(calculate_number,min_p)
            if particle_filter==1
                R = [mnoise_x(1,1)];
                mnoise1 = normrnd(mu_x(1),R(1,1));
                
                for j=1:point
                    d1 = normrnd(0,dnoise_x);
                    
                    %transition particles through state transition equation
                    Xm(1,j) = Xm1(1,j) + (1)*Xm1(2,j);
                    Xm(2,j) = Xm1(2,j) + d1;
                    
                    %get new measurement and update weights
                    Y(1,1) = Xm(1,j) + mnoise1;
                    
                    Wm(1,j) = Wm1(1,j)*exp(-(Y(1,1)-AoA_output(calculate_number))^2/2/R(1,1)^2);
                    
                end
                
                %normalize weights
                Wm(1,:)=Wm(1,:)/(sum(Wm(1,:)));
                
                %updating temp values for next iteration
                Xm1=Xm;
                Wm1=Wm;
                
                %desired output
                mean_p1(1,1)=0;
                for j=1:point
                    mean_p1(1,1) = mean_p1(1,1) + Xm(1,j)*Wm(1,j);
                end
                xout(calculate_number,1)=mean_p1(1,1);
                if calculate_number>4
                    AoA_output(calculate_number)=xout(calculate_number,1)
                end
                %     xout(calculate_number,2)=mean_p1(2,1);
                
                %check to see if we need to resample
                CV1=0;
                for j=1:point
                    CV1=CV1+(point*Wm(1,j)-1)^2;
                end
                CV1=CV1/point;
                
                %% For X-axis %%
                ESS1=point/(1+CV1);
                if(ESS1 < (0.5*point))
                    %resample
                    Q=cumsum(Wm1(1,:));
                    trand=rand(point+1);
                    Tsort=sort(trand);
                    Tsort(point+1)=1;
                    P=1;
                    k=1;
                    while(P<=point && k<=point)
                        if(Tsort(P) < Q(k))
                            Index(P)=k;
                            P=P+1;
                        else
                            k=k+1;
                        end
                    end
                    for P=1:point
                        Xm(1,P)=Xm1(1,Index(P));
                        Xm(2,P)=Xm1(2,Index(P));
                        Wm(1,P)=1/point;
                    end
                    Xm1(1,:)=Xm(1,:);
                    Xm1(2,:)=Xm(2,:);
                    Wm1(1,:)=Wm(1,:);
                end
                
            end
            %             figure(1);hold on;
            %             plot(cos(AoA_output(calculate_number)*pi/180)*dist(calculate_number)*100,sin(AoA_output(calculate_number)*pi/180)*dist(calculate_number)*100,'o')
        end

        
    end
    for calculate_number=1:180
        [min_error,p_num(calculate_number)]=min(abs(AoA(calculate_number,1:4)-ground_truth_AoA));
        best_AoA(calculate_number)=AoA(calculate_number,p_num(calculate_number));
    end
    
    
    %      cd("C:\Users\user\Desktop\????\20210605\output\zxc-1")
    %      save(point_num+'_removed_old.mat','AoA_output','AoA');
    
end