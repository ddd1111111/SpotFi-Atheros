clear;
close all;
init_parameter;
plot_map(210);
% figure(1);hold on;
% xlim([-50 50])
% ylim([-90 90])
% xlabel('ToF (ns)')
% ylabel('AoA (\circ)')
%
% figure(2);hold on;
while(1)
    pause(0.05)
    clk_time=clock;
    if clk_time(6)-fix(clk_time(6))<0.1
        calculate_number=calculate_number+1;
        cd("..\python_0")
        csv_tmp=csvread('csi_tmp.csv');
        rssi=csvread('rssi.csv');
        
        phase_csi_0=phase(csv_tmp(1:56));
        phase_csi_0=phase_csi_0-phase(csv_tmp(57)).*ones(1,56);%校正
        [csi_0_x,csi_0_y]=pol2cart(phase_csi_0,1);
        csi_0=[csi_0_x+i*csi_0_y];
        
        phase_csi_2=phase(csv_tmp(113:168));
        phase_csi_2=phase_csi_2-phase(csv_tmp(57)).*ones(1,56);%校正
        [csi_2_x,csi_2_y]=pol2cart(phase_csi_2,1);
        csi_2=[csi_2_x+i*csi_2_y];
        
        phase_csi_1=phase(csv_tmp(57:112));
        phase_csi_1=phase_csi_1-phase(csv_tmp(57)).*ones(1,56);
        [csi_1_x,csi_1_y]=pol2cart(phase_csi_1,1);
        csi_1=[csi_1_x+i*csi_1_y];
        
        
%         phase_csi_1=phase(csv_tmp(57:112));
%         phase_csi_0=phase(csv_tmp(1:56));
%         phase_csi_0=phase_csi_0-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56))-0.00184675654972853*[1:56]-(0.688313191228861);%-pi/2;%-0.0103*[1:56]+0.2192;
%         
%         [csi_0_x,csi_0_y]=pol2cart(phase_csi_0,1);%pol2cart:Polar to Cartesian Coordinates
%         csi_0=[csi_0_x+i*csi_0_y];
%         
%         phase_csi_2=phase(csv_tmp(113:168));
%         phase_csi_2=phase_csi_2-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56))+0.000208698745420363*[1:56]-(-0.285393881253289-pi/2);%-0.3-0.0019*[1:56]+0.0246;
%         
%         [csi_2_x,csi_2_y]=pol2cart(phase_csi_2,1);
%         csi_2=[csi_2_x+i*csi_2_y];
%         
%         phase_csi_1=phase(csv_tmp(57:112));
%         phase_csi_1=phase_csi_1-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56));
%         [csi_1_x,csi_1_y]=pol2cart(phase_csi_1,1);
%         csi_1=[csi_1_x+i*csi_1_y];
        
        figure(2)
        hold on
        axis([0, 56, -4, 4]);
        plot(1:56,phase(csi_0),'r')
        plot(1:56,phase(csi_1),'g')
        plot(1:56,phase(csi_2),'b')
        
        csv_tmp=[csi_0 csi_1 csi_2];
        
        save_csi(calculate_number,1:168)=csv_tmp;
        save_rssi(calculate_number,1:4)=rssi;
        
        for num=1:MUSIC_iter
            sample_csi_traceTmp_r(num,:)=csv_tmp(1,((num*168)-167):(num*168));
        end
        
        cd("..\SpotFi-Atheros_0")
        
        AoA(calculate_number,1:4)=[180 180 180 180];
        for p=1:2
            
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
            
            if ((phase(sample_csi_traceTmp_r(113))>=phase(sample_csi_traceTmp_r(57)))&&(phase(sample_csi_traceTmp_r(57))>=phase(sample_csi_traceTmp_r(1)))...
                    ||((phase(sample_csi_traceTmp_r(1))>=phase(sample_csi_traceTmp_r(57)))&&(phase(sample_csi_traceTmp_r(57))>=phase(sample_csi_traceTmp_r(113)))))
                
                %                 fc =5220*10^6;
                fc =2437*10^6;
                M = 3;    % number of rx antennas
                %fs = 40e6; % channel bandwidth
                fs = 20e6;
                c = 3e8;  % speed of light
                % d = 6.2e-2;  % distance between adjacent antennas in the linear antenna array
                d = 6.2e-2;
                %                 d = 2.85e-2;
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
                %                 AoA_check(calculate_number)= AoA(calculate_number,p);
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
                    for particle=particle_filter
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
                            dnoise_x = d_noise;
                            
                            % mu_x and mu_y are vectors of means for the measurement noise.
                            mu_x = measurement_noise;
                            
                            % mnoise_x, mnoise_y are vectors of standard deviations for the measurement noise.
                            mnoise_x = m_noise;
                            
                            Xm = repmat([X0;0],1,point);  % Array of state variable matrices at time t
                            Xm1 = repmat([X0;0],1,point); % Array of state variable matrices at time t-1
                            Wm = ones(2, point)/point;  % Array of weights at time t
                            Wm1 = ones(2, point)/point; % Array of weights at time t
                            Index = zeros(1,point); % Index, used for resampling
                            xout(1,1) = X0;
                        end
                    end
                    %                 end % if (init==1&&AoA(calculate_number,p)>=0)||(init==-1&&AoA(calculate_number,p)<=0)
                end
                if calculate_number>1
                    if AoA(calculate_number,p)<=0
                        AoA(calculate_number,p+2)=asin((sin(aoaEstimateMatrix(min_ToF_num,2)*pi/180)*pi+pi)/pi)*180/pi;
                    end
                    if AoA(calculate_number,p)>=0
                        AoA(calculate_number,p+2)=asin((sin(aoaEstimateMatrix(min_ToF_num,2)*pi/180)*pi-pi)/pi)*180/pi;
                    end
                end
            end
        end
        if calculate_number>1
            [no_used,min_p]=min(abs(AoA_output(calculate_number-1)-AoA(calculate_number,:)));
            AoA_output(calculate_number)=AoA(calculate_number,min_p);
            %                     if calculate_number>2
            %                         if ((abs(AoA(calculate_number-1,1)-AoA_output(calculate_number-2))>20)&&(abs(AoA(calculate_number-1,2)-AoA_output(calculate_number-2))>20)&&...
            %                                 (abs(AoA(calculate_number-1,3)-AoA_output(calculate_number-2))>20)&&(abs(AoA(calculate_number-1,4)-AoA_output(calculate_number-2))>20))
            %                             AoA_output(calculate_number-1)=AoA_output(calculate_number-2);
            %                         end
            %                     end
            for particle=particle_filter
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
            end
            %                     if abs(AoA_check(calculate_number)- AoA_output(calculate_number-1))<20
            %                         AoA_output(calculate_number)=AoA_check(calculate_number)
            rssi_amp(calculate_number)=(53-mean(rssi(4)))/18;  %P0=55;Gamma=2.0
            dist(calculate_number)=10^rssi_amp(calculate_number);
            figure(1);
            plot(cos(xout(calculate_number)*pi/180)*dist(calculate_number)*100,sin(xout(calculate_number)*pi/180)*dist(calculate_number)*100,'o')
            %                         error_dist(csi_num)=((cos(AoA_output(csi_num)*pi/180)*dist(csi_num)-2.7)^2+(sin(AoA_output(csi_num)*pi/180)*dist(csi_num)-1.2)^2)^(1/2);
            csvwrite('recv_csi.csv',csv_tmp);
            csvwrite('recv_rssi.csv',rssi);
            csvwrite('recv_AoA.csv',xout(calculate_number));
        end
    end
    %     test_point(calculate_number)=asin((sin(AoA_output(calculate_number)*pi/180)*pi-init*pi)/pi)*180/pi;
end
% end
% cd('C:\Users\user\Desktop\資料\20210314')
%  save('point4.mat','save_csi','save_rssi')