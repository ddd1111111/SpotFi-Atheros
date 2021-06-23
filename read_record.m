%----------------------------RESET-----------------------------------------
clear;
close all;
load("C:\Users\user\Desktop\20200605\save_data\recv_csi_data_5G_3M_1.mat")
format short g;
% init_check=0;
calculate_number=0;
Q_candidate = 1;%[220 440 660 880]%[1000 2000 3000 4000 5000];
MUSIC_iter = 10;
rssi=0;
sample_csi_traceTmp_r=zeros(1,168);
figuer_plot_range=[0, 6, -3, 3];
position_x(1)=0;
position_y(1)=0;
true_=0
real_location=[1,0]; %Tx(X,Y) location
for c_time=1:252
% for c_time=1:size(save_csi_data,2)
    for i=1:1
    tic   %time
    % pause(0.8)
    calculate_number=c_time;
    % cd("C:\Users\user\Desktop\20191205\Server\python")
    % csv_tmp=csvread('csi_tmp.csv');
    % for num=1:MUSIC_iter
    %     sample_csi_traceTmp_r(num,:)=csv_tmp(1,((num*168)-167):(num*168));
    % end
    % save_data{calculate_number}= sample_csi_traceTmp_r+
    %----------------------------RESHAPE_DATA-----------------------------------------
    rssi=save_rssi_data{1,(c_time)};
    sample_csi_traceTmp_r=save_csi_data{1,(c_time)};
    
    cd("C:\Users\user\Desktop\20200505\SpotFi-Atheros")
    target_antenna = 1;
    %----------------------------MUSIC_AOA-----------------------------------------
    for q = Q_candidate

        %fc = 5.63e9; % center frequency
        % fc = record{1}.channel*10^6;
        fc =5220*10^6;
        M = 3;    % number of rx antennas
        %fs = 40e6; % channel bandwidth
        fs = 20e6;
        c = 3e8;  % speed of light
        % d = 6.2e-2;  % distance between adjacent antennas in the linear antenna array
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
        paramRange.GridPts = [101 101 1]; % number of grid points in the format [number of grid points for ToF (Time of flight), number of grid points for angle of arrival (AoA), 1]
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
        paramRange.generateAtot = 0; %2 to 0 without saving mat file
        
        sample_csi_trace_m = sample_csi_traceTmp_r(q:q+MUSIC_iter-1,:);
        if i==1
            sample_csi_trace_m=sample_csi_trace_m;
            
        end
        if i==2
                sample_csi_trace_m=[-sample_csi_trace_m(1:56) ,sample_csi_trace_m(57:168)];
                
        end
        if i==3
                 sample_csi_trace_m=[sample_csi_trace_m(1:56) ,-sample_csi_trace_m(57:112),sample_csi_trace_m(113:168)];
                
        end
        if i==4
                sample_csi_trace_m=[sample_csi_trace_m(1:112),-sample_csi_trace_m(113:168)];
               
        end
        
        for ind = 1:MUSIC_iter
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
        
        % MUSIC algorithm for estimating angle of arrival
        % aoaEstimateMatrix is (nComps x 5) matrix where nComps is the number of paths in the environment. First column is ToF in ns and second column is AoA in degrees as defined in SpotFi paper
        aoaEstimateMatrix = backscatterEstimationMusic(sample_csi_trace_sanitized_m, M, N, c, fc,...
            T, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, MUSIC_iter, ones(2));

    end
    %----------------------------CALCULATE_DISTANCE-----------------------------------------
    rssi_amp(q)=(70-mean(rssi))/27;  %P0=70;Gamma=2.7
    dist(q)=10^rssi_amp(q);

    %%----------------------------FIND_AOA-----------------------------------------
    if(aoaEstimateMatrix(1,1)<aoaEstimateMatrix(2,1))
        position_x(calculate_number)=cos(aoaEstimateMatrix(1,2)*pi/180)*dist(q);  %position_x(1)=cos(aoaEstimateMatrix(1,2)*pi/180)*dist(q)
        position_y(calculate_number)=sin(aoaEstimateMatrix(1,2)*pi/180)*dist(q); %position_y(1)=sin(aoaEstimateMatrix(1,2)*pi/180)*dist(q)
        dist_error(q)=((cos(aoaEstimateMatrix(1,2)*pi/180)*dist(q)-1)^2+(sin(aoaEstimateMatrix(1,2)*pi/180)*dist(q))^2)^(1/2);
        mini(c_time,i)=(aoaEstimateMatrix(1,2));
    else
        position_x(calculate_number)=cos(aoaEstimateMatrix(2,2)*pi/180)*dist(q);   %position_x(1)=cos(aoaEstimateMatrix(2,2)*pi/180)*dist(q)
        position_y(calculate_number)=sin(aoaEstimateMatrix(2,2)*pi/180)*dist(q);   %position_y(1)=sin(aoaEstimateMatrix(2,2)*pi/180)*dist(q)
        dist_error(q)=((cos(aoaEstimateMatrix(2,2)*pi/180)*dist(q)-1)^2+(sin(aoaEstimateMatrix(2,2)*pi/180)*dist(q))^2)^(1/2);
        mini(c_time,i)=(aoaEstimateMatrix(2,2));
    end
     %----------------------------PLOT_FIGURE-----------------------------------------
    if (position_x(1)~=0&&position_y(1)~=0)
        figure(3);
        hold on;
        axis(figuer_plot_range);
        plot(0,0,'+');
        for plot_num=calculate_number %for plot_num=1:1
            plot( position_x(plot_num),position_y(plot_num),'o');
        end
        error_10_pack(calculate_number)= ((position_x(calculate_number)-real_location(1))^2+(position_y(calculate_number)-real_location(2))^2)^(1/2);
    end
%     if position_x(calculate_number)>0.5&&position_x(calculate_number)<1.5&&position_y(calculate_number)>-0.5&&position_y(calculate_number)<0.5
%         true_=true_+1
%     end
    toc
    end
end
% for j=1:181
%     [tmp(j,1:i) tmp1(j,1:i)]=sort(mini(j,1:i));
% end
% mean(tmp)
% rms(tmp)
% for i=1:181
% tmp(i)=mini(i,a(i))
% end
% for i=1:181
%     hold on;
%     axis(figuer_plot_range);
%     rssi=save_rssi_data{1,i};
%     rssi_amp(i)=(70-mean(rssi))/27
%     dist(i)=10^rssi_amp(i);
%     position_x(i)=cos(tmp(i)*pi/180)*dist(i)
%     position_y(i)=sin(tmp(i)*pi/180)*dist(i)
%     plot( position_x(i),position_y(i),'o');
% end