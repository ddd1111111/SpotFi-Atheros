clear;close all;
ground_truth_AoA=[-59.0362434679265,-53.1301023541560,-45,-33.6900675259798,-18.4349488229220,0,...
    18.4349488229220,33.6900675259798,45,53.1301023541560,45,36.8698976458440,26.5650511770780,...
    14.0362434679265,0,-14.0362434679265,-26.5650511770780,-36.8698976458440,-45,-51.3401917459099,...
    -45,-39.8055710922652,-35.5376777919744,-32.0053832080835,-29.0546040990771,-23.9624889745782,...
    -26.5650511770780,-29.7448812969422,-33.6900675259798,-38.6598082540901,-24.4439547804165,...
    -19.9831065219000,-15.2551187030578,-10.3048464687660,-5.19442890773481,0,5.19442890773481,...
    10.3048464687660,15.2551187030578,19.9831065219000,16.6992442339936,18.4349488229220,...
    23.1985905136482,26.5650511770780,30.9637565320735,11.3099324740202,9.46232220802562,...
    8.13010235415598,7.12501634890180,6.34019174590991,5.71059313749964,0,0,0,0,0,0];
for for_5pt=14
    if ground_truth_AoA(for_5pt)>0
        init=1;
    else
        init=-1;
    end
    point_num=string(for_5pt);
    load('C:\Users\user\Desktop\資料\20210605\output\zxc-1\'+point_num+'_find_pmusic_.mat')
    figure(1);
    subplot(3,1,1);hold on;
    plot(1:180,angle(tmp1),'r')
    plot(1:180,angle(tmp3),'b')
    xlabel('封包');ylabel('相位');legend('天線1第1個載波相位','天線3第1個載波相位');title(['第',num2str(for_5pt),'點(實際擺設角度:',num2str(ground_truth_AoA(for_5pt)),'度)']);
    
    subplot(3,1,2);hold on;
    plot(1:180,z(:,1),'r')
    plot(1:180,z(:,2),'g')
    plot(1:180,z(:,3),'b')
    plot(1:180,z(:,4),'k')
    xlabel('封包');ylabel('pmusic');legend('組合1Pmusic值','組合2Pmusic值','組合3Pmusic值','組合4Pmusic值');title('1,3重疊;2,4重疊');
    
    subplot(3,1,3);hold on;ylim([-90 90])
    plot(1:180,ground_truth_AoA(for_5pt)*ones(1,180),'m:','LineWidth',2)
    plot(1:180,AoA(:,1),'r')
    plot(1:180,AoA(:,2),'g')
    plot(1:180,AoA(:,3),'b')
    plot(1:180,AoA(:,4),'k')
    xlabel('封包');ylabel('估計角度');legend('實際角度','組合1','組合2','組合3','組合4')
    
    load('C:\Users\user\Desktop\資料\20210605\clear0\zxc-1\'+point_num+'.mat')
    figure(2);
    sub_carrier_number=1;
    for calculate_number=1:180
        csv_tmp=[reshape(file{calculate_number*10,1}.csi(1,1,1:56),1,56) reshape(file{calculate_number*10,1}.csi(2,1,1:56),1,56) reshape(file{calculate_number*10,1}.csi(3,1,1:56),1,56)];
        
        phase_csi_1=phase(csv_tmp(57:112));
        abs_csi_1=abs(csv_tmp(57:112));
        phase_csi_0=phase(csv_tmp(1:56));
        abs_csi_0=abs(csv_tmp(1:56));
        phase_csi_0=phase_csi_0-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56))-0.00184675654972853*[1:56]-(0.688313191228861);
        %         phase_csi_0=phase_csi_0-0.00184675654972853*[1:56]-(0.688313191228861);
        
        [csi_0_x,csi_0_y]=pol2cart(phase_csi_0,1);%pol2cart:Polar to Cartesian Coordinates
        csi_0=[csi_0_x+i*csi_0_y];
        
        phase_csi_2=phase(csv_tmp(113:168));
        abs_csi_2=abs(csv_tmp(113:168));
        phase_csi_2=phase_csi_2-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56))+0.000208698745420363*[1:56]-(-0.285393881253289)-pi/2;
        %         phase_csi_2=phase_csi_2+0.000208698745420363*[1:56]-(-0.285393881253289)-pi/2;
        
        [csi_2_x,csi_2_y]=pol2cart(phase_csi_2,1);
        csi_2=[csi_2_x+i*csi_2_y];
        
        phase_csi_1=phase(csv_tmp(57:112));
        phase_csi_1=phase_csi_1-phase(csv_tmp(56+sub_carrier_number)).*ones(1,56)-(phase_csi_1-phase_csi_1(sub_carrier_number).*ones(1,56));
        %         phase_csi_1=phase_csi_1;
        [csi_1_x,csi_1_y]=pol2cart(phase_csi_1,1);
        csi_1=[csi_1_x+i*csi_1_y];
        
        figure(2);
        subplot(1,2,1);
        xlim([1 56]);ylim([-2*pi 2*pi]);hold on;title('原始');xlabel('載波');ylabel('相位');
        plot(1:56,phase(csi_0),'r')
        plot(1:56,phase(csi_1),'g')
        plot(1:56,phase(csi_2),'b')
        %
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
        subplot(1,2,2)
        xlim([1 56]);ylim([-2*pi 2*pi]);hold on;title('挑選後');xlabel('載波');ylabel('相位');
        plot(1:56,phase(csv_tmp(1:56)),'r')
        plot(1:56,phase(csv_tmp(57:112)),'g')
        plot(1:56,phase(csv_tmp(113:168)),'b')
    end
end