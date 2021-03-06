clear;
close all;
load('C:\Users\user\Desktop\資料\20210605\output\zxc-1\7_p_num.mat')
load('C:\Users\user\Desktop\資料\20210605\clear0\zxc-1\7.mat');
load('C:\Users\user\Desktop\資料\20210605\clear0\zxc-1\7.mat');
load('C:\Users\user\Desktop\資料\20210605\output\zxc-1\7_best_removed.mat')
sub_carrier_number=1;
init=1;
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
    
    if p_num(calculate_number)==1
        figure(2);
        xlim([1 56]);ylim([-2*pi 2*pi]);hold on;title('最佳角度相位');xlabel('載波');ylabel('相位');
        plot(1:56,phase(csi_0),'r')
        plot(1:56,phase(csi_1),'g')
        plot(1:56,phase(csi_2),'b')
        tmp1(calculate_number)=csi_0(1)
        tmp3(calculate_number)=csi_2(1)
    end
    if p_num(calculate_number)==2
        figure(2);
        xlim([1 56]);ylim([-2*pi 2*pi]);hold on;title('最佳角度相位');xlabel('載波');ylabel('相位');
        plot(1:56,phase(-csi_0),'r')
        plot(1:56,phase(csi_1),'g')
        plot(1:56,phase(csi_2),'b')
        tmp1(calculate_number)=-csi_0(1)
        tmp3(calculate_number)=csi_2(1)
    end
    if p_num(calculate_number)==3
        figure(2);
        xlim([1 56]);ylim([-2*pi 2*pi]);hold on;title('最佳角度相位');xlabel('載波');ylabel('相位');
        plot(1:56,phase(-csi_0),'r')
        plot(1:56,phase(csi_1),'g')
        plot(1:56,phase(-csi_2),'b')
        tmp1(calculate_number)=-csi_0(1)
        tmp3(calculate_number)=-csi_2(1)
    end
    if p_num(calculate_number)==4        figure(2);
        xlim([1 56]);ylim([-2*pi 2*pi]);hold on;title('最佳角度相位');xlabel('載波');ylabel('相位');
        plot(1:56,phase(csi_0),'r')
        plot(1:56,phase(csi_1),'g')
        plot(1:56,phase(-csi_2),'b')
        tmp1(calculate_number)=csi_0(1)
        tmp3(calculate_number)=-csi_2(1)
    end
    
    
    figure(1);
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

figure(3);hold on
plot(1:180,angle(tmp1),'r')
plot(1:180,angle(tmp3),'b')
legend('天線1第1個載波相位','天線3第1個載波相位')

figure(4)
plot(1:180,best_AoA);xlabel('封包');ylabel('角度');ylim([-90 90]);
