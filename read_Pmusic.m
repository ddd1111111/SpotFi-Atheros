clear;
%zxc-1
ground_truth_AoA=[-59.0362434679265,-53.1301023541560,-45,-33.6900675259798,-18.4349488229220,0,...
    18.4349488229220,33.6900675259798,45,53.1301023541560,45,36.8698976458440,26.5650511770780,...
    14.0362434679265,0,-14.0362434679265,-26.5650511770780,-36.8698976458440,-45,-51.3401917459099,...
    -45,-39.8055710922652,-35.5376777919744,-32.0053832080835,-29.0546040990771,-23.9624889745782,...
    -26.5650511770780,-29.7448812969422,-33.6900675259798,-38.6598082540901,-24.4439547804165,...
    -19.9831065219000,-15.2551187030578,-10.3048464687660,-5.19442890773481,0,5.19442890773481,...
    10.3048464687660,15.2551187030578,19.9831065219000,16.6992442339936,18.4349488229220,...
    23.1985905136482,26.5650511770780,30.9637565320735,11.3099324740202,9.46232220802562,...
    8.13010235415598,7.12501634890180,6.34019174590991,5.71059313749964,0,0,0,0,0,0];

for threshold=15:-1:10
for for_5pt=1:57
    tmp=0;
    point_num=string(for_5pt);
    load('C:\Users\user\Desktop\資料\20210605\output\zxc-1\'+point_num+'_find_pmusic.mat');
    for calculate_num=1:180
        if min(z(calculate_num,:))>threshold
            tmp=[tmp calculate_num];
        end
    end
    AoA(tmp(2:size(tmp,2)),:);
    num(for_5pt)=size(tmp,2)-1;
    tmp_1_57(for_5pt,1:(size(tmp,2)-1))=tmp(2:size(tmp,2));
    %     for calculate_number=tmp(2:size(tmp,2))
    %         [min_error,p_num(calculate_number)]=min(abs(AoA(calculate_number,1:4)-ground_truth_AoA));
    %         best_AoA(calculate_number)=AoA(calculate_number,p_num(calculate_number));
    %     end
end
figure(1);
subplot(3,1,1);hold on;
plot(1:57,num);

tmp_min_error=0;
for for_5pt=1:57
    point_num=string(for_5pt);
    load('C:\Users\user\Desktop\資料\20210605\output\zxc-1\'+point_num+'_find_pmusic.mat');
    load('C:\Users\user\Desktop\資料\20210605\output\zxc-1\'+point_num+'_removed.mat');
    for num=1:size(tmp_1_57,2)
        if tmp_1_57(for_5pt,num)~=0
%             [min_error,p_num(num)]=min(abs(AoA(tmp_1_57(for_5pt,num),1:4)-ground_truth_AoA(for_5pt)));
%             best_AoA(num)=AoA(num,p_num(num));
%             tmp_min_error=[tmp_min_error min_error];
            min_error=abs(AoA_output(num)-ground_truth_AoA(for_5pt));
            tmp_min_error=[tmp_min_error min_error];
        end
    end
end
subplot(3,1,2);hold on;
cdfplot(tmp_min_error(2:size(tmp_min_error,2)))
size_min_error(threshold)=size(tmp_min_error,2)-1;
end
subplot(3,1,1);
legend('threshold=10','threshold=11','threshold=12','threshold=13','threshold=14','threshold=15');xlabel('位置');ylabel('數量')

subplot(3,1,3);
plot(10:15,size_min_error(10:15));xlabel('threshold');ylabel('封包數量');

for point_number=1:57
    point_num=string(point_number);
    load('C:\Users\user\Desktop\資料\20210605\output\zxc-1\'+point_num+'_removed.mat');
    tmp((180*(point_number-1)+1):(180*point_number))=abs(AoA_output-ground_truth_AoA(point_number));
%     cdfplot(abs(best_AoA-ground_truth(point_number)));
end
subplot(3,1,2);
a3=cdfplot(tmp);
set(a3,'Linewidth',2)
% legend('threshold=10','threshold=11','threshold=12','threshold=13','threshold=14','threshold=15','best');xlabel('誤差角度');ylabel('CDF');xlim([0 90]);
legend('threshold=15','threshold=14','threshold=13','threshold=12','threshold=11','threshold=10','best');xlabel('誤差角度');ylabel('CDF');xlim([0 90]);