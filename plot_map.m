function plot_map(map)
if map==210
    original=[200,200];
    figure(1)
    hold on;
    axis equal tight;
    grid on;
    axis([original(1)-250, original(1)+510,original(2)-240, original(2)+420]);
    xlabel('cm');
    ylabel('cm');
    
    
    plot([original(1)+32,original(1)+312],[original(2)-107,original(2)-107],'k','LineWidth',2)
    plot([original(1)+32,original(1)+312],[original(2)-37,original(2)-37],'k','LineWidth',2)
    plot([original(1)+32,original(1)+312],[original(2)+33,original(2)+33],'k','LineWidth',2)
    plot([original(1)+32,original(1)+312],[original(2)+108,original(2)+108],'k','LineWidth',2)
    plot([original(1)+32,original(1)+312],[original(2)+178,original(2)+178],'k','LineWidth',2)
    plot([original(1)+32,original(1)+312],[original(2)+245,original(2)+245],'k','LineWidth',2)
    plot([original(1)+32,original(1)+312],[original(2)+315,original(2)+315],'k','LineWidth',2)
    
    plot([original(1)+32,original(1)+32],[original(2)-107,original(2)+33],'k','LineWidth',2)
    plot([original(1)+32,original(1)+32],[original(2)+108,original(2)+178],'k','LineWidth',2)
    plot([original(1)+32,original(1)+32],[original(2)+245,original(2)+315],'k','LineWidth',2)
    plot([original(1)+192,original(1)+192],[original(2)-107,original(2)+33],'k','LineWidth',2)
    plot([original(1)+192,original(1)+192],[original(2)+108,original(2)+178],'k','LineWidth',2)
    plot([original(1)+192,original(1)+192],[original(2)+245,original(2)+315],'k','LineWidth',2)
    plot([original(1)+312,original(1)+312],[original(2)-107,original(2)+33],'k','LineWidth',2)
    plot([original(1)+312,original(1)+312],[original(2)+108,original(2)+178],'k','LineWidth',2)
    plot([original(1)+312,original(1)+312],[original(2)+245,original(2)+315],'k','LineWidth',2)
    %--------------------------------
    plot([original(1)+357,original(1)+427],[original(2)-230,original(2)-230],'k','LineWidth',2)
    plot([original(1)+357,original(1)+427],[original(2)-70,original(2)-70],'k','LineWidth',2)
    plot([original(1)+357,original(1)+427],[original(2)-36,original(2)-36],'k','LineWidth',2)
    plot([original(1)+357,original(1)+427],[original(2)+84,original(2)+84],'k','LineWidth',2)
    plot([original(1)+375,original(1)+495],[original(2)+84,original(2)+84],'k','LineWidth',2)
    plot([original(1)+375,original(1)+495],[original(2)+154,original(2)+154],'k','LineWidth',2)
    plot([original(1)+375,original(1)+495],[original(2)+230,original(2)+230],'k','LineWidth',2)
    plot([original(1)+375,original(1)+495],[original(2)+300,original(2)+300],'k','LineWidth',2)
    
    plot([original(1)+357,original(1)+357],[original(2)-230,original(2)-70],'k','LineWidth',2)
    plot([original(1)+357,original(1)+357],[original(2)-36,original(2)+84],'k','LineWidth',2)
    plot([original(1)+427,original(1)+427],[original(2)-230,original(2)-70],'k','LineWidth',2)
    plot([original(1)+427,original(1)+427],[original(2)-36,original(2)+84],'k','LineWidth',2)
    plot([original(1)+375,original(1)+375],[original(2)+84,original(2)+154],'k','LineWidth',2)
    plot([original(1)+375,original(1)+375],[original(2)+230,original(2)+300],'k','LineWidth',2)
    plot([original(1)+495,original(1)+495],[original(2)+84,original(2)+154],'k','LineWidth',2)
    plot([original(1)+495,original(1)+495],[original(2)+230,original(2)+300],'k','LineWidth',2)
    %------------------
    plot([original(1)-250,original(1)+505],[original(2)+413,original(2)+413],'k','LineWidth',4)
    plot([original(1)-250,original(1)+505],[original(2)-230,original(2)-230],'k','LineWidth',4)
    
    plot([original(1)+505,original(1)+505],[original(2)-230,original(2)+413],'k','LineWidth',4)
    
    txt = {'Desk'};
    text(original(1)+90,original(2)+0,txt,'Color',[0.7, 0.7 ,0.7])
    text(original(1)+90,original(2)-70,txt,'Color',[0.7, 0.7 ,0.7])
    text(original(1)+230,original(2)+0,txt,'Color',[0.7, 0.7 ,0.7])
    text(original(1)+230,original(2)-70,txt,'Color',[0.7, 0.7 ,0.7])
    text(original(1)+90,original(2)+145,txt,'Color',[0.7, 0.7 ,0.7])
    text(original(1)+230,original(2)+145,txt,'Color',[0.7, 0.7 ,0.7])
    text(original(1)+90,original(2)+282,txt,'Color',[0.7, 0.7 ,0.7])
    text(original(1)+230,original(2)+282,txt,'Color',[0.7, 0.7 ,0.7])
    text(original(1)+410,original(2)+119,txt,'Color',[0.7, 0.7 ,0.7])
    text(original(1)+410,original(2)+265,txt,'Color',[0.7, 0.7 ,0.7])
    text(original(1)+365,original(2)-150,txt,'Color',[0.7, 0.7 ,0.7])
    text(original(1)+365,original(2)+20,txt,'Color',[0.7, 0.7 ,0.7])
    
        
    
    rx_table=original+[
        -200,50;450,-150;450,150;100,-50;350,350;-50,350;-200,-200;300,-200;200,200
        ];
    plot(rx_table(1,1),rx_table(1,2),'r','Marker','x');t1=text(rx_table(1,1)-80,rx_table(1,2)-20,'Rx1(zxc-1)');set(t1, 'Color',[1, 0 ,0]);
    plot(rx_table(2,1),rx_table(2,2),'r','Marker','x');t1=text(rx_table(2,1)-80,rx_table(2,2)-20,'Rx2(user)');set(t1, 'Color',[1, 0 ,0]);
    plot(rx_table(3,1),rx_table(3,2),'r','Marker','x');t1=text(rx_table(3,1)-80,rx_table(3,2)-20,'Rx3(zxc-0)');set(t1, 'Color',[1, 0 ,0]);
    plot(rx_table(4,1),rx_table(4,2),'r','Marker','x');t1=text(rx_table(4,1)-80,rx_table(4,2)-20,'Rx4(zxc-1)');set(t1, 'Color',[1, 0 ,0]);
    plot(rx_table(5,1),rx_table(5,2),'r','Marker','x');t1=text(rx_table(5,1)-80,rx_table(5,2)-20,'Rx5(user)');set(t1, 'Color',[1, 0 ,0]);
    plot(rx_table(6,1),rx_table(6,2),'r','Marker','x');t1=text(rx_table(6,1)-80,rx_table(6,2)-20,'Rx6(zxc-0)');set(t1, 'Color',[1, 0 ,0]);
    plot(rx_table(7,1),rx_table(7,2),'r','Marker','x');t1=text(rx_table(7,1)-80,rx_table(7,2)-20,'Rx7(zxc-0)');set(t1, 'Color',[1, 0 ,0]);
    plot(rx_table(8,1),rx_table(8,2),'r','Marker','x');t1=text(rx_table(8,1)-80,rx_table(8,2)-20,'Rx8(zxc-0)');set(t1, 'Color',[1, 0 ,0]);
    plot(rx_table(9,1),rx_table(9,2),'r','Marker','x');t1=text(rx_table(9,1)-80,rx_table(9,2)-20,'Rx9(zxc-0)');set(t1, 'Color',[1, 0 ,0]);
    %     t1=text(-15,-20,'Rx1')
    %     set(t1, 'Color',[1, 0 ,0]);
    
    %     t2=text(335,-20,'Rx2')
    %     set(t2, 'Color',[1, 0 ,0]);
    %
    %     t3=text(335,180,'Tx2')
    %     set(t3, 'Color',[0, 0 ,1]);
    
    %     t1=text(120,-155,'point1')
    %     set(t1, 'Color',[0, 0 ,1]);
    %     t1=text(360,-215,'point2')
    %     set(t1, 'Color',[0, 0 ,1]);
    
    %%Rx label
    %     plot(0,0,'r','Marker','x');
    %     plot(350,0,'r','Marker','x');
    %     plot(350,200,'b','Marker','o');
    %     plot(100,-140,'b','Marker','x');
    %     plot(380,-200,'b','Marker','x');
    
    %%Tx label
    tx_table=original+[
        -50,-200;-50,-150;-50,-100;-50,-50;-50,0;-50,50;-50,100;-50,150;-50,200;-50,250;
        0,250;0,200;0,150;0,100;0,50;0,0;0,-50;0,-100;0,-150;0,-200;
        50,-200;100,-200;150,-200;200,-200;250,-200;250,-150;200,-150;150,-150;100,-150;50,-150;
        350,-200;350,-150;350,-100;350,-50;350,0;350,50;350,100;350,150;350,200;350,250;
        300,200;250,200;150,200;100,200;50,200;50,100;100,100;150,100;200,100;250,100;
        300,100;300,50;250,50;200,50;150,50;100,50;50,50;
        ];
    for num=1:size(tx_table,1)
        plot(tx_table(num,1),tx_table(num,2),'b','Marker','o');p(num)=text(tx_table(num,1)-10,tx_table(num,2)-10,num2str(num));set(p(num), 'Color',[0.81,0.55,0.16]);
    end
    for num=1:size(tx_table,1)
        AoA_table(1,num)=atan2(tx_table(num,2)-rx_table(1,2),tx_table(num,1)-rx_table(1,1))/pi*180;%atan2((y2-y1),(x2-x1))
        AoA_table(2,num)=atan2(-(tx_table(num,1)-rx_table(2,1)),tx_table(num,2)-rx_table(2,2))/pi*180-90;
        AoA_table(3,num)=atan2(-(tx_table(num,1)-rx_table(3,1)),tx_table(num,2)-rx_table(3,2))/pi*180-90;
        AoA_table(4,num)=-atan2(tx_table(num,1)-rx_table(4,1),abs(tx_table(num,2)-rx_table(4,2)))/pi*180;%atan2((y2-y1),(x2-x1))
        AoA_table(5,num)=-atan2(abs(tx_table(num,2)-rx_table(5,2)),tx_table(num,1)-rx_table(5,1))/pi*180+90;
        AoA_table(6,num)=-atan2(abs(tx_table(num,2)-rx_table(6,2)),tx_table(num,1)-rx_table(6,1))/pi*180+90;
        AoA_table(7,num)=atan2(tx_table(num,2)-rx_table(7,2),tx_table(num,1)-rx_table(7,1))/pi*180;%atan2((y2-y1),(x2-x1))
        AoA_table(8,num)=-atan2(tx_table(num,1)-rx_table(8,1),abs(tx_table(num,2)-rx_table(8,2)))/pi*180;%atan2((y2-y1),(x2-x1))
        AoA_table(9,num)=-atan2(abs(tx_table(num,2)-rx_table(9,2)),tx_table(num,1)-rx_table(9,1))/pi*180+90;
    end
    %%角度分布
%     figure(2);hold on;
%     for num=1:9
%         cdfplot(AoA_table(num,1:57))
%     end
%     legend('Rx1','Rx2','Rx3','Rx4','Rx5','Rx6','Rx7','Rx8','Rx9');title('57點角度分布');xlabel('角度');ylabel('CDF')
end
end