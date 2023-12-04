% ========================================================================
% This is the CRF_Calibration demo in the following paper:
% Camera Response Function Calibration Algorithm Based on Multi-exposure Images
% 2023, Acta Optica Sinica ,algorithm Version 1.0
% Copyright(c) 2023, Liuzheng Gao, Banglei Guan and Ang Su
% All Rights Reserved.
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is hereby
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
% Please refer to the following paper:
% L. Gao et al., "Camera Response Function Calibration Algorithm Based on
% Multi-exposure Images,2023" Acta Optica Sinica, In process
% Please kindly report any suggestions or corrections to gaoliuzheng@126.com

clear;clc;

% 0 - adjustable parameter
highestOrder=2;    % the highestOrder for polynomial can be adjusted 

% 1 - open image files
addpath("InputFile\");
imageName=dir("InputFile\");
imageNumber=length(imageName);
firstImage=imread(imageName(3).name); % for windows, the first two names are '.' and '..'; for macOS, first number should be 4
[H,W,D]=size(firstImage);
imageN=imageNumber-2;
Images=zeros(H,W,D,imageN);
Images(:,:,:,1)=firstImage;
for cx=2:imageN
    Images(:,:,:,cx)=imread(imageName(cx+2).name);
end

% 2 - sort images from dark to bright
tempAverage=zeros(1,imageN);
for cx=1:imageN
    tempAverage(cx)=sum(Images(1:2,:,1,cx),'all');
end
[~,Index]=sort(tempAverage);

% 3 - choose calibration data
lowerLimit=20;
upperLimit=230; 
localDiff=1;
Data=cell(imageN-1,D);
for cx=1:imageN-1
    for cy=1:D
        
        choseFlag=zeros(2,256);
        for cm=2:H-1
            for cn=2:W-1
                if Images(cm,cn,cy,cx)<Images(cm,cn,cy,cx+1) && Images(cm,cn,cy,cx)>=lowerLimit && Images(cm,cn,cy,cx+1)<=upperLimit
                    Dx1=abs(Images(cm,cn,cy,cx+1)-Images(cm-1,cn,cy,cx+1));
                    Dx2=abs(Images(cm,cn,cy,cx+1)-Images(cm+1,cn,cy,cx+1));
                    Dy1=abs(Images(cm,cn,cy,cx+1)-Images(cm,cn-1,cy,cx+1));
                    Dy2=abs(Images(cm,cn,cy,cx+1)-Images(cm,cn+1,cy,cx+1));
                    dx1=abs(Images(cm,cn,cy,cx)-Images(cm-1,cn,cy,cx));
                    dx2=abs(Images(cm,cn,cy,cx)-Images(cm+1,cn,cy,cx));
                    dy1=abs(Images(cm,cn,cy,cx)-Images(cm,cn-1,cy,cx));
                    dy2=abs(Images(cm,cn,cy,cx)-Images(cm,cn+1,cy,cx));
                    
                    if    Dx1<=localDiff && Dx2<=localDiff && Dy1<=localDiff && Dy2<=localDiff ...
                            && dx1<=localDiff && dx2<=localDiff && dy1<=localDiff && dy2<=localDiff 
                       
                        choseFlag(1,Images(cm,cn,cy,cx+1)+1)=choseFlag(1,Images(cm,cn,cy,cx+1)+1)+Images(cm,cn,cy,cx);
                        choseFlag(2,Images(cm,cn,cy,cx+1)+1)=choseFlag(2,Images(cm,cn,cy,cx+1)+1)+1;
                    end
                end
            end
        end
        
        tempX=1:256;
        index=choseFlag(2,:)~=0;
        
        tempY=choseFlag(1,index)./choseFlag(2,index);
        tempZ=tempX(index)-1;
       
        Data{cx,cy}=[tempY;tempZ]';
    end
end

% 4 - calibrate CRF
loopMaxNum=50;
figure(1);
showX=0:0.05:1;

oldDelta=100;
newDelta=0;
maxDelta=0.001;
Parameters=zeros(D,highestOrder+1);
RatioSum=zeros(D,imageN-1);
for cx=1:D
    subplot(1,D,cx);box on;
    A=zeros(highestOrder*(imageN-1)+1,highestOrder);
    B=zeros(highestOrder*(imageN-1)+1,1);

    A(end,:)=1; B(end)=1;
    oldRatio=zeros(1,imageN-1);
    newRatio=zeros(1,imageN-1);
    for cy=1:imageN-1
        tempSum=mean(Data{cy,cx}(:,1)./Data{cy,cx}(:,2));
        oldRatio(cy)=tempSum;
    end
    %oldRatio(:)=0.5;
    disp(oldRatio);
    loopCont=1;    
    while loopCont<loopMaxNum     

        Cn=zeros(1,highestOrder);
        for cy=1:imageN-1
            for cm=1:highestOrder
                for cn=1:length(Data{cy,cx})
                    temp=power(Data{cy,cx}(cn,1)/255,cm)-power(Data{cy,cx}(cn,2)/255,cm)*oldRatio(cy);
                    for cz=1:highestOrder
                        Cn(cz)=Cn(cz)+temp*(power(Data{cy,cx}(cn,1)/255,cz)-power(Data{cy,cx}(cn,2)/255,cz)*oldRatio(cy));
                    end
                end
                A(cm+(cy-1)*highestOrder,:)=Cn;
                Cn=zeros(1,highestOrder);
            end
        end
        X=(A\B)';
        X=[flip(X) 0];
        newDelta=0;
        for cy=1:imageN-1
            tempY=polyval(X,Data{cy,cx}/255);
            tempSum=sum(tempY);
            newRatio(cy)=mean(tempY(:,1)./tempY(:,2));
            newDelta=newDelta+power(tempSum(1)-tempSum(2)*oldRatio(cy),2);            
        end        
        oldDelta=newDelta;
        oldRatio=newRatio;
        showY=polyval(X,showX);
        plot(showX,showY);
        pause(0.05);  % pause 50ms for showing the changing process of the curves 
        showString=['Channel:' num2str(cx) ' LoopNum:' num2str(loopCont)];
        disp(showString);
        disp('   Ratio for adjacent exposure')
        disp(newRatio);
        if newDelta<maxDelta*(imageN-1)
            break;
        end
        loopCont=loopCont+1;
    end
    Parameters(cx,:)=X;
    RatioSum(cx,:)=newRatio;
    disp('-------------------------------')
    disp('   Polynomial coefficients: Cn, Cn-1, ...,C0');   
    disp(X);
    disp('===========================================================')
end
%% curves
figure(2);hold on;box on;set(gcf,'color','white');
curveIndex={'r*-','g','b.-'};
for cx=1:D
    showE=polyval(Parameters(cx,:),showX);
    plot(showX,showE(1,:),curveIndex{cx},LineWidth=1);
end
legend('R-channel','G-channel','B-channel',Location='northwest');
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel('Pixel value / 255');
ylabel('Relative radiation')
%% Exposure ratio
figure(3);hold on;box on;set(gcf,'color','white');
showRatioX=1:imageN-1;
curveRatioIndex={'ro-','g','b*-'};
for cx=1:D    
    plot(showRatioX,RatioSum(cx,:),curveRatioIndex{cx},LineWidth=1);
end
ylim([0.4,0.6]);
legend('R-channel','G-channel','B-channel',Location='northwest');
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel('Adjacent image pairs');
ylabel('Exposure ratio');
