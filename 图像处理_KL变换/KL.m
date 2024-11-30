clear;clc;

% 读取图像与处理
tif=double(imread("tmpicture.tif"));

[M,N,bands]=size(tif);
tif=reshape(tif, [], bands);

tif=(tif-mean(tif))./std(tif);

% 求解协方差矩阵
cov_matrix=cov(tif);

% 求解特征值D和特征向量X
[X,D]=eig(cov_matrix);
D=diag(D);

% 变换
new_tif=tif*X;
new_tif=reshape(new_tif,M,N,[]);

% 展示变化后的6个分量

figure
for i=1:bands
    subplot(2,3,i)
    imshow(new_tif(:,:,bands-i+1),[])
    title("特征值："+num2str(D(bands-i+1)))
end


% 取前K组最大的特征值对应的特征向量
k=3;
[~,index]=sort(D,'descend');
P=X(:,index(1:k));

% 旋转变换
new_tif=tif*P;

% 绘制特征值曲线
figure
hold on
xlabel("特征值序号")
ylabel("特征值（标准化）")
plot(D(index), ...
    LineStyle="-", ...
    Color="r", ...
    LineWidth=2, ...
    Marker="o")
hold off

% 将KL变换得到的三个分量进行假彩色增强输出
figure
imshow(reshape(new_tif,M,N,[]),[])
% plot(new_tif(:,1),new_tif(:,2),"*")






