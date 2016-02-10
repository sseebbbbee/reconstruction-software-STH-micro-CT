function [ s,n ] = make_simulated_images( image_width,image_height,y_0,x_0,nr_im,stand)

s=zeros(image_width*image_height,nr_im);
n=zeros(image_width*image_height,nr_im);
I=6*ones(image_height,image_width);
signal_width=3;
sig = zeros(image_width,image_height);
sig(floor(image_height/2)-signal_width+y_0:floor(image_height/2)+signal_width+y_0,floor(image_width/2)-signal_width+x_0:floor(image_width/2)+signal_width+x_0)=1;

for ii=1:nr_im

I_noise =I+ stand*randn(image_width); %randn has mean 0 and std 1
n(:,ii)=I_noise(:);
I_noise2 =I+ stand* randn(image_width);
I_Sig = I_noise2+ 15*sig;
s(:,ii)=I_Sig(:);
end
%figure, imagesc(I_noise),colormap(gray)
mean(I_noise(:))
std(I_noise(:))
yx=I_Sig(floor(image_height/2)-signal_width+y_0:floor(image_height/2)+signal_width+y_0,floor(image_width/2)-signal_width+x_0:floor(image_width/2)+signal_width+x_0);
a=mean(yx(:))
xy=I_noise(floor(image_height/2)-signal_width+y_0:floor(image_height/2)+signal_width+y_0,floor(image_width/2)-signal_width+x_0:floor(image_width/2)+signal_width+x_0);
b=mean(xy(:))
c=std(xy(:))
CNR=(a-b)/c

figure, imagesc(I_noise),colormap(gray),set(gca, 'CLim', [-160, 240]);
figure, imagesc(I_Sig),colormap(gray),set(gca, 'CLim', [-160, 240]);

end

