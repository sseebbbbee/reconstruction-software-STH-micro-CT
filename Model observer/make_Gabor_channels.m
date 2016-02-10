function [ U] = make_Gabor_channels(theta,f_c,w_s,beta,x_0,y_0,image_width,image_height)
% creating the Gabor channels
% w_s is the channel width
% f_c is the central frequency
% theta is the orientation
% beta is the phase factor.
% x_0 and y_0 is the center of the gabor function
xmax = floor(image_width/2);
ymax= floor(image_height/2);
xmin = -xmax;
ymin = -ymax;
[x,y] = meshgrid(xmin:xmax,ymin:ymax);

ii=1;
U=zeros(image_width*image_height,length(theta)*length(f_c)*length(beta));
for t=1:length(theta)
    for p=1:length(f_c)
        
        for b=1:length(beta)
            gabor= exp(-(4*log(2).*((x-x_0).^2+(y-y_0).^2)/(w_s(p)^2))).*cos(2*pi*f_c(p).*((x-x_0).*cos(theta(t))+(y-y_0).*sin(theta(t)))+beta(b));
            U(:,ii)=gabor(:);
            ii=ii+1;
        end
        
    end
end
%figure,imagesc(gabor)
end

