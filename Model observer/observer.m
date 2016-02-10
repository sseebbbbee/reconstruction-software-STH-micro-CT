clc, clear all, close all

stand=1; %standard deviation of the noise
image_width=129;
image_height=129;
% x_0=0;
% y_0=5;
nr_im=150;
x_pos=[0];
y_pos=[0];
s=zeros(image_width*image_height,nr_im,length(x_pos)*length(y_pos));
n=zeros(image_width*image_height,nr_im,length(x_pos)*length(y_pos));
sig_location=zeros(length(x_pos)*length(y_pos),2);
ii=1;
for xx=1:length(x_pos)
    for yy=1:length(y_pos)
        x_0=x_pos(xx);
        y_0=y_pos(yy);
        sig_location(ii,:)=[x_0,y_0];
        [ s(:,:,ii),n(:,:,ii) ] = make_simulated_images( image_width,image_height,y_0,x_0,nr_im,stand);
        ii=ii+1;
    end
end
%%

%parameters to design the Gabor channels
theta=[0,2*pi/5,4*pi/5,6*pi/5,8*pi/5];
f_c=[3/128,3/64,3/32,3/16];
w_s=[56.48, 28.24, 14.12,7.06];
beta=[0,pi/2];
x_search_loc=[-1:1:1];
y_search_loc= [-1:1:1];
t_s=zeros(length(x_search_loc)*length(y_search_loc),nr_im,size(sig_location,1));
t_n=zeros(length(x_search_loc)*length(y_search_loc),nr_im,size(sig_location,1));
search_location=zeros(length(x_search_loc)*length(x_search_loc),2);

jj=0;
for xx= 1:length(x_search_loc);
    for yy= 1:length(y_search_loc);
        x_0=x_search_loc(xx);
        y_0=y_search_loc(yy);
        
        jj=jj+1;
        search_location(jj,:)=[x_0,y_0];
        [U] = make_Gabor_channels(theta,f_c,w_s,beta,x_0,y_0,image_width,image_height);
        for ii=1:size(sig_location,1)
            [ t_s(jj,:,ii), t_n(jj,:,ii)] = CHO( s(:,:,ii), n(:,:,ii), U );
        end
    end
end

[t_s_max,ind_s]=max(t_s);
[t_n_max,ind_n]=max(t_n);
std_t_n=std(t_n_max);
random_numbers= -1 + (2).*rand(1,nr_im); %creating random number -1 to 1
a=6;
t_s_max=t_s_max + a*std_t_n*random_numbers;
t_n_max=t_n_max + a*std_t_n*random_numbers;
for  ii=1:size(sig_location,1)
    [task_SNR(ii),tpf(:,ii),fpf(:,ii),AUC(ii)] = Observer_results( t_s_max(:,:,ii),t_n_max(:,:,ii) );
end
task_SNR_m=mean(task_SNR);
AUC_m=mean(AUC);
tpf_m=mean(tpf,2);
fpf_m=mean(fpf,2);

figure,
plot(fpf_m,tpf_m);xlabel('FPF');ylabel('TPF')
AUC_m2 = -trapz(fpf_m,tpf_m);



