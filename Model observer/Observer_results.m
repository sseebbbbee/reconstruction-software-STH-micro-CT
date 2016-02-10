function [task_SNR,tpf,fpf,AUC] = Observer_results( t_s,t_n )
% calculating the performance parameters
task_SNR=sqrt( ((mean(t_s) - mean(t_n)).^2) / (.5*var(t_s) + .5*var(t_n)));
if max(t_s(:)) >=max(t_n(:))
    max_threshold=max(t_s(:));
else
    max_threshold=max(t_n(:));
end
if min(t_s(:)) <=min(t_n(:))
    min_threshold=min(t_s(:));
else
    min_threshold=min(t_n(:));
end

thresholds=linspace(min_threshold,max_threshold,length(t_s)+length(t_n));
tpf = zeros(length(thresholds),1);
fpf = zeros(length(thresholds),1);



for ii = 1:length(thresholds)
    tpf(ii)=sum(t_s >= thresholds(ii))/length(t_s);
    fpf(ii)=sum(t_n >= thresholds(ii))/length(t_n);
end


AUC = -trapz(fpf,tpf);

end

