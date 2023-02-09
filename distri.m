for i = 1:1000
    freq = normrnd(12,1);
    T = 60/freq*2;
    x = [1:60]';
    tau = rand*2*pi;
    y_sin = sin((x)/T*2*pi+tau);
    [pxx(:,i),f]=pwelch(y_sin,ones(1,60),[],60,120);
end

pxx_prop = pxx(2:31,:)./sum(pxx(2:31,:),1);

% 
% figure;clear g
% g = gramm('x',2:2:60,'y',[mean(pxx_prop,2)'./2;normpdf(2:2:60,12,1)],'color',{'sim','true'});
% g.geom_line();
% g.set_names('x','frequency','y','power','color','data');
% g.draw();
% g.export('file_name','psd_est_norm','file_type','png','width',10,'height',7.5);


%% generate distribution under null hypothesis
prop_null = 1/30;
prop_r = mean(pxx_prop,2);
% prop_r = zeros(30,1)*;
% prop_r(6) = 1;



for n = 1:20
    pxx_rp=[];
for i = 1:1000
    noise = randn(60,n);
    [pxx_rp(:,:,i),f]=pwelch(noise,ones(1,60),[],60,120);
end
f = f(2:31);
band = (f>=10 & f<=14);

pxx_prop_rp = pxx_rp(2:31,:,:)./sum(pxx_rp(2:31,:,:),1);

lambda_rp = squeeze(sum(sum(pxx_prop_rp(band,:,:)./prop_r(band),1),2) - sum(sum(pxx_prop_rp(band,:,:)./prop_null,1),2));

%figure;histogram(lambda_rp);

%% simulate data
amplevel = [0:0.05:0.5];
for amp = 1:length(amplevel)
    for rp = 1:1000
        noise = randn(60,n);
        freq = 12;%normrnd(12,1);
        T = 60/freq*2;
        x = [1:60]';
        tau = rand*2*pi;
        y_sin = sin((x)/T*2*pi+tau).*amplevel(amp)+noise;
        [pxx,f]=pwelch(y_sin,ones(1,60),[],60,120);
    
        pxx_prop = pxx(2:31,:)./sum(pxx(2:31,:),1);
        lambda = squeeze(sum(sum(pxx_prop(band,:)./prop_r(band),1),2) - sum(sum(pxx_prop(band,:)./prop_null,1),2));
        
        p = sum(lambda>lambda_rp)./1000;
        if p<0.05
            H_chi2(rp,amp,n) = 1;
        else
            H_chi2(rp,amp,n) = 0;
        end
    end
end
end

H = permute(H_chi2,[2 1 3]);
H = reshape(H,[size(H,1) size(H,2)*size(H,3)]);
H = H';
H_t = H+rand*0.01;
subjnum = repelem(1:n,1,1000);

figure;clear g
g = gramm('x',amplevel,'y',H,'lightness',subjnum);
g.stat_summary('geom',{'line','errorbar'},'type','sem','setylim','yes');
g.set_names('x','amplitude','y','power','lightness','sample number');
g.set_color_options("legend","merge");
g.geom_hline('yintercept',0.8,'style','k--');
g.draw();
g.export('file_name','power_shapetest','file_type','png','width',15,'height',10);
