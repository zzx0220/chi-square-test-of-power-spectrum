
% deciding simulation parameters
fs = 120; %Hz
T = 60; %timepoints
kkk = 4;

n = 1;H_chi2=[];H_perm=[];H_perm2=[];
amplevel = [0:0.05:0.5];
for rp = 1:1000
    for amp = 1:length(amplevel)
        for n = 1:20
            K = kkk * ones(T,n); % sample number at each timepoint

            noise = randn(60,n,kkk).*2;
            
            T_period = 12;
            x = [1:60]';
            tau = rand(1,n)*2*pi;
            %y_sin = sin((x)/5*2*pi+tau) + sin((x)/T_period*2*pi+tau).*amplevel(amp)+noise;
            
            y_sin = sin((x)/T_period*2*pi+tau).*amplevel(amp)+noise;
            %% chi2 test
            % calculate theta0 of each subject
           
            [pxx,f] = pwelch(mean(y_sin,3),ones(1,60),[],120,120);
        
            pxx = pxx(2:(end-1),:);
            f = f(2:(end-1));

            sigma2_est = var(y_sin,0,3);

            lambda = pxx*T*fs.*sum(K-1,1)./(2*sum(1./K,1).*sum((K-1).*sigma2_est,1));
            
            lambda_total = sum(2*lambda(10,:));
            p = 1-chi2cdf(lambda_total,2*n);
            if p<0.05
                H_chi2(rp,amp,n) = 1;
            else
                H_chi2(rp,amp,n) = 0;
            end
        end
    end
end

squeeze(mean(H_chi2));
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
%g.export('file_name','power_chi2test_samplebased_doublesignal','file_type','png','width',15,'height',10);
g.export('file_name','power_chi2test_samplebased_singlesignal','file_type','png','width',15,'height',10);

