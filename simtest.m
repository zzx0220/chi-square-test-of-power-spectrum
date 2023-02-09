
n = 1;H_chi2=[];H_perm=[];H_perm2=[];
amplevel = [0.5 0.6 0.7 0.8 0.9 1];
for rp = 1:100
    for amp = 1:length(amplevel)
        noise = randn(60,n);
        
        T = 10;
        x = [1:60]';
        tau = 3;%rand(1,n)*2*pi;
        y_sin = sin((x)/T*2*pi+tau).*amplevel(amp)+noise;
    
%         %% F test
%         [pxx,f] = pwelch(y_sin,ones(1,60),[],60,120);
%         pxx = pxx(2:31,:);
%         f = f(2:31);
%         
%         pxx = pxx./sum(pxx,1)/2*60;
%         lambda = sum(pxx(6,:));
%         p = 1-fcdf(lambda,2,60);
%         if p<0.05
%             H_f(rp,amp) = 1;
%         else
%             H_f(rp,amp) = 0;
%         end

        %% chi2 test
        % calculate theta0 of each subject
        [pxx,f] = pwelch(y_sin,ones(1,60),[],60,120);
        pxx = pxx(2:31,:);
        f = f(2:31);
        
        pxx = pxx./mean(pxx,1)*2;
        lambda = sum(pxx(6,:));
        p = 1-chi2cdf(lambda,2*n);
        if p<0.05
            H_chi2(rp,amp) = 1;
        else
            H_chi2(rp,amp) = 0;
        end
% 
%         %% permutation test 1
%         y_sin_me = mean(y_sin,2);
%         [pxx1,f] = pwelch(y_sin_me,ones(1,60),[],60,120);
%         pxx1 = pxx1(2:31);
%         f = f(2:31);
%         lambda = pxx1(6);
% 
%         lambda_rp=[];
%         for perm = 1:1000
%             y_sin_sim = Shuffle(y_sin);
%             y_sin_sim_me = mean(y_sin_sim,2);
%             [pxx_rp,f] = pwelch(y_sin_sim_me,ones(1,60),[],60,120);
%             pxx_rp = pxx_rp(2:31);
%             f = f(2:31);
%             lambda_rp(perm) = pxx_rp(6,:);
%         end
%         p = sum(lambda<lambda_rp)/1000;
%         if p<0.05
%             H_perm(rp,amp) = 1;
%         else
%             H_perm(rp,amp) = 0;
%         end

        %% permutation test 2
%         [pxx2,f] = pwelch(y_sin,ones(1,60),[],60,120);
%         pxx2 = pxx2(2:31,:);
%         f = f(2:31);
%         lambda = mean(pxx2(6,:));
% 
%         lambda_rp=[];
%         for perm = 1:1000
%             y_sin_sim = Shuffle(y_sin);
%             [pxx_rp,f] = pwelch(y_sin_sim,ones(1,60),[],60,120);
%             pxx_rp = pxx_rp(2:31,:);
%             f = f(2:31);
%             lambda_rp(perm) = mean(pxx_rp(6,:));
%         end
%         p = sum(lambda<lambda_rp)/1000;
%         if p<0.05
%             H_perm2(rp,amp) = 1;
%         else
%             H_perm2(rp,amp) = 0;
%         end

    end
end

power_f = mean(H_f,1);
power_chi2 = mean(H_chi2,1);
power_perm = mean(H_perm,1);
power_perm2 = mean(H_perm2,1);

acf = resp_sm;
back_delay_frames = [36:95];
Fs = 120; %Hz
T = 1/Fs; %period
L = size(acf,1);
hanningWin = hanning(L,'periodic'); % use periodic method as suggested by MATLAB for further fft analysis
n =256;