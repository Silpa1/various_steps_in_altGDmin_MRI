clc;
clear all;close all;
global n1 n2 n mk m S2 q kk;

%filenames=  {'Pincat.mat','brain_T2_T1.mat','speech_seq.mat','Cardiac_ocmr_data.mat','lowres_speech.mat','FB_ungated.mat'};
filenames={'FB_ungated.mat'};


[fid,msg] = fopen('Comparison.txt','wt');
fprintf(fid, '%s(%s) & %s & %s &  %s &  %s    \n','Dataset','Radial','Only Mean','AltGDMin','Mean+AltGDMin','AltGDMin MRI');
for jj = 1:1:numel(filenames)
    S = load(filenames{jj});
    X_image=cell2mat(struct2cell(S));
    [~,name,~] = fileparts(filenames{jj});
    radial=[4,8,16];
    [n1,n2,q]=size(X_image);
    n=n1*n2;
    X_star=reshape(X_image,[n,q]);
    for ii=1:1:length(radial)
        [mask]=goldencart(n1,n2,q,radial(ii));
        mask = fftshift(fftshift(mask,1),2);
        Samp_loc=double(find(logical(mask)));
        mask3=reshape(mask,[n1*n2, q]);
        mk=[];
        for i=1:1:q
            mk(i)=length(find(logical(mask3(:,i))));
            S2(1:mk(i),i)=double(find(logical(mask3(:,i))));
        end
        m=max(mk);
        Y=zeros(m,q);
        for k=1:1:q
            ksc = reshape( fft2( reshape(X_star(:,k), [n1 n2]) ), [n,1]) ;
            Y(1:mk(k),k)=double(ksc(S2(1:mk(k),k)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Mean  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        [zbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,10);
        Xhat=repmat(zbar_hat,1,q);
        Xhat_mean=reshape(Xhat,[n1,n2,q]);
        Time_mean=toc;
        Error_mean=RMSE_modi(Xhat_mean,X_image);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   AltGDmin  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m=max(mk);
        Y=Y(1:m,1:q);
        T=70;
        tic;
        [U0]=initAltGDMin(Y);
        [Uhat2, Bhat2]=AltGDmin(T,U0,Y);
        X_hat=reshape(Uhat2*Bhat2,[n1, n2,q]);
        Time_GD=  toc;
        Error_GD=RMSE_modi(X_hat,X_image);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Mean + AltgdMin  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T=70;
        tic;
        [zbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,10);
        Ytemp=reshape(Afft(zbar_hat),[m,q]);
        Ybar=Y-Ytemp;
        [U0]=initAltGDMin(Ybar);
        [Uhat2, Bhat2]=AltGDmin(T,U0,Ybar);
        X_hat=Uhat2*Bhat2;
        Xhat_GD_mean=X_hat+zbar_hat;
        Xhat_GD_mean=reshape(Xhat_GD_mean,[n1, n2,q]);
        Time_GD_mean=  toc;
        Error_GD_mean=RMSE_modi(Xhat_GD_mean,X_image);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% altGDmin_MRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m=max(mk);
        L=[];
        T=70;
        tic; 
        [zbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,10);
        Ytemp=reshape(Afft(zbar_hat),[m,q]);
        Ybar=Y-Ytemp;
        [U0]=initAltGDMin(Ybar);
        [Uhat2, Bhat2]=AltGDmin(T,U0,Ybar);
        X_hat=Uhat2*Bhat2;
        
        Yhat_hat=Y-Afft(X_hat+zbar_hat);
        Ehat=[];
        for kk=1:1:q
            Ehat(:,kk)=cgls_modi(@Afft_modi,@At_modi, Yhat_hat(:,kk) ,0,1e-36,3);
        end
        Zhat=X_hat+zbar_hat+Ehat;
        Zhat_MRI=reshape(Zhat,[n1, n2,q]);
        
        Time_MRI=  toc;
        Error_MRI=RMSE_modi(Zhat_MRI,X_image);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(fid, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial(ii),Error_mean,Time_mean,Error_GD,Time_GD,Error_GD_mean,Time_GD_mean,Error_MRI,Time_MRI);
    end
end
fclose(fid);

function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end
