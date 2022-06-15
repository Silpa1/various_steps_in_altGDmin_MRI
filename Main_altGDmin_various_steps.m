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
    X_mat=reshape(X_image,[n,q]);
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
            ksc = reshape( fft2( reshape(X_mat(:,k), [n1 n2]) ), [n,1]) ;
            Y(1:mk(k),k)=double(ksc(S2(1:mk(k),k)));
        end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Mean  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        [Xbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,10);
        Xhat=repmat(Xbar_hat,1,q);
        Xhat_bar=reshape(Xhat,[n1,n2,q]);
        Time_mean=toc;
        Error_mean=RMSE_modi(Xhat_bar,X_image);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   AltGDmin  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T=70;
        tic;
        [Uhat]=initAltGDMin(Y);
        [Uhat2, Bhat2]=GDMin_wi(T,Uhat,Y);
        xT=Uhat2*Bhat2;
        L(:,1:q)=xT;
        Xhat_GD=xT;
        Xhat_GD=reshape(Xhat_GD,[n1, n2,q]); 
        Time_GD=  toc;
        Error_GD=RMSE_modi(Xhat_GD,X_image);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Mean + AltgdMin  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T=70;
        tic;
        [Xbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,10);
        Ybar_hat=Afft(Xbar_hat);
        Ybar_hat=reshape(Ybar_hat,[m,q]);
        Yinter=Y-Ybar_hat;
        [Uhat]=initAltGDMin(Yinter);
        [Uhat2, Bhat2]=GDMin_wi(T,Uhat,Yinter);
        xT=Uhat2*Bhat2; 
        Xhat_GD_mean1=xT+Xbar_hat;
        Xhat_GD_mean=reshape(Xhat_GD_mean1,[n1, n2,q]);
        Time_GD_mean=  toc;
        Error_GD_mean=RMSE_modi(Xhat_GD_mean,X_image);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% altGDmin_MRI2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T=70;
        tic;
        [Xbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,10);
        Ybar_hat=Afft(Xbar_hat);
        Ybar_hat=reshape(Ybar_hat,[m,q]);
        Yinter=Y-Ybar_hat;
        [Uhat]=initAltGDMin(Yinter);
        [Uhat2, Bhat2]=GDMin_wi(T,Uhat,Yinter);
        xT=Uhat2*Bhat2;
        L(:,1:q)=xT+Xbar_hat;
        Ymec=Y-Afft(L);
        E_mec=[];
        for kk=1:1:q
            E_mec(:,kk)=cgls_modi(@Afft_modi,@At_modi, Ymec(:,kk) ,0,1e-36,3);
        end
        Xhat_GD_MEC1=L+E_mec;
        Xhat_GD_MEC=reshape(Xhat_GD_MEC1,[n1, n2,q]);
        
        Time_MRI=  toc;
        Error_MRI=RMSE_modi(Xhat_GD_MEC,X_image);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(fid, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial(ii),Error_mean,Time_mean,Error_GD,Time_GD,Error_GD_mean,Time_GD_mean,Error_MRI,Time_MRI);
    end
    
    
    
    
    
    %     p=[filename,radial;Error_GD_Sparse;Time_GD_Sparse];
    %     fid=fopen('AltGDMin_LSM.txt','w');
    %     fprintf(fid,'  %s(%s)     %s  \n','name','radial','AltGDMin');
    %     %     fprintf(fid,'  %s(%s)  &  %s  &  %s  &   %s  \n',name,'radial','k-t-SLR','L+S','AltGDMin');
    %     fprintf(fid,'  %s(%d)     & %8.4f (%2.4f) \n',p);
    %     fclose(fid);
end
fclose(fid);

function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end