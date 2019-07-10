function strf=sub_revCOR_AM(spdata,stim,tlag)

stim=stim';
% tlag=100;%25;
seg=size(stim,1);     % samples
Fchan=size(stim,2);   % frequency channel
N=size(stim,2)/Fchan;
stimN=1;

strf=[];
for cnt0=1:N
    strf0=0;
    stim0=stim(:,(1:Fchan)+Fchan*(cnt0-1));
    for cnt=1:stimN
        stim1=stim0((1:seg)+seg*(cnt-1),:);
        dsum=spdata((1:seg)+seg*(cnt-1));
        dsum=dsum(:)-mean(dsum);   %new2
        for abc = 1:Fchan
            stimrow = stim1(:,abc);
            stimrow=stimrow-mean(stimrow); %new2
            %strftemp(abc,:)=xcorr(stimrow,dsum,25);
            strftemp(abc,:)=real(ifft(conj(fft(stimrow)).*fft(dsum)));
            %strftemp(abc,:)=strftemp(abc,:)/sqrt(mean(stimrow.^2));
        end
%         strftemp = strftemp*(2/sqrt(mean(mean(stim1.^2)))/stimN);
%         strftemp=strftemp/(seg);
        strf0=strf0+strftemp;
    end
    strf0 = strf0*(2/sqrt(mean(mean(stim0.^2)))/stimN);
    strf0=strf0/(seg*stimN);
    %strf0=fliplr(strf0(:,1:tlag));
    strf0=strf0(:,1:tlag);
    strf0=strf0-mean(strf0(:));
    strf=cat(1,strf,strf0);
end