%% Generating a Switching Signal

function swi_t = Swi_signal(swi_seq,t)
% Input arguments: switching time sequence, time
% Output arguments: switching signal

swi_t=t;
swi_temp=zeros(1,length(swi_seq)); 
for r=1:length(swi_seq)
    swi_temp(r)=find(abs(swi_t-swi_seq(r))<1e-4);
end
swi_t(swi_temp)=0; % set the switching time to 0
temp=0;
for s=1:length(t)
    if(swi_t(s)==0)
        temp=temp+1;
        if(mod(temp,2)==1)
            swi_t(s)=1; % initial subsystem
        else
            swi_t(s)=2;
        end
    else
        swi_t(s)=swi_t(s-1);
    end
end
swi_t=[t;swi_t];