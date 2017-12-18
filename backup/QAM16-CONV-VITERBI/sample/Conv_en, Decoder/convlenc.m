%Convolutional Encoder ; input=1 bit -> output=2 bits with 3 memory elements, Code Rate=1/2
function [encoded_sequence]=convlenc(message)


%TEST MESSAGES
% message=[1 0 1 0 1 1 1 0 0 0 1 1 0 1 1 0 0];%prb 0-1
% message=[0 0 1 0 1 0 1 0 1 0 0 1 1 0 1 0 0];
% message=[1 1 1 0 1 0 1 1 0 1 0 0 1 0 1 0 0];%prb 0-1
% message=[0 0 0 0 1 0 0 1 0 1 0 1 0 1 0 0 0];
% message=[1 0 1 0 1 0 1 0 0 1 0 1 0 1 1 0 0];%prb 0-1
% message=[0 1 1 0 1 0 1 0 1 1 1 1 0 0 1 0 0];
% message=[1 0 1 1 1 0 0 0 1 0 1 0 1 0 0 0 0];%prb 0-1
% message=[0 1 0 1 0 0 1 1 0 1 1 0 0 1 1 0 0];
% message=[1 0 1 1 0 1 0 1 0 0 1 1 0 1 0 0 0];%prb 0-1
% message=[1 0 1 1 0 1 0 1 0 0 1 0 1 1 0 0 0];%prb 0-1
% message=[0 0 1 0 0 1 1 1 0 0 1 0 1 0 1 0 0];
% message=[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0];%prb 0-1
% message=[1 0 1 1 0 0 1 1 0 0 1 1 0 0 1 0 0];%prb 0-1
% message=[0 1 1 1 0 0 1 1 0 0 1 1 0 0 1 0 0];
% message=[1 0 0 0 0 0 1 1 1 0 1 1 1 1 1 0 0];%prb 0-1

enco_mem=[0 0 0];   %# of memory elements=3
encoded_sequence=zeros(1,(length(message))*2);
       
       enco_mem(1,3)=enco_mem(1,2);
       enco_mem(1,2)=enco_mem(1,1);
       enco_mem(1,1)=message(1,1);

       temp=xor(enco_mem(1),enco_mem(2));      
       o1=xor(temp,enco_mem(3));             %generator polynomial=111
       o2=xor(enco_mem(1),enco_mem(3));      %generator polynomial=101
       encoded_sequence(1,1)=o1;
       encoded_sequence(1,2)=o2;

msg_len=length(message);
c=3;
for i=2:msg_len
         
       enco_mem(1,3)=enco_mem(1,2);
       enco_mem(1,2)=enco_mem(1,1);
       if(i<=msg_len)
       enco_mem(1,1)=message(1,i);
       else
       enco_mem(1,1)=0;
       end
              
       temp=xor(enco_mem(1),enco_mem(2));    
       o1=xor(temp,enco_mem(3));
       o2=xor(enco_mem(1),enco_mem(3));
       
       encoded_sequence(1,c)=o1;    %o1 generating polynomial(1,1,1)
       c=c+1;
       encoded_sequence(1,c)=o2;    %o2 generating polynomial(1,0,1)
       c=c+1;
end

