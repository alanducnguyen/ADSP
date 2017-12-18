%Returns the hamming distance b/w two 2 bit codes
function [HD]=hamm_dist(A,B)
HD=0;
for i=1:2
    a=bitget(A, i);
    b=bitget(B, i);
    
    if(a~=b)
        HD=HD+1;
    end
end