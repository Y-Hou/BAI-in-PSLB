function [index,CAid] = LCA(X,w,b,CAid,Ainv_XYt)
%LCA Linear Context Alignment

for i=1:length(CAid)
    if LCD(X,w,b,Ainv_XYt,CAid{i}) == 0
        if Ainv_XYt(2)*CAid{i}(2)<0
            break;
        end
        CAid{i} = Ainv_XYt;
        index = i;
        return;
    end
end
index = length(CAid) + 1;
CAid{index} = Ainv_XYt;

end

