function [Sset,SetLb,SetUb,ReactionSet,MetaboliteSet] = getMinRxnSet(S,lb,ub)
% Computes a minimum reaction set with updated bounds
%
% Inputs:  stoichiometric matrix (S), lower (lb) and upper (ub) flux bounds
% Outputs: minimal reaction set and bounds
%
%%%%%%%%%%%%%%%%%%%%%% Lars Nielsen UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = full(S);
[NoMetabolites,NoReactions] = size(S);
ReactionSet   = zeros(NoReactions,4);
MetaboliteSet = zeros(NoMetabolites,1);
K = null(S);
for i=1:NoReactions
    if isempty(find(K(i,:)))
        ReactionSet(i,1)=-1;
        ReactionSet(i,2)= 0;
        ReactionSet(i,3)= 0;
        ReactionSet(i,4)= 0;            
    else
        [C,idx]=max(abs(K(i,:)));
        ReactionSet(i,2)=K(i,idx);
        K(i,:)=K(i,:)./K(i,idx);
    end
end
ReactionSetNo = 1;
for i=1:NoReactions
    if ~ReactionSet(i,1)
        ReactionSet(i,1) = ReactionSetNo;
        CurrentSetIdx = i;
        for j=i+1:NoReactions
            if min(norm(K(i,:)-K(j,:)),norm(K(i,:)+K(j,:)))<1e-10
                ReactionSet(j,1) = ReactionSetNo;
                CurrentSetIdx = [CurrentSetIdx,j];
            end
        end
        CurrentSet = ReactionSet(CurrentSetIdx,:);
        [MinAbsWeight,idx] = min(abs(CurrentSet(:,2)));
        ReactionSet(CurrentSetIdx,2) = CurrentSet(:,2)/CurrentSet(idx,2);        
        CurrentLb = -1e9;
        CurrentUb = 1e9;
        for j = CurrentSetIdx
            if ReactionSet(j,2)>0
                CurrentLb = max(CurrentLb,lb(j)/ReactionSet(j,2));
                CurrentUb = min(CurrentUb,ub(j)/ReactionSet(j,2));
            else
                CurrentLb = max(CurrentLb,ub(j)/ReactionSet(j,2));
                CurrentUb = min(CurrentUb,lb(j)/ReactionSet(j,2));
            end
        end
        if abs(CurrentUb-CurrentLb)<1e-9
            ReactionSet(CurrentSetIdx,1)=-2;
            ReactionSet(CurrentSetIdx,3)=CurrentLb;
            ReactionSet(CurrentSetIdx,4)=CurrentUb;            
        elseif abs(CurrentUb) <1e-9 % Change direction
            ReactionSet(CurrentSetIdx,2) = -ReactionSet(CurrentSetIdx,2);
            ReactionSet(CurrentSetIdx,3)= 0;
            ReactionSet(CurrentSetIdx,4)=-CurrentLb;            
            ReactionSetNo = ReactionSetNo + 1;
        else
            ReactionSet(CurrentSetIdx,3)=CurrentLb;
            ReactionSet(CurrentSetIdx,4)=CurrentUb;            
            ReactionSetNo = ReactionSetNo + 1;
        end
    end
end
ReactionSetNo = max(ReactionSet(:,1));
Sset = zeros(NoMetabolites,ReactionSetNo);
SetLb = zeros(ReactionSetNo,1);
SetUb = zeros(ReactionSetNo,1);
for col = 1:ReactionSetNo
    CurrentSetIdx = find(ReactionSet(:,1)==col);
    for j= CurrentSetIdx
        Sset(:,col) = Sset(:,col) + S(:,j)*ReactionSet(j,2);
    end
    SetLb(col) = ReactionSet(CurrentSetIdx(1),3);
    SetUb(col) = ReactionSet(CurrentSetIdx(1),4);
end
Sset(find(abs(Sset)<1e-9))=0;
CurrentSetMetaboliteNo = 1;
for row = 1:NoMetabolites
    if isempty(find(Sset(row,:)))
        MetaboliteSet(row) = -1;
    else
        MetaboliteSet(row) = CurrentSetMetaboliteNo;
        CurrentSetMetaboliteNo = CurrentSetMetaboliteNo+1;
    end
end
Sset = Sset(find(MetaboliteSet>0),:);        
    