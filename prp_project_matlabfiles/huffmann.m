function codebook = huffmann(symprobvec,lensymprobvec)

% initialising
codebook = zeros(lensymprobvec,3);
%special case if d=1,give 0,0,0
if lensymprobvec==1
    return;
end

%codebook(k,1) should store symbols
codebook(:,1) = symprobvec(:,1);

% sum of 2 childs at the end is given to parent[55555][555105]-[55510]
symprobvec(lensymprobvec-1,2) = symprobvec(lensymprobvec-1,2) + symprobvec(lensymprobvec,2);

% defining some wierd value thats not asci and decreasing length to n-1
% from n
symprobvec(lensymprobvec-1,1) = -1;

lensymprobvec = lensymprobvec -1;
symprobvec = symprobvec(1:lensymprobvec,:);

% sorting from to high to low[11,10,555]
symprobvec = sortrows(symprobvec,2,'descend');
 
% finding at which length the code broke actually we know it as d-1 in
% previous but we dont know after sorted,k=1

for i=1:lensymprobvec
    if symprobvec(i,1)==-1
        k=i;
        symprobvec(i,1)=128;
        break;
    end
end

% calling recursive function
codebook1 = huffmann(symprobvec,lensymprobvec);

% denoting the rows using recursived one [10,5,5,5];n
codebook(lensymprobvec,3)=codebook1(k,3)+1;
codebook(lensymprobvec+1,3)=codebook1(k,3)+1;
codebook(lensymprobvec,2)=1+2*codebook1(k,2);
codebook(lensymprobvec+1,2)=2*codebook1(k,2);

for ii=1:lensymprobvec-1
    if ii<k
        codebook(ii,2)=codebook1(ii,2);
        codebook(ii,3)=codebook1(ii,3);
    end
    % after sorting we mixed the p_n,p_n-1 prob to one and we sorted again
    % because of this we shouldnot consider its symbol so we will find its
    % position k and all index<k it will remain same for index>=k its the
    % after code rom codebook
    if ii>=k
        codebook(ii,2)=codebook1(ii+1,2);
        codebook(ii,3)=codebook1(ii+1,3);
    end
end


end

