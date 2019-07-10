% label pixels in a 3d matrix of size T*N*N by 26-connectivity 
T = 10;
N = 1024;
data = zeros(T,N,N); % predefine the array
for i = 1:T % load the data
    load("blk_"+string(i)+".mat");
    data(i,:,:) = re;
end
% ==============================================
% To deal with the periodic boundary and the edge elements,
% we simply add two layers: 1. 1 row and 1 column for periodicity
% 2. 2 rows, 2 columns and 1 T-dimensional layer margins for the simplicity
% of code
ndata = zeros(T+1,N+3,N+3); % new data matrix
nlabel = zeros(T+1,N+3,N+3); % new label matrix
eqdict = zeros(2,N*N); % dictionary that stores the equivalence of labels
f9 = zeros(3,3);%  the 9 pixels in 3x3 matrix that is part of 13 surrounding
f3 = zeros(1,3);%  the 3 pixels in 1x3 row vector that is part of 13

ndata(2:end,3:N+2,3:N+2) = data;
ndata(2:end,2,3:N+2) = data(:,N,:);
ndata(2:end,3:N+2,2) = data(:,:,N);
ndata(2:end,2,2) = data(:,N,N);

labelno = 0;

% the first pass
for i = 2:T
    for j = 2:N+2
        for k = 2:N+2
            if ndata(i,j,k) == 1   % we only label the pixel with value 1
                f9 = reshape(nlabel(i-1,j-1:j+1,k-1:k+1),[1,9]); % check the 13 surroundings 
                f3 = reshape(nlabel(i,j-1,k-1:k+1),[1,3]);
                f1 = nlabel(i,j,k-1);
                eqlabel = unique([f9(f9~=0),f3(f3~=0),f1(f1~=0)]);% store all the surrounding labels 
                if sum(eqlabel) ~= 0 % if there is any surrounding label
                    mini = min(eqdict(2,eqlabel)); % assign the label with smallest index
                    nlabel(i,j,k) = mini;
                    for m = 1:length(eqlabel)
                        eqdict(2,eqlabel(m)) = mini; % store the equivalence of labels
                    end
                    
                else             % if no surrounding label, create new label in the dictionary
                    labelno = labelno + 1; 
                    nlabel(i,j,k) = labelno;
                    eqdict(:,labelno) = [labelno;labelno];
                end
            end
        end
        if ndata(i,j,N+2) ~=0   % keep the consistence of the edge pixels for periodicity
            maxo2 = max(nlabel(i,j,2),nlabel(i,j,N+2));
            mino2 = min(nlabel(i,j,2),nlabel(i,j,N+2));
            for u = 1:labelno
                if eqdict(2,u) == maxo2
                    eqdict(2,u) = mino2;
                end
            end
        end
    end
    for v = 2:N+2      % keep the consistence of the edge pixels for periodicity
        if ndata(i,N+2,v) ~=0
            maxo2 = max(nlabel(i,N+2,v),nlabel(i,2,v));
            mino2 = min(nlabel(i,N+2,v),nlabel(i,2,v));
            for u = 1:labelno
                if eqdict(2,u) == maxo2
                    eqdict(2,u) = mino2;
                end
            end
        end
    end
end

% reassign the index of labels to make sure there is no empty between
% indices
eqdict = reshape(eqdict(eqdict~=0),2,[]);
row1 = unique(eqdict(2,:));
row2 = 1:length(row1);
for i = 1:labelno
    eqdict(2,i) = row2(row1==eqdict(2,i));
end

% the second pass
for i = 2:T
    for j = 2:N+2
        for k = 2:N+2
            if ndata(i,j,k) == 1
                nlabel(i,j,k) = eqdict(2,nlabel(i,j,k));
            end
        end
    end
end

label = nlabel(2:end,3:N+2,3:N+2);

%==========================================
% below is an example of how to select certain label of group
groupno = 3;
sele = reshape(label == groupno,T,N,N);
aviobj=VideoWriter('example_group.avi');
open(aviobj);

for a = 1:10
slice = reshape(sele(a,:,:),N,N);
imagesc([0.5:N-0.5],[0.5:N-0.5],slice);
colormap([1 1 1;0 0 0]);
set(gca,'YDir','normal');
grid minor
F(a) = getframe;
for i = 1:30
writeVideo(aviobj,F(a));
end
end
movie(F,1,2)
close(aviobj);
