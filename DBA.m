%*******************************************************************************
 % Copyright (C) 2013 Francois PETITJEAN, Ioannis PAPARRIZOS
 % This program is free software: you can redistribute it and/or modify
 % it under the terms of the GNU General Public License as published by
 % the Free Software Foundation, version 3 of the License.
 % 
 % This program is distributed in the hope that it will be useful,
 % but WITHOUT ANY WARRANTY; without even the implied warranty of
 % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 % GNU General Public License for more details.well
 % 
 % You should have received a copy of the GNU General Public License
 % along with this program.  If not, see <http://www.gnu.org/licenses/>.
 %*****************************************************************************/ 

% function average = DBA(sequences)
% 	index=randi(length(sequences),1);
% 	average=repmat(sequences{index},1);
% 	for i=1:15
% 		average=DBA_one_iteration(average,sequences);
% 	end
% end

function average = DBA(sequences)
    %sequences is a cell array of time series
    average = repmat(sequences{medoidIndex(sequences)},1);
	for i=1:15
		average=DBA_one_iteration(average,sequences);
	end
end

function sos = sumOfSquares(s,sequences)
    sos = 0.0;
    for i=1:length(sequences)
        dist = dtw(s,sequences{i});
        sos = sos + dist * dist;
    end
end

function score = dtw(S,T)
    costM = zeros(size(S, 1),size(T, 1));
    
    costM(1,1) = sum((S(1, :) - T(1, :)).^2);
    for i=2:size(S, 1)
        costM(i,1)= costM(i-1,1)+ sum((S(i,:)-T(1,:)).^2);
    end
    for i=2:size(T, 1)
        costM(1,i)= costM(1,i-1)+ sum((S(1,:)-T(i,:)).^2);
    end
    for i=2:size(S, 1)
        for j=2:size(T, 1)
            costM(i,j)=min(min(costM(i-1,j-1),costM(i,j-1)),costM(i-1,j))+sum((S(i,:)-T(j,:)).^2);
        end
    end
    score = sqrt(costM(size(S,1),size(T,1)));
end

function index = medoidIndex(sequences) 
	index = -1;
    lowestInertia = Inf;
	for i=1:length(sequences)
        tmpInertia = sumOfSquares(sequences{i},sequences);
        if (tmpInertia < lowestInertia)
            index = i;
            lowestInertia = tmpInertia;
        end
    end
end


function average = DBA_one_iteration(averageS,sequences)

	tupleAssociation = cell (1, size(averageS,1));
	for t=1:size(averageS,1)
		tupleAssociation{t}=[];
	end

	costMatrix = zeros(1000,1000);
	pathMatrix = zeros(1000,1000);

	for k=1:length(sequences)
	    sequence = sequences{k};
	    costMatrix(1,1) = distanceTo(averageS(1,:),sequence(1,:));
	    pathMatrix(1,1) = -1;
	    for i=2:size(averageS,1)
		costMatrix(i,1) = costMatrix(i-1,1) + distanceTo(averageS(i),sequence(1));
		pathMatrix(i,1) = 2;
	    end
	    
	    for j=2:size(sequence,1)
		costMatrix(1,j) = costMatrix(1,j-1) + distanceTo(sequence(j,:),averageS(1,:));
		pathMatrix(1,j) = 1;
	    end
	    
	    for i=2:size(averageS,1)
		for j=2:size(sequence,1)
		    indiceRes = ArgMin3(costMatrix(i-1,j-1),costMatrix(i,j-1),costMatrix(i-1,j));
		    pathMatrix(i,j)=indiceRes;
		    
		    if indiceRes==0
			res = costMatrix(i-1,j-1);
		    elseif indiceRes==1
			res = costMatrix(i,j-1);
		    elseif indiceRes==2
			res = costMatrix(i-1,j);
		    end
		    
		    costMatrix(i,j) = res + distanceTo(averageS(i,:),sequence(j,:));
		    
		end
	    end
	    
	    i=size(averageS,1);
	    j=size(sequence,1);
	    
	    while(true)
		tupleAssociation{i}(end+1,:) = sequence(j,:);
		if pathMatrix(i,j)==0
		    i=i-1;
		    j=j-1;
		elseif pathMatrix(i,j)==1
		    j=j-1;
		elseif pathMatrix(i,j)==2
		    i=i-1;          
		else
		    break;
		end
	    end
	    
	end

	for t=1:size(averageS,1)
	   averageS(t, :) = mean(tupleAssociation{t}, 1);
	end
	   
	average = averageS;

end

function value = ArgMin3(a,b,c)

	if (a<b)
	    if (a<c)
		value=0;
		return;
	    else
		value=2;
		return;
	    end
	else
	    if (b<c)
		value=1;
		return;
	    else
		value=2;
		return;
	    end
	end

end


function dist = distanceTo(a,b)
    dist=(a-b)*(a-b)';
end

function ex = test()
    sequences = {};
    sequences{100}=[];
    for i=1:100
        length = randi(100);
        sequences{i}=rand(1,length);
    end
    mean=DBA(sequences);
end
 

	    
