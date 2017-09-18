%Mehreen Ali
%mehreen.ali@helsinki.fi

function CIndex2 = civalue(y,yhat)

 ci = zeros(1,size(y,2));

   for k= 1:size(y,2) % col
    s = 0; % score sum
    n = 0; %  eligible pair of scores
    for i= 1:size(y,1) % row
        for j= i+1:size(y,1) % next row
            if (i < j)
                if(y(i,k) > y(j,k))
                    s = s + (yhat(i,k) > yhat(j,k)) + 0.5 * (yhat(i,k) == yhat(j,k));
                    n = n + 1;
                elseif (y(i,k) < y(j,k))
                    s = s + (yhat(i,k) < yhat(j,k)) + 0.5 * (yhat(i,k) == yhat(j,k));
                    n = n + 1;
                end
            end
        end
    end
   ci(k) = s / n;
   end

CIndex2 = nanmean(ci);