function [t] = fitterms(pred_data, pbest, feats)

t = table;
D = size(feats,2);

for i=1:D
    t.( strcat(feats{i}, '_GAM') ) = sigmoid(pred_data(:,i), pbest(1+(i-1)*2 : i*2));
    t.( feats{i} ) = pred_data(:,i);
end

end