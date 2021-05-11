function MAP = getMAP(posterior, avg, vals)

for d = avg
   dts = mean(posterior,d);
end
dts = squeeze(dts);

for t = 1:size(dts,2)
   msk = dts(:,t) == max(dts(:,t));
   MAP(t) = vals(msk);
end

end