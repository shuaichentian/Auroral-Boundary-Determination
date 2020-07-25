function phi = SDF(init_a)
%-- converts a mask to a SDF
  phi=bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;
  phi = -double(phi); % 
end