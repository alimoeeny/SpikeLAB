function r = normalize(a)
  mx = max(a);
  mi = min(a);
  ra = mx - mi;
  for i = 1:length(a)
      a(i) = (a(i) - mi) / ra;
  end
  
  r = a;
  