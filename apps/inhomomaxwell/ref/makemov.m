function makemov(dir, meqn, s, e, ymin, ymax)

# dir is name of directory: same as run_name
# load X coordinates
fn = sprintf('%s/frame.x', dir);
X = load(fn);
xmin = X(1);
xmax = X(end);

# set axis for all plots
axis([xmin,xmax, ymin,ymax]);

# make plots one by one
for p = s:e
  # load data
  fn = sprintf('%s/frame.q%d', dir, p);
  fd = load(fn);
  # do something with data here
  q = fd(:,meqn);  
  # make plot
  plot(X,q); xlabel('X');
  # save to file
  fn = sprintf('%s/plot%04d.png', dir, p);
  print(fn,'-dpng','-color');
end

axis;
