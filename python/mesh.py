import distmesh as dm
fd = lambda p: dm.ddiff(dm.drectangle(p,-1,1,-1,1), dm.dcircle(p,0,0,0.5))
fh = lambda p: 0.05+0.3*dm.dcircle(p,0,0,0.5)
p, t = dm.distmesh2d(fd, fh, 0.05, (-1,-1,1,1), [(-1,-1),(-1,1),(1,-1),(1,1)])
