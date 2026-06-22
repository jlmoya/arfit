
demopath = get_absolute_file_path("arfit.dem.gateway.sce");

subdemolist=['Introduction', 'ardem.sce'];


subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
