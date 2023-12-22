# Domain setup
global zmin
global zmax
global nelements
global h_stem
global h_leaf
global dz_bottom
global dz_top
# Domain setup


land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = (dz_bottom, dz_top),
)
canopy_domain = ClimaLSM.Domains.obtain_surface_domain(land_domain)

compartment_midpoints = [h_leaf / 2]#compartment_midpoints = [h_stem / 2, h_stem + h_leaf / 2];
compartment_surfaces = [zmax, h_leaf]
#compartment_surfaces = [zmax, h_stem, h_stem + h_leaf];



