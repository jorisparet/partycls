import numpy

def cartesian_to_spherical(r_xyz):
    """
    Cartesian to spherical coordinates conversion.
    """
    r_sph = numpy.empty(r_xyz.shape)
    xy = r_xyz[:,0]**2 + r_xyz[:,1]**2
    r_sph[:,0] = numpy.sqrt(xy + r_xyz[:,2]**2) # module
    r_sph[:,1] = numpy.arctan2(r_xyz[:,1], r_xyz[:,0]) # longitutde
    r_sph[:,2] = numpy.arctan2(numpy.sqrt(xy), r_xyz[:,2]) # latitude
    return r_sph

def pbc(r, box):
    """
    Apply periodic boundary conditions to vector `r`.
    """
    for i in range(len(box)):
        r[...,i] = numpy.where(r[...,i] >  box[i]/2, r[...,i] - box[i], r[...,i])
        r[...,i] = numpy.where(r[...,i] < -box[i]/2, r[...,i] + box[i], r[...,i])
    return r