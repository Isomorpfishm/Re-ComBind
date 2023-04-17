from utils import *



"""
This section computes scores for hydrogen bond (PNAS)
"""

def fdist(x, dist_cut, dist_bound):
    if x > dist_cut:
        return 0
    elif (dist_bound < x) and (x <= dist_cut):
        return (dist_cut-x)/(dist_cut-dist_bound)
    else:
        return 1
    
def fangle(x, angle_cut, angle_bound):
    if (angle_bound <= x) and (x < angle_cut):
        return (x-angle_bound)/(angle_cut-angle_bound)
    elif x < angle_bound:
        return 0
    else:
        return 1

def _piecewise_hbond(data_dist, data_angle, settings):
    dist_cut = settings['hbond_dist_cut']
    angle_cut = settings['hbond_angle_cut']
    dist_bound = settings['hbond_dist_bound']
    angle_bound = settings['hbond_angle_bound']

    hbond_score_dist, hbond_score_angle = [], []

    for i in range(len(data_dist)):
        hbond_score_dist.append(fdist(data_dist[i], dist_cut, dist_bound))
        hbond_score_angle.append(fangle(data_angle[i], angle_cut, angle_bound))
 
    return hbond_score_dist, hbond_score_angle
    


"""
This section computes scores for hydrophobic contacts (PNAS)
"""

def gdist(x, dist_cut, dist_bound, vdw):
    if (dist_bound*vdw < x) and (x <= dist_cut*vdw):
        return (dist_cut*vdw-x)/((dist_cut-dist_bound)*vdw)
    elif x <= dist_bound*vdw:
        return 1
    else:
        return 0

def _piecewise_hphob(data_dist, settings, vdw=3.4):
    dist_cut = settings['hphob_dist_cut']
    dist_bound = settings['hphob_dist_bound']

    hphob_score = []

    for i in range(len(data_dist)):
        hphob_score.append(gdist(data_dist[i], dist_cut, dist_bound, vdw))
 
    return hphob_score
    


"""
This section computes scores for salt bridges (PNAS)
"""

def hdist(x, dist_cut, dist_bound):
    if x > dist_cut:
        return 0
    elif (dist_bound < x) and (x <= dist_cut):
        return (dist_cut-x)/(dist_cut-dist_bound)
    else:
        return 1

def _piecewise_sbridge(data_dist, settings):
    dist_cut = settings['sbridge_dist_cut']
    dist_bound = settings['sbridge_dist_bound']

    sbridge_score = []

    for i in range(len(data_dist)):
        sbridge_score.append(hdist(data_dist[i], dist_cut, dist_bound))
 
    return sbridge_score
    


"""
This section computes scores for pi-pi stacking
"""

def idist(x, cut, bound):
    if x > cut:
        return 0
    elif (bound < x) and (x <= cut):
        return (cut-x)/(cut-bound)
    else:
        return 1
    
def iangle(dist, angle, dist_cut, dist_bound, angle_cut, angle_bound):
    if angle > 90.0:
        angle = 180.0 - angle
    else:
        pass

    if (angle_bound <= angle) and (angle <= angle_cut):
        if dist >= 5.0:
            return (dist_cut - dist)/(dist_cut - 5.0)
        else:
            return dist/dist_cut
    elif (0 <= angle) and (angle < 30.0):
        if dist <= dist_bound:
            return dist/dist_bound
        else:
            return (dist_cut - dist)/(dist_cut - dist_bound)
    else:
        return (dist_cut-dist)/dist_cut


def _piecewise_pistack(data_dist, data_angle, data_angle_1, data_angle_2, settings):
    dist_cut = settings['pistack_dist_cut']
    dist_bound = settings['pistack_dist_bound']
    angle_bound = settings['pistack_angle_bound']
    angle_cut = settings['pistack_angle_cut']

    pistack_score_dist, pistack_score_angle = [], []

    for i in range(len(data_dist)):
        pistack_score_dist.append(idist(x=data_dist[i], cut=dist_cut, bound=dist_bound))
        pistack_score_angle.append(iangle(dist=data_dist[i], angle=data_angle[i], dist_cut=dist_cut, dist_bound = dist_bound, angle_cut=angle_cut, angle_bound=angle_bound))
 
    return pistack_score_dist, pistack_score_angle



"""
This section computes scores for pi-cation interaction
"""

# def jdist(x, cut, bound):
#     if x <= bound:
#         return x/bound
#     elif (bound < x) and (x <= cut):
#         return 1
#     else:
#         return (6.5-x)/(6.5-cut)

def jdist(x, cut, bound):
    if x <= bound:
        return 1
    else:
        return (6.5-x)/(6.5-bound)
    
def jangle(x, cut, bound):
    if x > 90.0:
        return (180.0-x)/90.0
    else:
        return x/90.0


def _piecewise_pication(data_dist, data_angle, settings):
    dist_cut = settings['pication_dist_cut']
    dist_bound = settings['pication_dist_bound']
    angle_bound = settings['pication_angle_bound']
    angle_cut = settings['pication_angle_cut']

    pication_score_dist, pication_score_angle = [], []

    for i in range(len(data_dist)):
        pication_score_dist.append(jdist(data_dist[i], dist_cut, dist_bound))
        pication_score_angle.append(jangle(data_angle[i], angle_cut, angle_bound))
 
    return pication_score_dist, pication_score_angle
    
    
    
"""
This section computes scores for halogen bond
"""

# def kdist(x, dist_cut, dist_bound):
#     if x > dist_cut:
#         return 0
#     elif (dist_bound < x) and (x <= dist_bound+0.3):
#         return 1
#     elif (dist_bound+0.3 < x) and (x <= dist_cut):
#         return (dist_cut-x)/(dist_cut-dist_bound-0.3)
#     else:
#         return 1

def kdist(x, dist_cut, dist_bound):
    if x < dist_bound:
        return 1
    else:
        return (dist_cut-x)/(dist_bound)

def kangle(x, y, angle_cut, angle_bound):
    if (90.0 <= y) and (y < angle_bound):
        s2 = (angle_bound-y)/(angle_bound-90.0)
    elif angle_bound <= y:
        s2 = 1
    else:
        s2 = 0

    return s2

def _piecewise_Xbond(data_dist, data_angle_haan, data_angle_ahhn, settings):
    dist_cut = settings['Xbond_dist_cut']
    angle_cut = settings['Xbond_angle_cut']
    dist_bound = settings['Xbond_dist_bound']
    angle_bound = settings['Xbond_angle_bound']

    Xbond_score_dist, Xbond_score_angle = [], []

    for i in range(len(data_dist)):
        Xbond_score_dist.append(kdist(data_dist[i], dist_cut, dist_bound))
        Xbond_score_angle.append(kangle(data_angle_haan[i], data_angle_ahhn[i], angle_cut, angle_bound))
 
    return Xbond_score_dist, Xbond_score_angle


