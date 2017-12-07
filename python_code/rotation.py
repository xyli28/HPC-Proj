import math 
import numpy as np

def random_quaternion(rand=None):
    """Retuen uniform random unit quaternion.

    rand: array like or None
	Three independent random  variables uniformally distributed 
	between 0 and 1.
    """
    if rand is None:
        rand = np.random.rand(3)
    else:
        assert len(rand) == 3
    r1 = math.sqrt(1.0 - rand[0])
    r2 = math.sqrt(rand[0])
    pi2 = math.pi * 2.0
    t1 = pi2 * rand[1]
    t2 = pi2 * rand[2]
    return np.array([np.sin(t1)*r1,np.cos(t1)*r1,
		     np.sin(t2)*r2,np.cos(t2)*r2])

def quaternion_matrix(quaternion):
    """Return rotation from quaternion."""

    q = np.array(quaternion, dtype=np.float64, copy=True)  
    n = np.dot(q, q)  
    if n < np.finfo(float).eps:
        return np.identity(3)
    q *= math.sqrt(2.0/n)
    q = np.outer(q,q)
    return np.array([
        [1.0-q[2,2]-q[3,3],    q[1,2]-q[0,3],    q[1,3]+q[0,2]],
	[    q[1,2]+q[0,3],1.0-q[1,1]-q[3,3],    q[2,3]-q[0,1]],
	[    q[1,3]-q[0,2],    q[2,3]+q[0,1],1.0-q[1,1]-q[2,2]]])

def random_rotation_matrix(rand=None):
    """Return uniform random rotation matrix.

    rand: array like
        Three independent random variables that are uniformly distributed
        between 0 and 1 for each returned quaternion.
    """
    return quaternion_matrix(random_quaternion(rand))

# Calculate rotation matrix using Rodrigues' rotation formula
# from axis and angle
def get_rotation_matrix(axis, angle):
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)

def get_random_R():
    # Original vector that we want to rotate
    orig = np.array([0., 0., 1.])
    
    # First we rotate around the original vector:
    angle1 = np.random.random() * 2 * np.pi
    
    # Get the first rotation matrix around original axis
    R1 = get_rotation_matrix(orig, angle1)
    
    # Then, find new orientation for the original axis by picking a random
    # point on the sphere
    vec = np.random.normal(size=3)
    vec /= np.linalg.norm(vec)
    
    # Get angle by taking the dot product between original vector and
    # new vector
    angle2 = np.arccos(np.dot(orig, vec)) # or np.arccos(vec[2])
    # Get axis for rotation by taking dot product between original vector
    # and new vector
    axis = np.cross(orig, vec)
    axisnorm = np.linalg.norm(axis)
    
    # If the outer product is too small, then that means the axis is aligned
    # with the original vector -- no rotation
    if np.abs(axisnorm) <= 1e-10:
        axis = np.sign(np.dot(orig, axis)) * orig
    else:
        axis /= axisnorm
    
    # Get the second rotation matrix
    R2 = get_rotation_matrix(axis, angle2)
    
    # Overall rotation matrix is composition of first two rotation matrices
    R = np.dot(R1, R2)
    return R

