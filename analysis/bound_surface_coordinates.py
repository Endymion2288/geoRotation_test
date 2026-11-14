"""Coordinate diagnostics for recorded ACTS bounds, without track transport."""
import numpy as np


ROTATE_F_TO_A = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])


def surface_position(bound, transform, surface_type):
    """Return native global position and d(position)/d(loc0,loc1).

    ACTS SurfaceType: Disc=2, Plane=4, Curvilinear=6. The caller must supply
    native type metadata; a 4x4 transform alone cannot specify local geometry.
    """
    transform = np.asarray(transform)
    if surface_type in (4, 6):
        local = np.r_[bound[:2], 0.]
        derivative = np.array([[1., 0.], [0., 1.], [0., 0.]])
    elif surface_type == 2:
        radius, phi = bound[:2]
        cp, sp = np.cos(phi), np.sin(phi)
        local = np.array([radius*cp, radius*sp, 0.])
        derivative = np.array([[cp, -radius*sp], [sp, radius*cp], [0., 0.]])
    else:
        raise ValueError('Position conversion requires supported native surface type: '+str(surface_type))
    return transform[:3, :3] @ local+transform[:3, 3], transform[:3, :3] @ derivative


def direction_f(bound, variant):
    """Global unit direction in FASER frame and its bound-coordinate Jacobian."""
    phi, theta = bound[2:4]
    cp, sp, ct, st = np.cos(phi), np.sin(phi), np.cos(theta), np.sin(theta)
    rotation = np.eye(3) if variant == 'beforeGeo' else ROTATE_F_TO_A.T
    derivative = np.zeros((3, 6))
    derivative[:, 2:4] = rotation @ np.array([[-sp*st, cp*ct], [cp*st, sp*ct], [0., -st]])
    return rotation @ np.array([cp*st, sp*st, ct]), derivative


def rotate_tangent(bound, covariance, source_surface, target_surface):
    """Map a tangent Gaussian under a proper global rotation of one surface.

    The full local bases are checked first. Local coordinates remain identical
    for any native surface type; only global spherical angles are transformed.
    This is a coordinate Jacobian, not finite-distribution equivalence.
    """
    rotation = np.eye(4)
    rotation[:3, :3] = ROTATE_F_TO_A
    np.testing.assert_allclose(rotation @ source_surface, target_surface, rtol=0, atol=2e-9)
    result = np.asarray(bound, dtype=float).copy()
    phi, theta = result[2:4]
    cp, sp, ct, st = np.cos(phi), np.sin(phi), np.cos(theta), np.sin(theta)
    direction = ROTATE_F_TO_A @ np.array([cp*st, sp*st, ct])
    x, y, z = direction
    radius = np.hypot(x, y)
    if radius == 0.:
        raise ValueError('Rotated angular chart is singular')
    result[2:4] = np.arctan2(y, x), np.arctan2(radius, z)
    derivative = ROTATE_F_TO_A @ np.array([[-sp*st, cp*ct], [cp*st, sp*ct], [0., -st]])
    norm2 = float(direction @ direction)
    gradient = np.array([[-y/radius**2, x/radius**2, 0.],
        [z*x/(radius*norm2), z*y/(radius*norm2), -radius/norm2]])
    jacobian = np.eye(6)
    jacobian[2:4, 2:4] = gradient @ derivative
    return result, jacobian @ covariance @ jacobian.T
