# RVT code

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1 as mpl

# Parameters (non dimensionalized)
a = 1  # related to distance from midpoint to inflection point
b = 1  # height of hill (b is negative for valley, positive for ridge)
pr = 0.25  # Poisson's ratio
rg = 1  # density times gravity
ts = -1  # far field stress (negative is compression, positive is tension)

uinc = 0.01  # horizontal step distance
umin = -3  # horizontal starting position
umax = 3  # horizontal ending position

vinc = 0.01  # vertical step distance
vmin = -4  # vertical starting position
vmax = 0  # vertical ending position

if b < 0:
    if a <= abs(b):
        raise Exception("While b < 0, the condition a > abs(b) must be met.")

# ******************************************************************************
# Beginning of RVT routine
#
# ******************************************************************************

# small radius around point w = -ia
rad = 1e-4

so = rg * b
po = 1j
ai = po * a

# Phi @ w = -ia for gravity sol.
phia = -so * ((4 * a + b) * (1 - pr) + b) / (8 * (1 - pr) * (2 * a + b))

# Phu @ w = -ia for tectonic sol.
phiat = -(ts * b) / (4 * (2 * a + b))

# d2phia is the 2nd derivative of phi @ w = -ia for tectonic sol.
d2phia = -(ts * b * (4 * a + b) * (b - 12 * a))
d2phia = d2phia / (2 * a * (2 * a + b) * (4 * a + b) ** 3)

ulist = np.arange(umin, umax, uinc)
vlist = np.arange(vmin, vmax, vinc)

[U, V] = np.meshgrid(ulist, vlist)

w = U + 1j * V
r = np.sqrt(np.square(U) + np.square(V + a))

# conformal map
z = w + (a * b) / (w - ai)

X = np.real(z)
Y = np.imag(z)
dz = (np.square(w - ai) - a * b) / np.square(w - ai)

# A(w)
aw1 = po * (4 * a + b) / (8 * (w - ai))
aw2 = a * b * (w - 3 * ai) / (8 * (1 - pr) * np.power(w - ai, 3))
aw = -aw1 - aw2

# phi(w) for gravity solution
phi = -aw * so / dz + a * b * phia / (dz * np.square(w - ai))

# A0(w)
awt = -(a * b * ts) / (2 * np.square(w - ai))

# phi0(w)
phit = -awt / dz + a * b * phiat / (dz * np.square(w - ai))

# sim(sigmaxx + sigmayy) for tectonic and gravity sol.
ssumt = 4 * np.real(phit) + ts
ssum = ssumt + 4 * np.real(phi) + rg * Y / (1 - pr)

# 1st derivative of phi(w)
d1aw1 = po * (4 * a + b) / (8 * np.square(w - ai))
d1aw2 = 2 * a * b / (8 * (1 - pr) * np.power(w - ai, 3))
d1aw3 = 6 * po * a**2 * b / (8 * (1 - pr) * np.power(w - ai, 4))
d1aw = d1aw1 + d1aw2 - d1aw3
d1phi = -so * d1aw / dz - (2 * a * b * (phi + phia)) / (dz * np.power(w - ai, 3))

# 2nd derivative of phi(w) to be used @ w = -ia
d2phi1 = 2 * phi / np.square(w - ai)
d2phi2 = 4 * d1phi / (w - ai)
d2phi3 = so * po * a**2 * b / (2 * (1 - pr) * np.power(w - ai, 5))
d2phi = -(d2phi1 + d2phi2 + d2phi3) / dz

# B(w)
bw1 = po * (4 * a + b) / (8 * (w - ai))
bw2 = (1 - 2 * pr) * a * b * (w - 3 * ai) / (8 * (1 - pr) * np.power(w - ai, 3))
bw = -so * (bw1 + bw2)

psi1 = w * d1phi + bw + phi
d1phi1 = -(a * b * ts * (4 * a + b) * (w - ai))
d1phi2 = 2 * (2 * a + b) * np.square(np.square(w - ai) - a * b)
d1phit = d1phi1 / d1phi2
bwt = -awt

# part of (sigmayy - sigmaxx + 2isigmaxy)
psi1t = w * d1phit + bwt + phit

# test on closeness of w to -ia
# if w is near -ia, the Taylor expansion about -ia is used
psi2 = np.zeros(U.shape, dtype="complex_")
psi2t = np.zeros(U.shape, dtype="complex_")

psi2[r < rad] = 0.5 * a * b * d2phi[r < rad]
psi2t[r < rad] = 0.5 * a * b * d2phia

psi2[r >= rad] = a * b * d1phi[r >= rad] / (w[r >= rad] + ai) - a * b * (
    phi[r >= rad] - phia
) / np.square(w[r >= rad] + ai)
psi2t[r >= rad] = a * b * d1phit[r >= rad] / (w[r >= rad] + ai) - a * b * (
    phit[r >= rad] - phiat
) / np.square(w[r >= rad] + ai)

psi = -(psi1 + psi2) / dz
psit = -(psi1t + psi2t) / dz
zbar = np.conjugate(z)

str1 = 2 * (zbar * d1phi / dz + psi)
str1t = 2 * (zbar * d1phit / dz + psit)

# Sums and differences of sigmas
dift = np.real(str1t) - ts
dif = dift + np.real(str1) + rg * Y * (1 - 2 * pr) / (1 - pr)
sigxy = np.imag(str1) + np.imag(str1t)
sigxy = sigxy / 2
sigx = (ssum - dif) / 2
sigy = (ssum + dif) / 2

# Principal stress directions
alfa = 0.5 * np.arctan2(2 * sigxy, sigx - sigy)
alfa = np.rad2deg(alfa)

# Principal stress magntudes
sig1 = ssum / 2 + np.sqrt(np.square(dif / 2) + np.square(sigxy))
sig2 = ssum / 2 - np.sqrt(np.square(dif / 2) + np.square(sigxy))

# Maximum shear stress
tmax = (sig1 - sig2) / 2

# Out of plane stress
sigz = pr * (sigx + sigy)

# ******************************************************************************
# End of RVT routine
#
# ******************************************************************************

# Plot the resultant stress fields
fig = plt.figure()
ax1 = fig.add_subplot(1, 3, 1)
cf1 = ax1.contourf(X, Y, sigx, 10, cmap="jet")
ax1.set_xlabel("X/b")
ax1.set_ylabel("Y/b")
ax1.set_title("$\sigma_{xx}$")
ax1.set_aspect("equal")

divider = mpl.make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf1, cax=cax)

ax2 = fig.add_subplot(1, 3, 2)
cf2 = ax2.contourf(X, Y, sigy, 10, cmap="jet")
ax2.set_xlabel("X/b")
ax2.set_ylabel("Y/b")
ax2.set_title("$\sigma_{yy}$")
ax2.set_aspect("equal")

divider = mpl.make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf2, cax=cax)

ax3 = fig.add_subplot(1, 3, 3)
cf3 = ax3.contourf(X, Y, sigxy, 10, cmap="jet")
ax3.set_xlabel("X/b")
ax3.set_ylabel("Y/b")
ax3.set_title("$\sigma_{xy}$")
ax3.set_aspect("equal")

divider = mpl.make_axes_locatable(ax3)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf3, cax=cax)

plt.show()
