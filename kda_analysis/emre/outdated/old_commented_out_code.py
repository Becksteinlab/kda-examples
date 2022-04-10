
# ===============================================================================
# ============== 3D Flux over k_AA_a and k_AA_s Code ============================
# ===============================================================================

# amin = 5e2
# amax = 1e3
# smin = 1e-1
# smax = 5e-1
# divisions = 1e2
# a_step = (amax - amin)/divisions
# s_step = (smax - smin)/divisions
# k_AA_a = np.arange(amin, amax+a_step, a_step)
# k_AA_s = np.arange(smin, smax+s_step, s_step)
#
# H_in = 10**(-6.5)  # M
# H_out = 10**(-7.5) # M
# D_in = 25e-9/np.sqrt(10)
# D_out = 25e-9*np.sqrt(10)
#
# def calc_k_AA_Z(x, y, funcs, kAAa=k_AA_a, kAAs=k_AA_s, Din=D_in, Dout=D_out):
#     Z = np.zeros_like(x)
#     for i, ka in enumerate(kAAa):
#         for j, ks in enumerate(kAAs):
#             norm = norm_lf(Din, Dout, ka, ks)
#             if norm == 0:
#                 norm = -1
#             vals = []
#             for func in funcs:
#                 val = func(Din, Dout, ka, ks)
#                 vals.append(val/norm)
#             vals_tot = np.sum(vals)
#             Z[j, i] = vals_tot
#     return Z
#
# X, Y = np.meshgrid(k_AA_a, k_AA_s)
# ZH = calc_k_AA_Z(X, Y, H8_lfs)
# ZD = calc_k_AA_Z(X, Y, D8_lfs)
#
# mag_vals = np.array([np.abs(np.max(ZH)), np.abs(np.max(ZD)), np.abs(np.min(ZH)), np.abs(np.min(ZD))])
#
# v_max = 1.05*np.max(mag_vals)
# v_min = -v_max
#
# def plot_3D_k_AA_flux(X, Y, Z, vmin, vmax, azimuth=-50, elev=40):
#     fig = plt.figure(figsize=(9, 6), tight_layout=True)
#     ax = fig.add_subplot(111, projection='3d')
#     cmap = plt.cm.bwr
#     norm = matplotlib.colors.Normalize(vmin, vmax)
#     zlabel = r'Flux ($s^{-1}$)'
#     ax.set_xlabel(r'$k_{anti}$')
#     ax.set_ylabel(r'$k_{sym}$')
#     ax.set_zlim(Z.min(), Z.max())
#     ax.view_init(elev=elev, azim=azimuth)
#     ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10, linewidth=0.8, color="black")
#     surf = ax.plot_surface(X, Y, Z, cmap=cmap, norm=norm, alpha=0.65)
#     cset = ax.contourf(X, Y, Z, 10, zdir='z', offset=Z.min(), cmap=cmap, norm=norm)
#     cb = fig.colorbar(surf, cmap=cmap, norm=norm, shrink=0.5, aspect=7)
#     cb.set_label(zlabel)

# ===============================================================================
# ============== ACTIVE TRANSPORT CODE ==========================================
# ===============================================================================

# def active_transport(A_in, A_out, A_flux):
#     # A_in, A_out, A_flux must be equal length arrays
#     Ar = A_in/A_out
#     i = 0
#     flux_indices = []
#     for flux, ar in zip(A_flux, Ar):
#         if flux > 0:
#             if ar < 1:
#                 print("Active Transport for Tr, Flux = {}, {:.2e}".format(ar, flux))
#                 flux_indices.append(i)
#         elif flux < 0:
#             if ar > 1:
#                 print("Active Transport for Tr, Flux = {}, {:.2e}".format(ar, flux))
#                 flux_indices.append(i)
#         else:
#             print("No Active Transport for Tr, Flux = {}, {:.2e}".format(ar, flux))
#         i += 1
#     return flux_indices
#
# D_inside = D_in*np.ones(len(ZD))
# D_outside = D_out*np.ones(len(ZD))
# print("Ratio: {:.1f}".format(D_inside[0]/D_outside[0]))
# at_list = []
# for i in range(len(ZD[0])):
#     at = active_transport(D_inside, D_outside, ZD[i, :])
#     if at != []:
#         at_list.append([i, at])

# ===============================================================================
# ============== OLD CODE =======================================================
# ===============================================================================

# plot_cycle_flux(RAA, H8_cycles_valid, H_fluxes, H_flux, colors, species='H', raw=True)
# plot_cycle_flux(RAA, D8_cycles_valid, D_fluxes, D_flux, colors, species='D', raw=True)

# ===========================
# == Old Data: 3D plots =====
# ===========================

# min_power = 1
# max_power = 1
# step1 = 0.01
# step2 = 0.01
# powers = np.arange(-min_power, max_power + step1, step1)
# powers2 = 0.5*np.arange(-min_power, max_power + step2, step2)
# ratio = 10.0**powers
# ratio2 = 10.0**powers2
# k_AA_factor = 1e1
# D_factor = 25e-9
#
# k_AA_a = []
# k_AA_s = []
# for r in ratio:
#     k_AA_a.append(k_AA_factor*r)
#     k_AA_s.append(k_AA_factor/r)
#
# D_in = []
# D_out = []
# for r2 in ratio2:
#     D_in.append(D_factor/r2)
#     D_out.append(D_factor*r2)
#
# multiplier = 1
# multiplier2 = 1e-2
# k_AA_a = 1e-1*np.array(k_AA_a)
# k_AA_s = 1e1*np.array(k_AA_s)
# D_in = np.array(D_in)
# D_out = np.array(D_out)
#
# RAA = k_AA_a/k_AA_s
# Tr = D_in/D_out
#
# H_in = 10**(-6.5)  # M
# H_out = 10**(-7.5) # M
# Hr = 10*np.ones(len(Tr))

# ==================
# == Old 3D func ===
# ==================

# def calc_RAA_Tr_Z(x, y, funcs, raa=RAA, tr=Tr, Din=D_in, Dout=D_out, ka=k_AA_a, ks=k_AA_s):
#     Z = np.zeros_like(x)
#     for i in range(len(raa)):
#         for j in range(len(tr)):
#             norm = norm_lf(Din[j], Dout[j], ka[i], ks[i])
#             vals = []
#             for func in funcs:
#                 val = func(Din[j], Dout[j], ka[i], ks[i])
#                 vals.append(val/norm)
#             vals_tot = np.sum(vals)
#             Z[j, i] = vals_tot
#     return Z

# X, Y = np.meshgrid(RAA, Tr)
# ZH = calc_Z(X, Y, H8_lfs)
# ZD = calc_Z(X, Y, D8_lfs)

# =============================================
# == Old Net Cycle Flux By Cycle Data =========
# =============================================

# H_flux = []
# D_flux = []
# H_fluxes = []
# D_fluxes = []

# H_tfluxes = []
# D_tfluxes = []

# RAA = []
# Tr = []
# Hr = []

# H_in = 10**(-6.5)  # M
# H_out = 10**(-7.5) # M

# powers = np.arange(-1, 1, 5e-5)
# powers2 = 0.5*np.arange(-1, 1, 5e-5)
# ratio = 10.0**powers
# ratio2 = 10.0**(0.5*np.sin(powers2*3.1))
# factor = 1e1
# D_factor = 25e-9

# for r, r2 in zip(ratio, ratio2):
#     ka = factor*r
#     ks = factor/r
# #     print(ka, ks)
#     D_in = D_factor/r2
#     D_out = D_factor*r2
#     RAA.append(ka/ks)
#     Tr.append(D_in/D_out)
#     Hr.append(H_in/H_out)
#     norm = norm_lf(D_in, D_out, ka, ks)
#     H_vals = [Hfunc(D_in, D_out, ka, ks)/norm for Hfunc in H8_lfs]
#     D_vals = [Dfunc(D_in, D_out, ka, ks)/norm for Dfunc in D8_lfs]
#     H_flux.append(np.sum(H_vals))
#     D_flux.append(np.sum(D_vals))
#     H_fluxes.append(H_vals)
#     D_fluxes.append(D_vals)
#     H_tvals = [Htfunc(D_in, D_out, ka, ks)/norm for Htfunc in H8_tlfs]
#     D_tvals = [Dtfunc(D_in, D_out, ka, ks)/norm for Dtfunc in D8_tlfs]
#     H_tfluxes.append(H_tvals)
#     D_tfluxes.append(D_tvals)

# H_flux = np.array(H_flux)
# D_flux = np.array(D_flux)
# H_fluxes = np.array(H_fluxes).T
# D_fluxes = np.array(D_fluxes).T

# H_tfluxes = np.array(H_tfluxes).T
# D_tfluxes = np.array(D_tfluxes).T

# RAA = np.array(RAA)

# =====================
# == Misc plots =======
# =====================

# fig = plt.figure(figsize = (10, 7), tight_layout=True)
# ax = fig.add_subplot(111)
# ax.semilogx(RAA, H_flux, '-', lw=2, color="#A02020", label=r"$J_{H^{+}}$")
# # ax.semilogx(RAA, D_flux, '-', lw=2, color="#6BE35D", label=r"$J_{D}$")
# ax.axhline(y=0, xmin=0, xmax=1, ls='-', color='black')
# ax.axvline(x=1e0, ymin=0, ymax=1, ls='-', color='black')
# # ax.set_xlim(2244, 2248)
# ax.set_title("Proton and Drug Flux for EmrE")
# ax.set_ylabel(r"Flux ($s^{-1}$)")
# ax.set_xlabel(r"$R_{AA}$")
# ax.legend(loc='best')
# ax.grid(True)

# fig = plt.figure(figsize = (10, 7), tight_layout=True)
# ax = fig.add_subplot(111)
# # ax.semilogx(RAA, H_flux, '-', lw=2, color="#A02020", label=r"$J_{H^{+}}$")
# ax.semilogx(RAA, D_flux, '-', lw=2, color="#6BE35D", label=r"$J_{D}$")
# ax.axhline(y=0, xmin=0, xmax=1, ls='-', color='black')
# ax.axvline(x=1e0, ymin=0, ymax=1, ls='-', color='black')
# # ax.set_xlim(2244, 2248)
# ax.set_title("Proton and Drug Flux for EmrE")
# ax.set_ylabel(r"Flux ($s^{-1}$)")
# ax.set_xlabel(r"$R_{AA}$")
# ax.legend(loc='best')
# ax.grid(True)

# H_mag = np.abs(H_flux)
# D_mag = np.abs(D_flux)
# norm = H_mag + D_mag
# data_norm = (H_mag/D_mag)
# fig2 = plt.figure(figsize = (10, 7), tight_layout=True)
# ax = fig2.add_subplot(111)
# ax.loglog(RAA, data_norm, ls='-', lw=2, color="red", label=r"$J_{H^{+}}$ / $J_{D}$")
# # ax.axvline(x=1e0, ymin=0, ymax=1, ls='--', color='black')
# # ax.axhline(y=1e0, xmin=0, xmax=1, ls='--', color='black')
# # ax.set_xlim(0.5, 0.9)
# ax.set_title("Proton-Drug Flux Ratio")
# ax.set_ylabel("Flux Ratio")
# ax.set_xlabel(r"$R_{AA}$")
# ax.legend(loc='best')
# ax.grid(True)


# =====================
# == Rate mappings ====
# =====================
# k31 = k1
# k13 = k2
# k57 = k3
# k75 = k4
# k42 = k5
# k24 = k6
# k68 = k7
# k86 = k8
# k34 = k9
# k43 = k10
# k56 = k11
# k65 = k12
# k12 = k13
# k21 = k14
# k78 = k15
# k87 = k16
# k71 = k17
# k17 = k18
# k53 = k19
# k35 = k20
# k64 = k21
# k46 = k22
# k82 = k23
# k28 = k24
# k109 = k25
# k910 = k26
# k19 = k27
# k91 = k28
# k710 = k29
# k107 = k30
