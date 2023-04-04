# JESUS

import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt

# Functions

def disp_deltaT(i, A, j, alpha, P, x1, x2, thick, depth, dip):
    
    points = 1000
    i = float(i)
    A = float(A)
    j = float(j)
    alpha = float(alpha)
    P = float(P)
    x1 = float(x1)
    x2 = float(x2)
    thick = float(thick)
    depth = float(depth)
    dip = float(dip)
    
    x = np.linspace(x1, x2, points)
    A = np.deg2rad(A)
    i = np.deg2rad(i)
    j = np.deg2rad(j)
    alpha = np.deg2rad(alpha)
    dip = np.deg2rad(dip)
    
    def calc_b_beta(i, A, j, alpha):
        b = 1 - np.cos(i)*np.sin(A)
        beta = 1 - np.cos(j)*np.sin(alpha)
        return b, beta

    def calc_amplitude(P, b, beta, dip):
        return 2*P*b*beta*np.sin(dip)

    def calc_R1_R2_theta(x_position, thick, depth):
        R1 = np.sqrt(depth**2 + (x_position - thick/2)**2)
        R2 = np.sqrt(depth**2 + (x_position + thick/2)**2)
        theta = np.zeros(len(x_position))
        for idx in range(len(x_position)):
            if x_position[idx] >= -thick/2 and x_position[idx] <= thick/2:
                theta[idx] = np.pi - abs(np.arcsin(depth/(R1[idx]))) - abs(np.arcsin(depth/(R2[idx])))
            else:
                theta[idx] = abs(np.arcsin(depth/(R1[idx])) - np.arcsin(depth/(R2[idx])))
        return R1, R2, theta

    def calc_I_J(i, A, j, alpha):
        I = np.arctan2(np.tan(i), np.cos(A))
        J = np.arctan2(np.tan(j), np.cos(alpha))
        return I, J

    b, beta = calc_b_beta(i, A, j, alpha)
    amp = calc_amplitude(P, b, beta, dip)
    R1, R2, theta = calc_R1_R2_theta(x, thick, depth)
    I, J = calc_I_J(i, A, j, alpha)
    H = I + J - dip
    deltaT = amp * (theta * np.sin(H) + np.log(R2 / R1) * np.cos(H))
    
    plt.plot(x, deltaT, ls="--", lw=2)
    plt.xlabel("Profile Distance X (m)")
    plt.ylabel("Anomalous Feild Strength (nT)")
    fig = plt.gcf()
    fig.set_size_inches(8, 10.5)
    plt.tight_layout()
    plt.savefig("tempsave.png")
    
    disp_image = tk.PhotoImage(file="tempsave.png")
    disp_label = tk.Label(disp_frame, image=disp_image)
    disp_label.place(relx=0, rely=0, relw=1, relh=1)

def plot_data_points():
    
    filename = filedialog.askopenfilename(initialdir='/',title="Select .csv file", \
                                          filetypes=(("comma seperated values", "*.csv"), \
                                          ("all files","*.*")))
    data = np.loadtxt(filename, delimiter = ",")
    mag_readings = data[:,0]
    x_position = data[:,1]
    
    plt.plot(x_position, mag_readings, ls=".", lw=2, color="#000000")
    plt.xlabel("Profile Distance X (m)")
    plt.ylabel("Anomalous Feild Strength (nT)")
    fig = plt.gcf()
    fig.set_size_inches(8, 10.5)
    plt.tight_layout()
    plt.savefig("tempsave.png")
    
    disp_image = tk.PhotoImage(file="tempsave.png")
    disp_label = tk.Label(disp_frame, image=disp_image)
    disp_label.place(relx=0, rely=0, relw=1, relh=1)
    
def clear_curve():
    ax = plt.gca()
    
    ax.lines.pop(-1)
    plt.savefig("tempsave.png")
    
    disp_image = tk.PhotoImage(file="tempsave.png")
    disp_label = tk.Label(disp_frame, image=disp_image)
    disp_label.place(relx=0, rely=0, relw=1, relh=1)
    
def export_figure():
    filename = "ExportDikeMagAnoms.png"
    fig = plt.gcf()
    fig.set_tight_layout(False)
    plt.title("Magnetic Anomaly due to Flat Topped Dike Model")
    plt.tight_layout()
    plt.savefig(filename)
    
    
root = tk.Tk()

canvas = tk.Canvas(root, height=900, width=1000)
canvas.pack()

# Frames
button_frame = tk.Frame(root, bg="#000000")
button_frame.place(relx=0.025, rely=0.9, relw=0.95, relh=0.075)

var_frame = tk.Frame(root, bg="#ffffff", bd=5)
var_frame.place(relx=0.025, rely=0.125, relw=0.35, relh=0.775)

disp_frame = tk.Frame(root, bg="#ffffff", bd=5)
disp_frame.place(relx=0.375, rely=0.025, relw=0.6, relh=0.875)

title_frame = tk.Frame(root, bg="#e0e0eb")
title_frame.place(relx=0.025, rely=0.025, relw=0.35, relh=0.1)

# Title
title_label1 = tk.Label(title_frame, text="Magnetic Anomaloy Due To", bg="#e0e0eb")
title_label1.place(relx=0, rely=0, relw=1, relh=0.33)
title_label2 = tk.Label(title_frame, text="Flat Topped Dike Model", bg="#e0e0eb")
title_label2.place(relx=0, rely=0.33, relw=1, relh=0.33)
title_label3 = tk.Label(title_frame, text="Induced + Remnant Magnetization", bg="#e0e0eb")
title_label3.place(relx=0, rely=0.66, relw=1, relh=0.34)

# Buttons
Add_Curve = tk.Button(button_frame, bg="#7575a3", text="Add Curve", \
                      relief="raised", command=lambda: disp_deltaT(entry_i.get(), \
                      entry_A.get(), entry_j.get(), entry_alpha.get(), \
                      entry_P.get(), entry_x1.get(), entry_x2.get(), \
                      entry_thick.get(), entry_depth.get(), entry_dip.get()))
Add_Curve.place(relx=0, rely=0, relw=0.25, relh=1)

Import_Curve = tk.Button(button_frame, bg="#7575a3", text="Import Curve", \
                         relief="raised", command=lambda: plot_data_points())
Import_Curve.place(relx=0.25, rely=0, relw=0.25, relh=1)

Clear_Curve = tk.Button(button_frame, bg="#7575a3", text="Clear Last Curve", \
                        relief="raised", command=lambda: clear_curve())
Clear_Curve.place(relx=0.5, rely=0, relw=0.25, relh=1)

Export_Figure = tk.Button(button_frame, bg="#7575a3", text="Export", \
                          relief="raised", command=lambda: export_figure())
Export_Figure.place(relx=0.75, rely=0, relw=0.25, relh=1)

# User Inputs
label_A = tk.Label(var_frame, text="Angle between traverse and magnetic north", \
                  bg="#ffffff")
label_A.place(relx=0, rely=0, relw=1, relh=0.025)
unit_A = tk.Label(var_frame, text="degrees", bg="#ffffff")
unit_A.place(relx=0.8, rely=0.037, relw=0.2, relh=0.025)
entry_A = tk.Entry(var_frame, bg="#e6e6e6")
entry_A.place(relx=0, rely=0.025, relw=0.8, relh=0.05)

label_i = tk.Label(var_frame, text="Inclination of Earths magnetic feild", \
                  bg="#ffffff")
label_i.place(relx=0, rely=0.1, relw=1, relh=0.025)
unit_i = tk.Label(var_frame, text="degrees", bg="#ffffff")
unit_i.place(relx=0.8, rely=0.137, relw=0.2, relh=0.025)
entry_i = tk.Entry(var_frame, bg="#e6e6e6")
entry_i.place(relx=0, rely=0.125, relw=0.8, relh=0.05)

label_P = tk.Label(var_frame, text="Magnetic polarization of dike", \
                  bg="#ffffff")
label_P.place(relx=0, rely=0.2, relw=1, relh=0.025)
unit_P = tk.Label(var_frame, text="A/m*10^-9", bg="#ffffff")
unit_P.place(relx=0.75, rely=0.237, relw=0.25, relh=0.025)
entry_P = tk.Entry(var_frame, bg="#e6e6e6")
entry_P.place(relx=0, rely=0.225, relw=0.75, relh=0.05)

label_x = tk.Label(var_frame, text="Transverse endpoints (m)", \
                  bg="#ffffff")
label_x.place(relx=0, rely=0.3, relw=1, relh=0.025)
unit_thick = tk.Label(var_frame, text=" to ", bg="#ffffff")
unit_thick.place(relx=0.4, rely=0.337, relw=0.2, relh=0.025)
entry_x1 = tk.Entry(var_frame, bg="#e6e6e6")
entry_x1.place(relx=0, rely=0.325, relw=0.4, relh=0.05)
entry_x2 = tk.Entry(var_frame, bg="#e6e6e6")
entry_x2.place(relx=0.6, rely=0.325, relw=0.4, relh=0.05)

label_thick = tk.Label(var_frame, text="Thickness of dike", \
                  bg="#ffffff")
label_thick.place(relx=0, rely=0.4, relw=1, relh=0.025)
unit_thick = tk.Label(var_frame, text="m", bg="#ffffff")
unit_thick.place(relx=0.8, rely=0.437, relw=0.2, relh=0.025)
entry_thick = tk.Entry(var_frame, bg="#e6e6e6")
entry_thick.place(relx=0, rely=0.425, relw=0.8, relh=0.05)

label_depth = tk.Label(var_frame, text="Depth of dike", \
                  bg="#ffffff")
label_depth.place(relx=0, rely=0.5, relw=1, relh=0.025)
unit_depth = tk.Label(var_frame, text="m", bg="#ffffff")
unit_depth.place(relx=0.8, rely=0.537, relw=0.2, relh=0.025)
entry_depth = tk.Entry(var_frame, bg="#e6e6e6")
entry_depth.place(relx=0, rely=0.525, relw=0.8, relh=0.05)

label_dip = tk.Label(var_frame, text="Dip of dike", \
                  bg="#ffffff")
label_dip.place(relx=0, rely=0.6, relw=1, relh=0.025)
unit_dip = tk.Label(var_frame, text="degrees", bg="#ffffff")
unit_dip.place(relx=0.8, rely=0.637, relw=0.2, relh=0.025)
entry_dip = tk.Entry(var_frame, bg="#e6e6e6")
entry_dip.place(relx=0, rely=0.625, relw=0.8, relh=0.05)

label_alpha = tk.Label(var_frame, text="Angle between traverse and projection of remnant onto xy", \
                  bg="#ffffff")
label_alpha.place(relx=0, rely=0.7, relw=1, relh=0.025)
unit_alpha = tk.Label(var_frame, text="degrees", bg="#ffffff")
unit_alpha.place(relx=0.8, rely=0.737, relw=0.2, relh=0.025)
entry_alpha = tk.Entry(var_frame, bg="#e6e6e6")
entry_alpha.place(relx=0, rely=0.725, relw=0.8, relh=0.05)

label_j = tk.Label(var_frame, text="Inclination of remnant magnetization from horizontal", \
                  bg="#ffffff")
label_j.place(relx=0, rely=0.8, relw=1, relh=0.025)
unit_j = tk.Label(var_frame, text="degrees", bg="#ffffff")
unit_j.place(relx=0.8, rely=0.837, relw=0.2, relh=0.025)
entry_j = tk.Entry(var_frame, bg="#e6e6e6")
entry_j.place(relx=0, rely=0.825, relw=0.8, relh=0.05)

root.mainloop()
# # anlge between traverse (xaxis) and magnetic north
# A = 0

# # inclination of Earth's magnetic feild (positive in north, 90 at mag pole)
# i = 90

# # magnetic polarization of the dike, (== magnetic suceptability * TMI)
# P = 100

# # profile to transverse the dike over
# x = np.linspace(-200, 220, 1000)

# # thickness of dike in x direction
# thick = 20

# # depth of the top of the flat-topped dike
# depth = 50

# # dip of the dike, measured from the positive xaxis
# dip = 45

# # inclination of the remnant magnetisation from the horizontal [-180, 180]
# j = 30

# # angle between positive xaxis and projection of magnitisation on the hoizontal plane 
# # [-180, 180], alpha == A for purely induced magnitisation
# alpha = 0

# # anomalous magnetic feild
# T = calc_deltaT(i, A, j, alpha, P, x, thick, depth, dip)

# plt.plot(x, T, 'r-', lw=2, label = "T")
# xl1 = [0, thick]
# xl2 = [0,0]
# plt.plot(xl1, xl2, 'b-', lw=2, label = "Dike: dip = {}, thickness = {}, depth = {}".format(dip, thick, depth))
# plt.legend()
# plt.title("Flat topped dike model")
# plt.xlabel("Transverse Profile (m)")
# plt.ylabel("Anomalous Feild (nT)")