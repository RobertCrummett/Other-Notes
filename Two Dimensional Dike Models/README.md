# README.md for Flat Topped Dike Magnetic Anomaly Calculators
These codes utilize the relationships provided by Colin Reeves (2005), which was originally from Reford and Sumner, (1964).\
The relationship is calculated in Python. The libraries pyinstaller and tkinter have been used to convert the code into a GUI executable.

## EXPLANATION OF FIGURE
The x axis is a traverse taken over a flat topped dike, with the top centered at x = 0.
This model assumes the dike is infinitely deep and long in the y-direction. 
See the papers listed for further figures and explanations: Reeves (2005), Reeves (1984), Reford and Sumner (1964)

## FUNCTIONALITY 
There are four buttons:
ADD CURVE: This button takes the user input values and calculated the theoretical anomalous magnetic feild resulting from the dike.
IMPORT CURVE: For real data, the program can add this data to the figure for comparison. The data must be organized in a two columned
    .csv file, with the filtered anomalous mag readings in the first column, and the xpoisition along the traverse on the second.
CLEAR LAST CURVE: This button clears the last curve added.
EXPORT: This button saves the current figure as "ExportDikeMagAnoms.png"

## USER INPUTS
A = anlge between traverse (x axis) and magnetic north [degrees]
i = inclination of Earth's magnetic feild (positive in north, 90 at mag pole) [degrees]
P = magnetic polarization of the dike (magnetic suceptability * TMI) [units of TMI]
x = profile to transverse the dike over [m]
thick = thickness of dike along x axis [m]
depth = depth from flat surface to the top of the flat-topped dike [m]
dip = dip of the dike, measured from the positive x axis (0 is horizontal) [degrees]
j = inclination of the remnant magnetisation from the horizontal (-180, 180) [degrees]
alpha = angle between positive x axis and projection of magnitization onto the hoizontal plane (-180, 180) [degrees]