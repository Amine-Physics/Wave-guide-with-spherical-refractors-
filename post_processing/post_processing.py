import pandas as pd
import matplotlib.pyplot as plt  
import numpy as np

## plotting the intensity as a function of position 
df = pd.read_csv("../builds/graph_intensity.csv")
print(df.columns.tolist())
print(df.head())
intensities = df[" Intensities"].values 
x  = df["x "].values

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(x, intensities, "k", linewidth=1.5)
ax.set_xlabel("x")
ax.set_ylabel("Intensity")
ax.set_title("Transverse Intensity Across the Waveguide")
fig.savefig("Intensity.png", dpi=200, bbox_inches="tight")
plt.close(fig)

## Plotting the modal amplitude as a function of position 


df_a = pd.read_csv("../builds/graph_amplitude.csv"); 
print(df_a.columns.tolist());
abs_a  = df_a[" ABS"].values 
m = df_a["mode "].values[0] 
x = df_a[" x "].values 
Re = df_a[" Re "].values 
Im = df_a[" Im "].values 
a = Re + 1j * Im 
phase = np.angle(a); 
phase_unwrapped = np.unwrap(phase);
fig_2 , ax_2 = plt.subplots( 2, 2 , figsize = (10 , 8)) 

## ploting real and imaginary parts
ax_2[0,0].plot(x , Re , "b", label = "Re", linewidth = 1.5) 
ax_2[0,0].plot(x , Im , "r", label = "Im", linewidth = 1.5) 
ax_2[0,0].set_xlabel("x") 
ax_2[0,0].legend() 

## ploting amplitude 

ax_2[0,1].plot(x , abs_a , "k" , linewidth = 1.5) 
ax_2[0,1] .set_xlabel("x"); ax_2[0,1].set_ylabel("|a|") 

## ploting the phase
ax_2[1,0].plot(x , phase_unwrapped , "k" , linewidth = 1.5)
ax_2[1,0].set_xlabel("x"); ax_2[1,0].set_ylabel("phase (rad)") 

## plotting trajectory in complex plane 
ax_2[1,1].plot(Re , Im , "k" , linewidth = 1.5) 
ax_2[1,1].set_xlabel("Re") ; ax_2[1,1].set_ylabel("Im") 
ax_2[1,1].set_aspect("equal") 

fig_2.tight_layout() 
fig_2.savefig("Amplitude.png", dpi = 200 , bbox_inches = "tight")
plt.close(fig_2)



## ploting power spectrum 

df_ps = pd.read_csv("../builds/power_spectrum.csv") 
print(df_ps.columns.tolist())

modes = df_ps["mode "].values
quantities = [" R ", " T ", " port_ref_in ", " port_ref_out"]
labels = ["R", "T", "port reflection (in)", "port reflection (out)"]

n = len(quantities)
width = 0.15                       # thin bars, so they look like sticks
offsets = (np.arange(n) - (n - 1) / 2) * width * 1.3
xpos = np.arange(len(modes))       # one slot per mode

fig, ax = plt.subplots(figsize=(8, 5))
for k, (q, lab) in enumerate(zip(quantities, labels)):
    ax.bar(xpos + offsets[k], df_ps[q].values, width=width, label=lab)

ax.set_xticks(xpos)
ax.set_xticklabels([str(int(m)) for m in modes])
ax.set_xlabel("Mode m")
ax.set_ylabel("Power fraction")
ax.legend()
fig.tight_layout()
fig.savefig("power_spectrum.png", dpi=200)
plt.close(fig)

## compute balance  






























