# %%
import numpy as np
import matplotlib.pyplot as plt

def plot_modal_energy(input_file, output_file, xrange_zoom):
    data = np.genfromtxt(input_file, delimiter=",", comments="#")
    print("data shape: ",data.shape)

    t = data[:,0]
    E = data[:,1:]
    E_labels = [f"E_{i+1}" for i in range(E.shape[1])]
    E_colors = [f"C{i}" for i in range(E.shape[1])]
    
    E_total = np.sum(E, axis=1)

    fig = plt.figure(figsize=(8,6))

    plt.sca(plt.subplot2grid((2,2),(1,0), colspan=2))
    plt.axvline(xrange_zoom[0], linestyle="--", color="k", alpha=0.5, label="zoom area")
    plt.axvline(xrange_zoom[1], linestyle="--", color="k", alpha=0.5, label="zoom area")
    plt.plot(t, E_total, "k", label="total energy")
    plt.plot(t,E[:,:5], label=E_labels[:5])
    plt.xlabel("time / arb. units")
    plt.ylabel("energy per mode / arb. units")
    plt.legend(loc="center right")

    plt.sca(plt.subplot2grid((2,2),(0,0)))
    plt.plot(t,E[:,0], label=E_labels[0])
    plt.xlim(xrange_zoom)
    plt.xlabel("time / arb. units")
    plt.ylabel("energy per mode / arb. units")
    plt.legend()

    plt.sca(plt.subplot2grid((2,2),(0,1)))
    plt.plot(t,E[:,1:5], label=E_labels[1:5])
    plt.xlim(xrange_zoom)
    plt.xlabel("time / arb. units")
    plt.ylabel("energy per mode / arb. units")
    plt.legend()

    plt.tight_layout()
    plt.savefig(output_file)
    
print("Plot E2_2 ...")
plot_modal_energy("data/E2_2.csv", "images/E2_2.png", (0,1000))

print("Plot E2_3_alpha0.01 ...")
plot_modal_energy("data/E2_3_alpha0.01.csv", "images/E2_3_alpha0.01.png", (0,1000))
print("Plot E2_3_alpha0.1 ...")
plot_modal_energy("data/E2_3_alpha0.1.csv", "images/E2_3_alpha0.1.png", (0,1000))

print("Plot E2_4_alpha0.01 ...")
plot_modal_energy("data/E2_4_alpha0.01.csv", "images/E2_4_alpha0.01.png", (0,10000))
print("Plot E2_4_alpha0.1 ...")
plot_modal_energy("data/E2_4_alpha0.1.csv", "images/E2_4_alpha0.1.png", (0,10000))
print("Done.")
# %%
from scipy.integrate import cumulative_trapezoid


def plot_time_average(input_file, output_file):
    data = np.genfromtxt(input_file, delimiter=",", comments="#")
    print("data shape: ",data.shape)

    t = data[:,0]
    E = data[:,1:]
    E_timeaverage_labels = [f"<E_{i+1}>_time" for i in range(E.shape[1])]
    
    dt = np.diff(t)[0]
    assert len(np.unique(np.diff(t)))==1
    
    E_timeaverage = np.cumsum(E, axis=0)*(dt/t.reshape(-1,1))
    
    plt.figure(figsize=(8,6))
    plt.plot(t, E_timeaverage, label=E_timeaverage_labels)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("time / arb. units")
    plt.ylabel("time average of the modal energy / arb. units")
    plt.savefig(output_file)
    
print("Plot E2_4_average_alpha0.01 ...")
plot_time_average("data/E2_4_alpha0.01.csv", "images/E2_4_average_alpha0.01.png")
print("Plot E2_4_average_alpha0.1 ...")
plot_time_average("data/E2_4_alpha0.1.csv", "images/E2_4_average_alpha0.1.png")
    