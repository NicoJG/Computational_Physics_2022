# 1 Coupled harmonic oscillators

$$
h(t) = a \cos(2 \pi f t + \phi)
$$

## 1) 
- f=2, phi=0 => cosine with a lot of oscillations  
<img src="images/signal_f2_phi0.png" width="500" height="400">  

- f=1, phi=0 => cosine with half of the   oscillations  
<img src="images/signal_f1_phi0.png" width="500" height="400">  
- f=1, phi=pi/2 => -sine with half of the oscillations  
<img src="images/signal_f1_phi_pihalf.png" width="500" height="400">  

- overall with N=250 the oscillations are not very smooth  

## 2)
- using f=2, phi=0  
- we get a power sectrum with peaks only at f=2 and f=8
(same as figure 2)  
<img src="images/powerspectrum.png" width="500" height="400">  
- f=2 is what we expect because this was used to generate the data  
- f=8 is an artifact of the fft (how precise the oscillations are modeled => dt)  
- a cosine with f=8 would generate the same datapoints, but then there are 2 oscillations per dt  
- with dt = 0.01 instead of 0.1 we get a second peak at much higher frequencies:  
<img src="images/powerspectrum_dt0.01.png" width="500" height="400">  
<img src="images/signal_dt0.01.png" width="250" height="200">  

- real signals always have a double (mirrored) FFT (negative frequencies)  

## 3)  
- for f=2, phi=0, N=250 and shifted FFT frequencies:  
<img src="images/powerspectrum_shift.png" width="500" height="400">  
- this is more convenient, because now we have negative frequencies and not very large frequencies.  
- cosine is symmetric around 0 => cos(x)=cos(-x)  
- FFT is periodic in the frequencies  
- Frequency range: $-f_c \leq 0 \leq f_c$ with $f_c = 1/(2 \Delta t) = 5$  
- Frequency precision: $\Delta f = 1/(N \Delta t) = 1/t_\text{max} = 1/25 = 0.04$  
- N=253:  
<img src="images/powerspectrum_N253.png" width="500" height="400">  
- N=255:  
<img src="images/powerspectrum_N255.png" width="500" height="400">  
- with a different number of samples, the frequencies may "leak out" of the single point of the correct frequency => noise in the fft  
- Why?  

## 4)  
- now we generate $h(t) = a_1 \cos(2 \pi f_1 t + \phi_1) + a_2 \cos(2 \pi f_2 t + \phi_2)$ with $a_1 = a_2 = 1, \phi_1 = \phi_2 = 0$ but $f_1 = 2$ and $f_2 = 6$ (and again $N=250, \Delta t=0.1$)  
- Nyquist frequency is still $f_c = 1/(2 \Delta t) = 5$  
- power spectrum:  
<img src="images/powerspectrum_ex4.png" width="500" height="400">  
- signal:  
<img src="images/signal_ex4.png" width="500" height="400">  
- now we see an additional peak at $f=\pm 4$, even though we have $f_2 = 6$  
- our dt is to big to see $f=6$ => decrease dt to 0.05 so that $f_c = 10$:  
<img src="images/powerspectrum_ex4_dt0.05.png" width="500" height="400">  
<img src="images/signal_ex4_dt0.05.png" width="500" height="400">  


## 5)  
- $H = \sum_{i=1}^{3} \frac{m v_i^2}{2} + \sum_{i=0}^{3} \frac{\kappa}{2}(q_{i+1}-q_i)^2$    
- $m \frac{d^2}{dt^2} q_i(t) = \kappa [q_{i+1}(t) - 2q_i(t) + q_{i-1}(t)]$    
- Velocity Verlet:  
    - $v_i(t+\Delta t /2) = v_i(t) + \frac{1}{2}a_i(t) \Delta t$  
    - $q_i(t+\Delta t) = q_i(t) + v_i(t+\Delta t/2) \Delta t$  
    - calculate $a_i(t+\Delta t)$  
    - $v_i(t+\Delta t) = v_i(t + \Delta t /2) + \frac{1}{2} a_i(t+\Delta t) \Delta t$  
- with initial conditions $q_1(0) = 0.01 Å, q_2(0) = 0.005 Å, q_3(0) = -0.005 Å$ and $v_i(0)=0$:  
    - $N=1000$, $t_{max} = 0.25$:  
    <img src="E1_5_N1000_tmax0.25.png" width="500" height="750">  
    - $N=100$, $t_{max} = 0.25$:  
    <img src="E1_5_N100_tmax0.25.png" width="500" height="750">  
    - $N=1000$, $t_{max} = 1$:  
    <img src="E1_5_N1000_tmax1.png" width="500" height="750">  
- Energy is conserved but with less larger timesteps it oscillates  
- The powerspectrum should have a maximum of 3 peaks (we have 2 active modes)  


## 6) 3 coupled oscillators without walls  
<img src="E1_6_7_derivation.png" width="500" height="600">   
<img src="E1_6.png" width="500" height="750">   

- resonance frequencies: 0, (74-75.5) THz, (38.5-39.5) THz   
- experimentally found resonances:  
    - wavenumbers $\nu$: 1388, 2349, 667 1/cm  
    - frequency $f = c \nu$: 41.6112, 70.4212, 19.9962 THz   
https://en.wikipedia.org/wiki/Carbon_dioxide#Structure,_bonding_and_molecular_vibrations 

## 7)  
- analytical frequencies:  
    - f_1 = 0 THz  
    - f_2 = 74.7506 THz  
    - f_3 = 39.0508 THz  


<img src="E1_7_derivation_freq.png" width="500" height="400"> 