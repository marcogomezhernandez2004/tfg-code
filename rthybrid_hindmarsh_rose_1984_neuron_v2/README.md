### RTHybrid Hindmarsh-Rose (1984) neuron model

**Requirements:** None  
**Limitations:** None  

![RTHybrid Hindmarsh-Rose (1984) neuron model GUI](rthybrid_hindmarsh_rose_1984_neuron.png)


#### RTHybrid Hindmarsh-Rose (1984) neuron model
RTHybrid module for RTXI that implements the [Hindmarsh-Rose (1984) neuron model](https://doi.org/10.1098%2Frspb.1984.0024).

#### Input
1. input(0) - Isyn (nA) : Synaptic input current (in nA)
2. input(1) - Burst duration (s) : Burst duration, from the first spike of one burst to the first spike of the next (in seconds).

#### Output
1. output(0) - Vm (v) : Membrane potential (in V)
2. output(1) - Vm (mV) : Membrane potential (in mV)

#### Parameters
1. Burst duration (s) - Burst duration, from the first spike of one burst to the first spike of the next (in seconds). If you want to use the value that comes as input(1), type -1 in this field.
2. I - 
3. a - 
4. b - 
5. c - 
6. d - 
7. r - 
8. s - 
9. xr - 
10. Vh - 
11. x0 - 
12. y0 - 
13. z0 - 

#### States
1. v - 
2. s_points - 
3. dt - 
4. syn - 
