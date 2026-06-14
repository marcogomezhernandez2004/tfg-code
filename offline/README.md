# Offline chemical synapse model parameterization for biohybrid neural circuits

The goal of this program is to automate the offline parameterization of the **graded chemical synapse model by Golowasch et al.** [1] for biohybrid circuits, where biological neurons and computational models interact in real time (although in this program, the parameterization is performed offline). This model exhibits ***parameter sloppiness*** (multiple parameter combinations can yield similar outputs), making manual tuning difficult. It separates the synaptic current to enable modular and frequency-selective coupling:
- **Fast component**: presynaptic spikes (fast wave/high frequency).
- **Slow component**: presynaptic slow wave (low frequency).

## Key Features

- **Bayesian Optimization (BO)**: models the objective function using a **Gaussian Process (GP)**, configured with:
  - **Kernel**: **Squared Exponential Automatic Relevance Determination (SE-ARD)**; its hyperparameters are reoptimized using the **Resilient Backpropagation (Rprop)** algorithm.
  - **Initialization**: **Latin Hypercube Sampling (LHS)**.
  - **Acquisition function**: **Expected Improvement (EI)**; maximized via the **Subplex** algorithm.
- **Evaluation function** based on frequency decomposition of the presynaptic potential using a **low-pass filter (Butterworth)** to separate fast and slow components to use them as references. It independently assesses the **shape** (**Pearson correlation**) and **range** of each current component.
- **Parameter space handling**:
  - **Parameter sloppiness mitigation**:
    - **Reparameterizations**: `R = k2/k1`.
    - **Logarithmic scales**: conductances, `k1` and `R`.
    - **SE-ARD kernel**.
  - **Dynamic search bounds** adapted to each experiment's signal characteristics.

## Implementation Details

The codebase is written in **C++20**.

This implementation uses a **Runge-Kutta integrator** for solving the differential equations of the computational neuronal and synaptic models.

It operates on a presynaptic voltage recording (CSV) in a **unidirectional** setup (because the presynaptic register can not react dynamicly to currents) with a **Hindmarsh-Rose** [2] neuronal model as the postsynaptic neuron. The CSV can be changed for any other and the neuronal model can also be easily substituted by any other model of Neun (thanks to the use of templates). The provided CSV contains data from two Pyloric Dilator (PD) neurons of the stomatogastric ganglion of a crab, which exhibit the characteristic bursting behavior of Central Pattern Generator (CPG) neurons.

The user can **choose which current components** are used (calculated) and optimized to maximize efficiency.

The **presynaptic voltage is scaled** to the scale of the postsynaptic previous the synapse to make them compatible.

Because it does not wait for real time, this implementation is extremely fast. The optimization process typically requires **400 evaluations** (**50 initial samples** and **350 BO iterations**) to reach a good set of parameters, completing in around **30 seconds**.

### Compilation

You can compile the project using CMake:

```bash
mkdir build
cd build
cmake ..
make
```

### Executables

- `bo`: Main executable that runs the Bayesian Optimization process.
  - **Usage**: `./bo <csv_path> <column_idx> <csv_step (ms)> <start_time (ms)> <stabilization_time (ms)> <evaluation_time (ms)> <observation_time (ms)> <initial_samples> <iterations> <search_phase> <check_drift> <syn_model_step_factor> <syn_component> <cutoff_frequency (kHz)> <expected_i_min> <expected_i_max> <i_min> <i_max> <verbose> <output_yaml> <jsonl_history_file_path>`
  - **Arguments**:

| Argument | Description |
|---|---|
| `<csv_path>` | Path to the CSV file containing the presynaptic voltage recording. |
| `<column_idx>` | Column index of the voltage in the CSV (0-indexed). |
| `<csv_step (ms)>` | Temporal step size between CSV samples in milliseconds. |
| `<start_time (ms)>` | Time in milliseconds to start reading the recording from. |
| `<stabilization_time (ms)>` | Initial time to let the model stabilize before evaluating. |
| `<evaluation_time (ms)>` | Total duration in milliseconds for the synaptic evaluation. |
| `<observation_time (ms)>` | Time used to compute the base statistics of the signal for the scaling. |
| `<initial_samples>` | Number of initial Latin Hypercube Sampling (LHS) samples. |
| `<iterations>` | Number of Bayesian Optimization iterations to run after the initial samples. |
| `<search_phase>` | 1 to search for phase (excitatory) coupling, 0 for antiphase (inhibitory). |
| `<check_drift>` | 1 to enable signal drift correction, 0 to disable. |
| `<syn_model_step_factor>` | Sub-stepping factor for the synaptic model integration. |
| `<syn_component>` | Synaptic component to optimize (0 = fast, 1 = slow, 2 = both). |
| `<cutoff_frequency (kHz)>` | Cutoff frequency for the Butterworth low-pass filter in kHz to separate spikes from the slow wave in the presynaptic potential. |
| `<expected_i_min>` / `<expected_i_max>` | Expected minimum and maximum current bounds in nA. |
| `<i_min>` / `<i_max>` | Absolute clamping limits for the output current in nA. |
| `<verbose>` | 1 to enable verbose console output, 0 otherwise. |
| `<output_yaml>` | Path to save the best parameters found in YAML format. |
| `[jsonl_history_file_path]` | (Optional) Path to export the optimization history in JSONL format. |
- `consts_calculator`: Auxiliary executable that precalculates the dynamical constants of the neuronal model (minimum, maximum, and average points per burst for each integration step size). It prints to standard output in C++ format to be directly copied into the `utils.hpp` file. Executed only once during development or when neuronal parameters change.

## Main Requirements

Third-party requirements are listed below:

| Requirement | Description | Version | Installation |
|---|---|---|---|
| [Neun](https://github.com/GNB-UAM/Neun/tree/437c5b95898a0542f2c96a5344856d8dac12ad52) | Neuronal and synapse model library to build dynamical systems | commit `437c5b9` | Automatic |
| [Limbo](https://github.com/resibots/limbo/releases/tag/v2.1.0) | Bayesian Optimization and Gaussian Processes library | v2.1.0 | Automatic |
| [KFR](https://github.com/kfrlib/kfr/releases/tag/7.0.1) | Digital signal processing library (low-pass filters like Butterworth, vectorized operations) | 7.0.1 | Automatic |
| [nlohmann/json](https://github.com/nlohmann/json/releases/tag/v3.11.3) | JSON serialization library (for the BO history) | v3.11.3 | Automatic |
| [yaml-cpp](https://github.com/jbeder/yaml-cpp/releases/tag/yaml-cpp-0.9.0) | YAML file I/O library (for the result parameters) | 0.9.0 | Automatic |
| [NLopt](https://github.com/stevengj/nlopt/releases/tag/v2.7.1) | Optimization library (Subplex) | v2.7.1 | **Pre-installation required** |
| [Eigen3](https://eigen.tuxfamily.org) | C++ template library for linear algebra (required by Limbo) | | **Pre-installation required** |
| [Boost](https://www.boost.org) | C++ libraries (components: system, filesystem, thread) | | **Pre-installation required** |
| [TBB](https://github.com/oneapi-src/oneTBB) | Intel Threading Building Blocks for parallel execution | | **Pre-installation required** |

## Related Repositories

| Repository | Description |
|---|---|
| [tfg-data](https://github.com/marcogomezhernandez2004/tfg-data) | Data, plots, and scripts related to them. |
| [RTHybrid](https://github.com/GNB-UAM/RTHybrid/tree/a83071dcb4ac85f85b7beb2f7b7b5f68e785db22) (commit `a83071d`) | Standalone real-time neuronal model program by GNB. The presynaptic signal scaling algorithm is based on its implementation. |

## References

[1] J. Golowasch, M. Casey, L. F. Abbott, and E. Marder, “Network stability from activity-dependent regulation of neuronal conductances,” Neural Computation, vol. 11, pp. 1079–1096, 07 1999.

[2] J. L. Hindmarsh and R. M. Rose, “A model of neuronal bursting using three coupled first order differential equations,” Proceedings of the Royal Society of London. B. Biological Sciences, vol. 221, pp. 87–102, 03 1984.
