# Chemical synapse model parameterization for biohybrid neural circuits

This repository contains the source code developed for the Bachelor's Thesis (TFG) *"Parametrización del modelo de sinapsis química para circuitos neuronales biohíbridos"* (Chemical synapse model parameterization for biohybrid neural circuits).

The goal is to automate the parameterization of the **graded chemical synapse model by Golowasch et al.** [1] for biohybrid circuits, where biological neurons and computational models interact in real time. This model exhibits ***parameter sloppiness*** (multiple parameter combinations can yield similar outputs), making manual tuning difficult. It separates the synaptic current to enable modular and frequency-selective coupling:
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

Both implementations use a **Runge-Kutta integrator** for solving the differential equations of the computational neuronal and synaptic models.

The user can **choose which current components** are used (calculated) and optimized to maximize efficiency.

## Repository Structure

The repository is organized into the following main directories, each targeting a specific execution environment. The optimization process typically requires **400 evaluations** (**50 initial samples** and **350 BO iterations**) to reach a good set of parameters.

| Directory | Description |
|---|---|
| [`offline/`](offline/) | **Offline Implementation**: Operates on a presynaptic voltage recording (CSV) in a **unidirectional** setup with a neuronal model as the postsynaptic neuron. Because it does not wait for real time, this implementation is extremely fast, completing 400 evaluations in around **30 seconds**. |
| [`online/`](online/) | **Online Implementation**: An **RTXI** module (version 2.3) for real-time experimentation. Supports **bidirectional** interactions with two simultaneous synapses, compatible with live neurons and/or models. Since it runs in real time, it takes under **16 minutes** for the same number of evaluations. |
| [`rthybrid_hindmarsh_rose_1984_neuron_v2/`](rthybrid_hindmarsh_rose_1984_neuron_v2/) | RTXI module for the modified **Hindmarsh-Rose** neuronal model (**with *v_h* parameter**), corrected from the RTHybrid for RTXI repository of the GNB at UAM. |

> Each subdirectory contains its own self-contained `README.md` with specific descriptions and instructions.

## Main Requirements

Third-party requirements vary by version. The main ones are listed below (refer to each subdirectory's README for details):

| Requirement | Description | Version |
|---|---|---|
| [Neun](https://github.com/GNB-UAM/Neun/tree/437c5b95898a0542f2c96a5344856d8dac12ad52) | Neuronal and synapse model library to build dynamical systems (offline version) | commit `437c5b9` |
| [Limbo](https://github.com/resibots/limbo/releases/tag/v2.1.0) | Bayesian Optimization and Gaussian Processes library | v2.1.0 |
| [KFR](https://github.com/kfrlib/kfr/releases/tag/7.0.1) | Digital signal processing library (low-pass filters like Butterworth, vectorized operations) | 7.0.1 |
| [nlohmann/json](https://github.com/nlohmann/json/releases/tag/v3.11.3) | JSON serialization library (for the BO history) | v3.11.3 |
| [yaml-cpp](https://github.com/jbeder/yaml-cpp/releases/tag/yaml-cpp-0.9.0) | YAML file I/O library (offline version, for the result parameters) | 0.9.0 |
| [NLopt](https://github.com/stevengj/nlopt/releases/tag/v2.7.1) | Optimization library (Subplex) | v2.7.1 |
| [RTXI](https://github.com/RTXI/rtxi/releases/tag/v2.3) | Real-time experimentation platform (online version) | v2.3 |
| [Eigen3](https://eigen.tuxfamily.org) | C++ template library for linear algebra (required by Limbo) | |
| [Boost](https://www.boost.org) | C++ libraries (components: system, filesystem, thread) | |
| [TBB](https://github.com/oneapi-src/oneTBB) | Intel Threading Building Blocks for parallel execution | |

## Related Repositories

- **Project data**: [tfg-data](https://github.com/marcogomezhernandez2004/tfg-data) — Data, plots, and scripts related to them.
- **RTHybrid**: [RTHybrid](https://github.com/GNB-UAM/RTHybrid/tree/a83071dcb4ac85f85b7beb2f7b7b5f68e785db22) (commit `a83071d`) — Standalone real-time neuronal model program by GNB. The presynaptic signal scaling algorithm is based on its implementation.
- **RTHybrid for RTXI**: [rthybrid-for-rtxi](https://github.com/GNB-UAM/rthybrid-for-rtxi/tree/f13a015084b9819a02663833aaccbc8e86d36161) (commit `f13a015`) — RTXI modules that replicate RTHybrid functionality (neuronal models, burst analysis, amplitude scaling...), useful with the online synapse module.
- **mimic-signal**: [mimic-signal](https://github.com/RTXI/mimic-signal/tree/e15e26cef364575de9a0629591ac9e2218637226) (commit `e15e26c`) — RTXI module for applying gain and offset to a signal (for amplitude scaling), useful with the online synapse module.
- **plugin-template**: [plugin-template](https://github.com/RTXI/plugin-template/tree/b7b3a3b606cca17778cac8c5a4846d6df8a0733a) (commit `b7b3a3b`) — Official RTXI module template (online version). The legacy version of this template is used for compatibility with RTXI 2.3.

## Thesis Information

| | |
|---|---|
| **Author** | Marco Gómez Hernández |
| **Advisor** | Pablo Varona Martínez |
| **University** | Universidad Autónoma de Madrid (UAM) |
| **School** | Escuela Politécnica Superior (EPS) |
| **Degree** | Bachelor's in Computer Science (*Grado en Ingeniería Informática*) |
| **Research Group** | Grupo de Neurocomputación Biológica (GNB) |
| **Date** | May 2026 |

## License

Refer to the `LICENSE` files in each subdirectory.

## References

[1] J. Golowasch, M. Casey, L. F. Abbott, and E. Marder, “Network stability from activity-dependent regulation of neuronal conductances,” Neural Computation, vol. 11, pp. 1079–1096, 07 1999.
