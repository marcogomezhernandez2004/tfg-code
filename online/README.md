### Bidirectional Chemical Synapse BO

**Requirements:** libraries KFR DSP (compiled with -DCMAKE_POSITION_INDEPENDENT_CODE=ON), Limbo (included as a git submodule, run `git submodule update --init --recursive` from the root of the repository), nlohmann/json (placed in the `include/nlohmann/` directory, no need to do anything else), Eigen3, NLopt, and Boost (the last three are assumed to be installed in the system).
**Limitations:** The user should change the -lkfr_dsp_... flag from avx2 to the appropriate one for their system in the Makefile; and KFR DSP must be compiled with -DCMAKE_POSITION_INDEPENDENT_CODE=ON.

![Bidirectional Chemical Synapse BO GUI](bidirectional_chemical_synapse_BO.png)

<!--start-->
<p><b>Bidirectional Chemical Synapse BO</b><br>module for RTXI that implements a bidirectional chemical synapse model and runs Bayesian Optimization (BO) online to fit synaptic parameters based on captured voltage/current signals.</p>
<!--end-->

#### Input
1. input(0) - V_pre 1->2 (mV) : Presynaptic membrane potential 1->2
2. input(1) - V_post 1->2 (mV) : Postsynaptic membrane potential 1->2
3. input(2) - V_pre 2->1 (mV) : Presynaptic membrane potential 2->1
4. input(3) - V_post 2->1 (mV) : Postsynaptic membrane potential 2->1
5. input(4) - V_pre min 1->2 (mV) : Dynamic V_pre min 1->2
6. input(5) - V_pre max 1->2 (mV) : Dynamic V_pre max 1->2
7. input(6) - V_post min 1->2 (mV) : Dynamic V_post min 1->2
8. input(7) - V_post max 1->2 (mV) : Dynamic V_post max 1->2
9. input(8) - V_pre min 2->1 (mV) : Dynamic V_pre min 2->1
10. input(9) - V_pre max 2->1 (mV) : Dynamic V_pre max 2->1
11. input(10) - V_post min 2->1 (mV) : Dynamic V_post min 2->1
12. input(11) - V_post max 2->1 (mV) : Dynamic V_post max 2->1

#### Output
1. output(0) - Current 1->2 (nA) : Total synaptic current 1->2
2. output(1) - Current 2->1 (nA) : Total synaptic current 2->1

#### Parameters
1. BO initial samples - Number of initialization samples for BO
2. BO iterations - Number of BO iterations after initial sampling
3. BO evaluation time (ms) - Time to record signals per evaluation
4. BO stabilization time (ms) - Wait time after setting params before recording
5. BO search phase (1/0) - 1 = Enable, 0 = Disable
6. BO current min to achieve 1->2 (nA) - Target minimum current for direction 1->2
7. BO current max to achieve 1->2 (nA) - Target maximum current for direction 1->2
8. BO current min to achieve 2->1 (nA) - Target minimum current for direction 2->1
9. BO current max to achieve 2->1 (nA) - Target maximum current for direction 2->1
10. BO cutoff frequency 1 (kHz) - To separate the I_fast and I_slow for BO in synapse 1->2
11. BO cutoff frequency 2 (kHz) - To separate the I_fast and I_slow for BO in synapse 2->1
12. Dynamic V_pre min and max 1->2 (1/0) - 1 = Enable, 0 = Disable; necessary for BO
13. V_pre min 1->2 (mV) - Necessary for BO
14. V_pre max 1->2 (mV) - Necessary for BO
15. Dynamic V_post min and max 1->2 (1/0) - 1 = Enable, 0 = Disable; necessary for BO
16. V_post min 1->2 (mV) - Necessary for BO
17. V_post max 1->2 (mV) - Necessary for BO
18. Dynamic V_pre min and max 2->1 (1/0) - 1 = Enable, 0 = Disable; necessary for BO
19. V_pre min 2->1 (mV) - Necessary for BO
20. V_pre max 2->1 (mV) - Necessary for BO
21. Dynamic V_post min and max 2->1 (1/0) - 1 = Enable, 0 = Disable; necessary for BO
22. V_post min 2->1 (mV) - Necessary for BO
23. V_post max 2->1 (mV) - Necessary for BO
24. Current min 1->2 (nA) - Fixed output clamp min for current 1->2
25. Current max 1->2 (nA) - Fixed output clamp max for current 1->2
26. Current min 2->1 (nA) - Fixed output clamp min for current 2->1
27. Current max 2->1 (nA) - Fixed output clamp max for current 2->1
28. Verbose (1/0) - Enable/disable BO candidate evaluation logging
29. factor in dt (ms) = period (ms) * factor - Factor for calculating dt form the period; dt in ms
30. Use I_fast 1->2 (1/0) - 1 = Enable, 0 = Disable
31. Use I_slow 1->2 (1/0) - 1 = Enable, 0 = Disable
32. Use I_fast 2->1 (1/0) - 1 = Enable, 0 = Disable
33. Use I_slow 2->1 (1/0) - 1 = Enable, 0 = Disable
34. E_syn 1->2 (mV) - Synaptic reversal potential (1->2)
35. g_fast 1->2 (uS) - Fast conductance (1->2)
36. s_fast 1->2 (1/mV) - Fast sigmoid slope (1->2)
37. V_fast 1->2 (mV) - Fast sigmoid threshold (1->2)
38. g_slow 1->2 (uS) - Slow conductance (1->2)
39. k1 1->2 (1/ms) - Slow gating opening rate (1->2)
40. k2 1->2 (1/ms) - Slow gating closing rate (1->2)
41. s_slow 1->2 (1/mV) - Slow sigmoid slope (1->2)
42. V_slow 1->2 (mV) - Slow sigmoid threshold (1->2)
43. E_syn 2->1 (mV) - Synaptic reversal potential (2->1)
44. g_fast 2->1 (uS) - Fast conductance (2->1)
45. s_fast 2->1 (1/mV) - Fast sigmoid slope (2->1)
46. V_fast 2->1 (mV) - Fast sigmoid threshold (2->1)
47. g_slow 2->1 (uS) - Slow conductance (2->1)
48. k1 2->1 (1/ms) - Slow gating opening rate (2->1)
49. k2 2->1 (1/ms) - Slow gating closing rate (2->1)
50. s_slow 2->1 (1/mV) - Slow sigmoid slope (2->1)
51. V_slow 2->1 (mV) - Slow sigmoid threshold (2->1)

#### States
1. BO evaluations completed - Finishes when this is initial samples + iterations
2. I_fast 1->2 (nA) - Fast synaptic current component 1->2
3. I_slow 1->2 (nA) - Slow synaptic current component 1->2
4. I_fast 2->1 (nA) - Fast synaptic current component 2->1
5. I_slow 2->1 (nA) - Slow synaptic current component 2->1

