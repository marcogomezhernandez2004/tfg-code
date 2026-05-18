#include <algorithm>
#include <cmath>
using namespace kfr;

namespace EvaluationPrivateConfig
{
    // Factor de longitud de padding para filtfilt (proporción sobre fs/fc)
    static constexpr double PAD_LEN_FACTOR = 1.5;

    // Orden del filtro Butterworth pasa-bajos para separar i_fast de i_slow
    static constexpr int BUTTERWORTH_ORDER = 4;

    // Pesos relativos de i_fast e i_slow en la puntuación final (se promedian)
    static constexpr double I_FAST_WEIGHT = 0.5;
    static constexpr double I_SLOW_WEIGHT = 0.5;

    // Factor de reducción de la distancia máxima de rango esperada
    static constexpr double EXPECTED_I_MARGIN_RANGE_FACTOR = 0.8;
}

namespace EvaluationPrivateConstants
{
    // Valor inicial de m_slow (variable de gating lento) al resetear la sinapsis
    static constexpr double M_SLOW_INITIAL_VALUE = 0.0;
}

// Calcula la distancia máxima admisible para normalizar el error de rango.
// Extiende el rango esperado con un margen proporcional y aplica un factor de reducción
static double calculate_expected_i_dist_max(double expected_i_min,
                                            double expected_i_max)
{

    const double range = expected_i_max - expected_i_min;
    return (range + (range * BOPublicConfig::EXPECTED_I_MARGIN_FACTOR * 2.0)) *
           EvaluationPrivateConfig::EXPECTED_I_MARGIN_RANGE_FACTOR;
}

// Reescala un valor del rango fuente [src_min, src_min+src_range] al destino [dst_min, dst_min+dst_range]
static double rescale_to_target(double value,
                                double src_min, double src_range,
                                double dst_min, double dst_range)
{

    const double norm = (value - src_min) / safe_divisor(src_range);
    return dst_min + (norm * dst_range);
}

// Reescala sin offset: solo aplica la proporción entre rangos (para componentes centradas en cero)
static double rescale_to_target_no_offset(double value,
                                          double src_range,
                                          double dst_range)
{

    const double norm = value / safe_divisor(src_range);
    return norm * dst_range;
}

// Calcula la puntuación de forma mediante correlación de Pearson.
// Centra la señal, calcula r, lo normaliza a [0,1].
// Si search_phase=true invierte (busca anti-correlación -> fase)
static double pearson_score(univector<double> &sig,
                            const univector<double> &ref_sig_centered,
                            double ref_sig_factor,
                            bool search_phase)
{

    sig -= mean(sig); // Centra la señal candidata
    const double sig_factor = std::sqrt(sum(sqr(sig)));
    // Coeficiente de Pearson r in [-1, 1]
    const double r = sum(sig * ref_sig_centered) / safe_divisor(sig_factor * ref_sig_factor);
    // Normaliza a [0, 1]: (r+1)/2
    const double normalized = (r + 1.0) / 2.0;
    // En phase: queremos anti-correlación (devuelve 1 - normalized)
    return search_phase ? 1.0 - normalized : normalized;
}

// Calcula la puntuación de rango: mide cuánto se desvían los extremos observados de los esperados.
// Error = media de |obs_min - exp_min| y |obs_max - exp_max|, normalizado por dist_max
static double range_score(double observed_min, double observed_max,
                          double expected_min, double expected_max,
                          double pts_dist_max)
{

    const double error = (std::abs(observed_min - expected_min) +
                          std::abs(observed_max - expected_max)) *
                         0.5;
    const double normalized_error = error / safe_divisor(pts_dist_max);
    return 1.0 - normalized_error;
}

// Precalcula las señales de referencia y los rangos esperados de corriente.
// Filtra v_pre con Butterworth pasa-bajos para obtener:
//   - ref_i_fast = v_pre - filtrado (componente rápida = residuo de alta frecuencia)
//   - ref_i_slow = filtrado (componente lenta = baja frecuencia)
// Si se usan ambas componentes, los rangos se reescalan proporcionalmente.
ConstantEvaluationVals calc_constant_evaluation_vals(
    const univector_ref<const double> &v_pre_sig,
    double v_pre_min,
    double v_pre_max,
    double csv_step,
    double fc,
    bool use_i_fast,
    bool use_i_slow,
    double expected_i_min,
    double expected_i_max,
    bool search_phase)
{

    ConstantEvaluationVals result;

    const double fs = 1.0 / safe_divisor(csv_step); // Frecuencia de muestreo
    const size_t use_size = v_pre_sig.size();

    // Padding simétrico para filtfilt (reflexión en los bordes para evitar artefactos)
    const size_t effective_pad = std::min(use_size - 1,
                                          static_cast<size_t>(EvaluationPrivateConfig::PAD_LEN_FACTOR * fs / safe_divisor(fc)));
    univector<double> padded(use_size + (2 * effective_pad));

    // Copia la señal al centro del buffer con padding
    univector_ref<double> padded_seg = padded.slice(effective_pad, use_size);
    process(padded_seg, v_pre_sig);

    // Relleno reflejado en los bordes (como scipy.signal.filtfilt con padtype='odd')
    double *padded_ptr = padded.data();
    const double *v_pre_sig_ptr = v_pre_sig.data();
    const double left_edge_x2 = 2.0 * v_pre_sig_ptr[0];
    const double right_edge_x2 = 2.0 * v_pre_sig_ptr[use_size - 1];
    for (size_t i = 0; i < effective_pad; i++)
    {

        padded_ptr[effective_pad - 1 - i] = left_edge_x2 - v_pre_sig_ptr[i + 1];
        padded_ptr[use_size + effective_pad + i] = right_edge_x2 - v_pre_sig_ptr[use_size - 2 - i];
    }

    // Aplica filtro Butterworth pasa-bajos bidireccional (zero-phase) a fc frecuencia de corte
    filtfilt(padded, to_sos<double>(iir_lowpass(
                         butterworth(EvaluationPrivateConfig::BUTTERWORTH_ORDER), fc, fs)));

    const bool use_both = use_i_fast && use_i_slow;
    double i_dist_max = 0.0;
    double v_pre_range = 0.0, expected_i_range = 0.0;
    if (use_both)
    {
        // Cuando se usan ambas, los rangos de cada componente se reescalan proporcionalmente
        v_pre_range = v_pre_max - v_pre_min;
        expected_i_range = expected_i_max - expected_i_min;
    }
    else
    {
        // Con una sola componente, el rango esperado es el global
        i_dist_max = calculate_expected_i_dist_max(expected_i_min, expected_i_max);
    }

    if (use_i_fast)
    {
        // Referencia i_fast = v_pre original - componente filtrada (residuo de alta frecuencia)
        univector<double> &ref_i_fast_sig_centered = result.ref_i_fast_sig_centered;
        ref_i_fast_sig_centered = v_pre_sig - padded_seg;

        double &ref_i_fast_min = result.ref_i_fast_min;
        double &ref_i_fast_max = result.ref_i_fast_max;
        double &i_fast_dist_max = result.i_fast_dist_max;
        if (use_both)
        {
            // Reescala los extremos de la referencia al rango esperado de corriente
            // En phase: la corriente va invertida respecto al voltaje
            if (search_phase)
            {
                ref_i_fast_min = rescale_to_target_no_offset(-maxof(ref_i_fast_sig_centered), v_pre_range, expected_i_range);
                ref_i_fast_max = rescale_to_target_no_offset(-minof(ref_i_fast_sig_centered), v_pre_range, expected_i_range);
            }
            else
            {
                ref_i_fast_min = rescale_to_target_no_offset(minof(ref_i_fast_sig_centered), v_pre_range, expected_i_range);
                ref_i_fast_max = rescale_to_target_no_offset(maxof(ref_i_fast_sig_centered), v_pre_range, expected_i_range);
            }
            i_fast_dist_max = calculate_expected_i_dist_max(ref_i_fast_min, ref_i_fast_max);
        }
        else
        {

            ref_i_fast_min = expected_i_min;
            ref_i_fast_max = expected_i_max;
            i_fast_dist_max = i_dist_max;
        }

        // Centra y precalcula el factor de Pearson para la referencia i_fast
        ref_i_fast_sig_centered -= mean(ref_i_fast_sig_centered);
        result.ref_i_fast_sig_factor = std::sqrt(sum(sqr(ref_i_fast_sig_centered)));
    }

    if (use_i_slow)
    {
        // Referencia i_slow = componente filtrada (baja frecuencia de v_pre)
        univector<double> &ref_i_slow_sig_centered = result.ref_i_slow_sig_centered;
        ref_i_slow_sig_centered = padded_seg;

        double &ref_i_slow_min = result.ref_i_slow_min;
        double &ref_i_slow_max = result.ref_i_slow_max;
        double &i_slow_dist_max = result.i_slow_dist_max;
        if (use_both)
        {
            // Reescala los extremos de la referencia lenta al rango esperado
            if (search_phase)
            {
                // En phase: invierte y reescala desde -v_pre_max
                const double v_pre_min_to_use = -v_pre_max;
                ref_i_slow_min = rescale_to_target(-maxof(ref_i_slow_sig_centered), v_pre_min_to_use, v_pre_range,
                                                   expected_i_min, expected_i_range);
                ref_i_slow_max = rescale_to_target(-minof(ref_i_slow_sig_centered), v_pre_min_to_use, v_pre_range,
                                                   expected_i_min, expected_i_range);
            }
            else
            {
                ref_i_slow_min = rescale_to_target(minof(ref_i_slow_sig_centered), v_pre_min, v_pre_range,
                                                   expected_i_min, expected_i_range);
                ref_i_slow_max = rescale_to_target(maxof(ref_i_slow_sig_centered), v_pre_min, v_pre_range,
                                                   expected_i_min, expected_i_range);
            }
            i_slow_dist_max = calculate_expected_i_dist_max(ref_i_slow_min, ref_i_slow_max);
        }
        else
        {

            ref_i_slow_min = expected_i_min;
            ref_i_slow_max = expected_i_max;
            i_slow_dist_max = i_dist_max;
        }

        // Centra y precalcula el factor de Pearson para la referencia i_slow
        ref_i_slow_sig_centered -= mean(ref_i_slow_sig_centered);
        result.ref_i_slow_sig_factor = std::sqrt(sum(sqr(ref_i_slow_sig_centered)));
    }

    return result;
}

// Evalúa un candidato sináptico:
// 1. Configura la sinapsis con los parámetros del candidato
// 2. Simula paso a paso: la neurona presináptica inyecta v_pre de la señal escalada,
//    la postsináptica se integra con la corriente sináptica como entrada
// 3. Recoge i_fast e i_slow generados
// 4. Calcula las puntuaciones de rango y forma (Pearson) ponderadas
template <typename Integrator, typename NeuronType,
          ResetStateFunc<NeuronType> ResetStateFuncType,
          GetVFunc<NeuronType> GetVFuncType>
ChemicalSynapseEvaluation evaluate_candidate(
    const ChemicalSynapseParams &candidate,
    ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapse,
    NeuronType &model_neur,
    const ScaledSigResult &scaled_result,
    EvaluationISigBuffers &buffers,
    size_t v_pre_sig_start_idx,
    const ConstantEvaluationVals &constant_evaluation_vals,
    bool use_i_fast,
    bool use_i_slow,
    bool search_phase,
    double i_min,
    double i_max,
    ResetStateFuncType reset_state_neur,
    GetVFuncType get_v_neur)
{

    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    const ChemicalSynapseParams &params = candidate;

    // Configura los parámetros de la sinapsis con los valores del candidato
    synapse.set(ChemicalSynapsisType::Esyn, params.e_syn);

    double *i_fast_sig_ptr = nullptr;
    if (use_i_fast)
    {
        i_fast_sig_ptr = buffers.i_fast_sig.data();

        synapse.set(ChemicalSynapsisType::gfast, params.g_fast);
        synapse.set(ChemicalSynapsisType::sfast, params.s_fast);
        synapse.set(ChemicalSynapsisType::Vfast, params.v_fast);
    }

    double *i_slow_sig_ptr = nullptr;
    if (use_i_slow)
    {
        i_slow_sig_ptr = buffers.i_slow_sig.data();

        synapse.set(ChemicalSynapsisType::gslow, params.g_slow);
        synapse.set(ChemicalSynapsisType::Vslow, params.v_slow);
        synapse.set(ChemicalSynapsisType::sslow, params.s_slow);
        synapse.set(ChemicalSynapsisType::k1, params.k1);
        synapse.set(ChemicalSynapsisType::k2, params.k2);
    }

    // Resetea m_slow a 0 y el estado de la neurona postsináptica
    synapse.set(ChemicalSynapsisType::mslow, EvaluationPrivateConstants::M_SLOW_INITIAL_VALUE);

    reset_state_neur(model_neur);

    constexpr auto i_enum = ChemicalSynapsisType::i;
    constexpr auto i_fast_enum = ChemicalSynapsisType::ifast;
    constexpr auto i_slow_enum = ChemicalSynapsisType::islow;

    const size_t total_size = scaled_result.sig.size();
    const size_t points_factor = scaled_result.points_factor;
    const double dt = scaled_result.dt;
    const double *v_pre_sig_ptr = scaled_result.sig.data();
    const double *interpolated_points_ptr = scaled_result.interpolated_points.data();

    size_t interp_pts_counter = 0;
    size_t v_pre_sig_idx = 0;

    // Fase de estabilización: avanza la simulación sin registrar las corrientes
    // para que el estado de la sinapsis y neurona postsináptica se estabilicen
    for (; v_pre_sig_idx < v_pre_sig_start_idx; v_pre_sig_idx++)
    {
        synapse.step(dt, v_pre_sig_ptr[v_pre_sig_idx], get_v_neur(model_neur));

        const double i_val = std::clamp(synapse.get(i_enum), i_min, i_max);

        // Inyecta -i (convención de corriente: entrada negativa a la neurona)
        model_neur.add_synaptic_input(-i_val);
        model_neur.step(dt);
        // Sub-pasos de interpolación dentro de cada muestra
        for (size_t k = 1; k < points_factor; k++, interp_pts_counter++)
        {

            synapse.step(dt, interpolated_points_ptr[interp_pts_counter], get_v_neur(model_neur));
            const double i_interp_val = std::clamp(synapse.get(i_enum), i_min, i_max);
            model_neur.add_synaptic_input(-i_interp_val);
            model_neur.step(dt);
        }
    }

    // Fase de evaluación: avanza y registra las corrientes i_fast e i_slow
    size_t syn_sig_idx = 0;

    for (; v_pre_sig_idx < total_size - 1; v_pre_sig_idx++, syn_sig_idx++)
    {
        synapse.step(dt, v_pre_sig_ptr[v_pre_sig_idx], get_v_neur(model_neur));
        const double i_val = std::clamp(synapse.get(i_enum), i_min, i_max);
        model_neur.add_synaptic_input(-i_val);
        model_neur.step(dt);

        // Registra las corrientes sinápticas de este paso
        if (use_i_fast)
            i_fast_sig_ptr[syn_sig_idx] = synapse.get(i_fast_enum);
        if (use_i_slow)
            i_slow_sig_ptr[syn_sig_idx] = synapse.get(i_slow_enum);
        // Sub-pasos interpolados (no se registran, solo se simula para mantener estado)
        for (size_t k = 1; k < points_factor; k++, interp_pts_counter++)
        {
            synapse.step(dt, interpolated_points_ptr[interp_pts_counter], get_v_neur(model_neur));
            const double i_interp_val = std::clamp(synapse.get(i_enum), i_min, i_max);
            model_neur.add_synaptic_input(-i_interp_val);
            model_neur.step(dt);
        }
    }

    // Último punto (sin sub-pasos posteriores)
    synapse.step(dt, v_pre_sig_ptr[v_pre_sig_idx], get_v_neur(model_neur));
    const double i_val = std::clamp(synapse.get(i_enum), i_min, i_max);
    model_neur.add_synaptic_input(-i_val);
    model_neur.step(dt);
    if (use_i_fast)
        i_fast_sig_ptr[syn_sig_idx] = synapse.get(i_fast_enum);
    if (use_i_slow)
        i_slow_sig_ptr[syn_sig_idx] = synapse.get(i_slow_enum);

    // --- Cálculo de puntuaciones ---
    constexpr double I_FAST_WEIGHT = EvaluationPrivateConfig::I_FAST_WEIGHT;
    constexpr double I_SLOW_WEIGHT = EvaluationPrivateConfig::I_SLOW_WEIGHT;

    double i_range_score_accum = 0.0;
    double i_shape_score_accum = 0.0;
    double total_weight = 0.0;

    if (use_i_fast)
    {

        univector<double> &i_fast_sig = buffers.i_fast_sig;
        // Puntuación de rango: compara [min,max] observados con los esperados
        const double i_fast_range_score = range_score(
            minof(i_fast_sig),
            maxof(i_fast_sig),
            constant_evaluation_vals.ref_i_fast_min,
            constant_evaluation_vals.ref_i_fast_max,
            constant_evaluation_vals.i_fast_dist_max);

        // Puntuación de forma: correlación de Pearson con la referencia
        const double i_fast_shape_score = pearson_score(
            i_fast_sig,
            constant_evaluation_vals.ref_i_fast_sig_centered,
            constant_evaluation_vals.ref_i_fast_sig_factor,
            search_phase);

        i_range_score_accum += I_FAST_WEIGHT * i_fast_range_score;
        i_shape_score_accum += I_FAST_WEIGHT * i_fast_shape_score;
        total_weight += I_FAST_WEIGHT;
    }

    if (use_i_slow)
    {

        univector<double> &i_slow_sig = buffers.i_slow_sig;
        const double i_slow_range_score = range_score(
            minof(i_slow_sig),
            maxof(i_slow_sig),
            constant_evaluation_vals.ref_i_slow_min,
            constant_evaluation_vals.ref_i_slow_max,
            constant_evaluation_vals.i_slow_dist_max);

        const double i_slow_shape_score = pearson_score(
            i_slow_sig,
            constant_evaluation_vals.ref_i_slow_sig_centered,
            constant_evaluation_vals.ref_i_slow_sig_factor,
            search_phase);

        i_range_score_accum += I_SLOW_WEIGHT * i_slow_range_score;
        i_shape_score_accum += I_SLOW_WEIGHT * i_slow_shape_score;
        total_weight += I_SLOW_WEIGHT;
    }

    // Promedio ponderado de ambas componentes
    return ChemicalSynapseEvaluation(i_range_score_accum / total_weight,
                                     i_shape_score_accum / total_weight);
}
