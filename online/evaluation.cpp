#include <chrono>
#include <cmath>
#include <kfr/all.hpp>
#include "bidirectional_chemical_synapse_BO.h"
#include <concepts>
using namespace kfr;

namespace EvaluationPrivateConfig
{
    // Orden del filtro Butterworth pasa-bajos para separar i_fast de i_slow
    static constexpr int BUTTERWORTH_ORDER = 4;

    // Pesos relativos de i_fast e i_slow en la puntuación de cada dirección
    static constexpr double I_FAST_WEIGHT = 0.5;
    static constexpr double I_SLOW_WEIGHT = 0.5;
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

// Puntuación de forma mediante correlación de Pearson.
// Centra la señal, calcula r, normaliza a [0,1].
// Si search_phase: invierte el resultado (busca anti-correlación -> fase)
template <typename T>
    requires std::same_as<T, univector<double>> || std::same_as<T, univector_ref<double>>
static double pearson_score(univector<double> &sig,
                            const T &ref_sig_centered,
                            double ref_sig_factor,
                            unsigned int search_phase)
{

    sig -= mean(sig); // Centra la señal candidata
    const double sig_factor = std::sqrt(sum(sqr(sig)));
    // Coeficiente de Pearson r in [-1, 1]
    const double r = sum(sig * ref_sig_centered) / safe_divisor(sig_factor * ref_sig_factor);
    // Normaliza a [0,1]
    const double normalized = (r + 1.0) / 2.0;
    return search_phase ? 1.0 - normalized : normalized;
}

// Puntuación de rango: mide la desviación de los extremos observados respecto a los esperados
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

// Evalúa las señales de una sola dirección sináptica (1->2 o 2->1):
// 1. Filtra v_pre con Butterworth pasa-bajos (filtfilt con padding reflectivo)
// 2. Obtiene referencia i_fast (residuo alta frecuencia) e i_slow (componente filtrada)
// 3. Si se usan ambas, reescala los rangos esperados proporcionalmente
// 4. Calcula puntuaciones de rango y forma (Pearson) ponderadas por componente
static ChemicalSynapseEvaluation evaluate_sigs_one_direction(
    univector<double> &v_pre_sig,
    univector<double> &i_fast_sig,
    univector<double> &i_slow_sig,
    size_t effective_pad,
    double fs,
    double fc,
    unsigned int use_i_fast,
    unsigned int use_i_slow,
    unsigned int search_phase,
    double expected_i_min,
    double expected_i_max,
    double i_dist_max,
    univector<double> &padded_buff)
{

    const size_t use_size = v_pre_sig.size();

    // Copia v_pre al centro del buffer con padding
    univector_ref<double> padded_seg = padded_buff.slice(effective_pad, use_size);
    process(padded_seg, v_pre_sig);

    // Relleno reflejado en los bordes (como scipy filtfilt padtype='odd')
    double *padded_ptr = padded_buff.data();
    const double *v_pre_sig_ptr = v_pre_sig.data();
    const double left_edge_x2 = 2.0 * v_pre_sig_ptr[0];
    const double right_edge_x2 = 2.0 * v_pre_sig_ptr[use_size - 1];
    for (size_t i = 0; i < effective_pad; i++)
    {

        padded_ptr[effective_pad - 1 - i] = left_edge_x2 - v_pre_sig_ptr[i + 1];
        padded_ptr[use_size + effective_pad + i] = right_edge_x2 - v_pre_sig_ptr[use_size - 2 - i];
    }

    // Aplica filtro Butterworth pasa-bajos bidireccional (zero-phase)
    filtfilt(padded_buff, to_sos<double>(iir_lowpass(
                              butterworth(EvaluationPrivateConfig::BUTTERWORTH_ORDER), fc, fs)));

    const bool use_both = use_i_fast && use_i_slow;

    constexpr double I_FAST_WEIGHT = EvaluationPrivateConfig::I_FAST_WEIGHT;
    constexpr double I_SLOW_WEIGHT = EvaluationPrivateConfig::I_SLOW_WEIGHT;

    double i_range_score_accum = 0.0;
    double i_shape_score_accum = 0.0;
    double total_weight = 0.0;

    double v_pre_min = 0.0, v_pre_max = 0.0, v_pre_range = 0.0, expected_i_range = 0.0;
    if (use_both)
    {
        // Cuando se usan ambas, los rangos de cada componente se reescalan proporcionalmente
        v_pre_min = minof(v_pre_sig);
        v_pre_max = maxof(v_pre_sig);
        v_pre_range = v_pre_max - v_pre_min;
        expected_i_range = expected_i_max - expected_i_min;
    }

    if (use_i_fast)
    {
        // Referencia i_fast = v_pre original - filtrado (residuo de alta frecuencia)
        // Reutiliza v_pre_sig como buffer para la referencia fast
        univector<double> &ref_i_fast_sig_aux = v_pre_sig;
        ref_i_fast_sig_aux -= padded_seg;

        double ref_i_fast_min, ref_i_fast_max, i_fast_dist_max;
        if (use_both)
        {
            // Reescala los extremos al rango esperado de corriente
            if (search_phase)
            {
                ref_i_fast_min = rescale_to_target_no_offset(-maxof(ref_i_fast_sig_aux), v_pre_range, expected_i_range);
                ref_i_fast_max = rescale_to_target_no_offset(-minof(ref_i_fast_sig_aux), v_pre_range, expected_i_range);
            }
            else
            {
                ref_i_fast_min = rescale_to_target_no_offset(minof(ref_i_fast_sig_aux), v_pre_range, expected_i_range);
                ref_i_fast_max = rescale_to_target_no_offset(maxof(ref_i_fast_sig_aux), v_pre_range, expected_i_range);
            }
            i_fast_dist_max = calculate_expected_i_dist_max(ref_i_fast_min, ref_i_fast_max);
        }
        else
        {
            ref_i_fast_min = expected_i_min;
            ref_i_fast_max = expected_i_max;
            i_fast_dist_max = i_dist_max;
        }

        // Centra y calcula el factor de Pearson para la referencia i_fast
        ref_i_fast_sig_aux -= mean(ref_i_fast_sig_aux);
        const double ref_i_fast_factor = std::sqrt(sum(sqr(ref_i_fast_sig_aux)));

        const double i_fast_range_score = range_score(
            minof(i_fast_sig),
            maxof(i_fast_sig),
            ref_i_fast_min,
            ref_i_fast_max,
            i_fast_dist_max);

        const double i_fast_shape_score =
            pearson_score(i_fast_sig,
                          ref_i_fast_sig_aux,
                          ref_i_fast_factor,
                          search_phase);

        i_range_score_accum += I_FAST_WEIGHT * i_fast_range_score;
        i_shape_score_accum += I_FAST_WEIGHT * i_fast_shape_score;
        total_weight += I_FAST_WEIGHT;
    }

    if (use_i_slow)
    {
        // Referencia i_slow = componente filtrada (baja frecuencia de v_pre)
        // Reutiliza padded_seg como buffer para la referencia slow
        univector_ref<double> &ref_i_slow_sig_aux = padded_seg;

        double ref_i_slow_min, ref_i_slow_max, i_slow_dist_max;
        if (use_both)
        {
            // Reescala los extremos de la referencia lenta al rango esperado
            if (search_phase)
            {
                const double v_pre_min_to_use = -v_pre_max;
                ref_i_slow_min = rescale_to_target(-maxof(ref_i_slow_sig_aux), v_pre_min_to_use, v_pre_range,
                                                   expected_i_min, expected_i_range);
                ref_i_slow_max = rescale_to_target(-minof(ref_i_slow_sig_aux), v_pre_min_to_use, v_pre_range,
                                                   expected_i_min, expected_i_range);
            }
            else
            {
                ref_i_slow_min = rescale_to_target(minof(ref_i_slow_sig_aux), v_pre_min, v_pre_range,
                                                   expected_i_min, expected_i_range);
                ref_i_slow_max = rescale_to_target(maxof(ref_i_slow_sig_aux), v_pre_min, v_pre_range,
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

        // Centra y calcula el factor de Pearson para la referencia i_slow
        ref_i_slow_sig_aux -= mean(ref_i_slow_sig_aux);
        const double ref_i_slow_factor = std::sqrt(sum(sqr(ref_i_slow_sig_aux)));

        const double i_slow_range_score = range_score(
            minof(i_slow_sig),
            maxof(i_slow_sig),
            ref_i_slow_min,
            ref_i_slow_max,
            i_slow_dist_max);

        const double i_slow_shape_score =
            pearson_score(i_slow_sig,
                          ref_i_slow_sig_aux,
                          ref_i_slow_factor,
                          search_phase);

        i_range_score_accum += I_SLOW_WEIGHT * i_slow_range_score;
        i_shape_score_accum += I_SLOW_WEIGHT * i_slow_shape_score;
        total_weight += I_SLOW_WEIGHT;
    }

    // Promedio ponderado de ambas componentes
    return ChemicalSynapseEvaluation(i_range_score_accum / total_weight,
                                     i_shape_score_accum / total_weight);
}

// ==========================================
//  evaluate_candidate: sincronización NRT -> RT y evaluación bidireccional
// ==========================================
// Este método se llama desde el hilo NRT (dentro de la BO) para cada candidato:
// 1. Comprueba si se ha solicitado parar
// 2. Espera a que el hilo RT haya leído el índice actual del double-buffer
// 3. Escribe los nuevos parámetros en el buffer alterno y los publica (release)
// 4. Espera la estabilización (sleep)
// 5. Activa la recogida de señales por el hilo RT (RT_storing = true)
// 6. Espera a que el hilo RT llene los buffers (RT_storing vuelve a false)
// 7. Evalúa las señales recogidas por dirección y promedia
ChemicalSynapseEvaluation BidirectionalChemicalSynapseBO::evaluate_candidate(
    const Candidate &candidate,
    double fs,
    size_t effective_pad_12,
    size_t effective_pad_21,
    EvaluationPadBuffers &pad_buffers,
    double i_dist_max_12,
    double i_dist_max_21,
    size_t &curr_synapse_idx)
{

    // Comprueba cancelación antes de empezar
    if (stop_BO.load(std::memory_order_relaxed))
        throw StopEvaluation();

    const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
    const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;

    if (!use_syn_12 && !use_syn_21)
    {

        QMetaObject::invokeMethod(this, "set_evaluations_completed", Qt::QueuedConnection,
                                  Q_ARG(double, evaluations_completed + 1));
        return ChemicalSynapseEvaluation{0.0, 0.0};
    }

    const std::chrono::duration<double, std::milli> active_wait_duration(BOPublicConfig::ACTIVE_WAIT_MS);

    // --- Paso 1: espera a que el hilo RT haya leído el índice actual ---
    // Necesario para no sobrescribir parámetros que el RT aún está usando
    if (!wait_until_RT_read_idx_or_stop(curr_synapse_idx))
        throw StopEvaluation();

    // --- Paso 2: escribe los nuevos parámetros en el buffer alterno ---
    const size_t new_synapse_idx = 1 - curr_synapse_idx;

    copy_selected_synapse_params(params_12[new_synapse_idx],
                                 candidate.params_12,
                                 use_i_fast_12,
                                 use_i_slow_12);

    copy_selected_synapse_params(params_21[new_synapse_idx],
                                 candidate.params_21,
                                 use_i_fast_21,
                                 use_i_slow_21);

    // Publica el nuevo índice con release (hace visibles las escrituras al hilo RT)
    synapse_idx.store(new_synapse_idx, std::memory_order_release);
    curr_synapse_idx = new_synapse_idx;

    // --- Paso 3: estabilización (espera antes de empezar a medir) ---
    if (stabilization_time > 0.0)
    {
        std::this_thread::sleep_for(std::chrono::duration<double, std::milli>(stabilization_time));
    }

    // --- Paso 4: activa la recogida de señales por el hilo RT ---
    storing_idx = 0;

    // Release: el hilo RT ve storing_idx = 0 y empieza a llenar los buffers
    RT_storing.store(true, std::memory_order_release);

    // --- Paso 5: espera a que el hilo RT llene todos los elementos ---
    // El hilo RT pone RT_storing = false cuando storing_idx >= num_elements
    while (RT_storing.load(std::memory_order_acquire))
    {
        if (stop_BO.load(std::memory_order_relaxed))
        {
            RT_storing.store(false, std::memory_order_relaxed);
            throw StopEvaluation();
        }
        std::this_thread::sleep_for(active_wait_duration);
    }

    // --- Paso 6: evalúa las señales recogidas por dirección ---
    double i_range_score_accum = 0.0;
    double i_shape_score_accum = 0.0;
    double num_directions = 0.0;

    if (use_syn_12)
    {

        ChemicalSynapseEvaluation score_12 = evaluate_sigs_one_direction(
            v_pre_sig_12,
            i_fast_sig_12,
            i_slow_sig_12,
            effective_pad_12,
            fs,
            fc_1,
            use_i_fast_12,
            use_i_slow_12,
            search_phase,
            expected_i_min_12, expected_i_max_12,
            i_dist_max_12,
            pad_buffers.padded_buff_12);
        i_range_score_accum += score_12.i_range_score;
        i_shape_score_accum += score_12.i_shape_score;
        num_directions++;
    }

    if (use_syn_21)
    {

        ChemicalSynapseEvaluation score_21 = evaluate_sigs_one_direction(
            v_pre_sig_21,
            i_fast_sig_21,
            i_slow_sig_21,
            effective_pad_21,
            fs,
            fc_2,
            use_i_fast_21,
            use_i_slow_21,
            search_phase,
            expected_i_min_21, expected_i_max_21,
            i_dist_max_21,
            pad_buffers.padded_buff_21);
        i_range_score_accum += score_21.i_range_score;
        i_shape_score_accum += score_21.i_shape_score;
        num_directions++;
    }

    // Promedia las puntuaciones de ambas direcciones
    const double i_range_score = i_range_score_accum / num_directions;
    const double i_shape_score = i_shape_score_accum / num_directions;
    QMetaObject::invokeMethod(this, "set_evaluations_completed", Qt::QueuedConnection,
                              Q_ARG(double, evaluations_completed + 1));

    return ChemicalSynapseEvaluation(i_range_score, i_shape_score);
}
