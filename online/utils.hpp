#ifndef EVALUATION_UTILS_H
#define EVALUATION_UTILS_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <Eigen/Core>
#include <kfr/all.hpp>

namespace BOPublicConfig
{
    // Factor de reducción de la distancia máxima de rango esperada
    inline constexpr double EXPECTED_I_MARGIN_RANGE_FACTOR = 0.8;

    // Tiempo de espera activa (ms) entre comprobaciones de sincronización NRT <-> RT
    inline constexpr double ACTIVE_WAIT_MS = 10.0;

    // Margen de voltaje añadido a los extremos de v_pre para v_fast y v_slow online
    // (compensa posible variabilidad de la señal real)
    inline constexpr double V_MARGIN_FACTOR = 0.125;

    // v_fast_min = v_pre_min + rango * V_FAST_MIN_FACTOR
    inline constexpr double V_FAST_MIN_FACTOR = 0.25;

    // Porcentaje de margen sobre el rango esperado de corriente para ampliar la zona de búsqueda de g
    inline constexpr double EXPECTED_I_MARGIN_FACTOR = 0.5;

    // s_fast in [S_FAST_MIN_FACTOR, S_FAST_MAX_FACTOR] / rango_v_pre
    inline constexpr double S_FAST_MIN_FACTOR = 4.12;
    inline constexpr double S_FAST_MAX_FACTOR = 13.72;
    // s_slow in [S_SLOW_MIN_FACTOR, S_SLOW_MAX_FACTOR] / rango_v_pre
    inline constexpr double S_SLOW_MIN_FACTOR = 1.72;
    inline constexpr double S_SLOW_MAX_FACTOR = 3.43;

    // k1 in [K1_MIN_FACTOR * fc_adim, K1_MAX_FACTOR * fc_adim] (proporcional a la frecuencia de corte adimensional)
    inline constexpr double K1_MAX_FACTOR = 29.845125;

    inline constexpr double K1_MIN_FACTOR = 29.845125e-4;

    // e_syn: distancia proporcional al rango de v_post (similar a offline)
    inline constexpr double E_SYN_FAR_TERM = 3.86;
    inline constexpr double E_SYN_NEAR_TERM = 0.2;

    // g_min = g_max * G_MIN_FACTOR
    inline constexpr double G_MIN_FACTOR = 0.001;

    // R = k2/k1 in [R_MIN, R_MAX]
    inline constexpr double R_MAX = 30.0;
    inline constexpr double R_MIN = 0.001;

    // Puntuación de rango muy mala; se asigna si el resultado no es finito
    inline constexpr double VERY_BAD_RANGE_SCORE = -1e6;
}

namespace BOPublicConstants
{
    // Valor mínimo positivo para log de valores cercanos a cero
    inline constexpr double SMALL_LOG = std::numeric_limits<double>::min();
    // Umbral mínimo para divisores; evita divisiones por cero
    inline constexpr double SMALL_DIVISOR = std::numeric_limits<double>::epsilon();
    inline constexpr double NEGATIVE_SMALL_DIVISOR = -SMALL_DIVISOR;
}

// Excepción lanzada cuando el hilo NRT detecta que se ha solicitado parar la BO
struct StopEvaluation
{
};

// Sigmoide sináptica: 1 / (1 + exp(s * (v_threshold - v_pre)))
inline double chemical_sigmoid(double s,
                               double v_threshold,
                               double v_pre)
{
    return 1.0 / (1.0 + std::exp(s * (v_threshold - v_pre)));
}

// Calcula la distancia máxima admisible para normalizar el error de rango
inline double calculate_expected_i_dist_max(double expected_i_min,
                                            double expected_i_max)
{
    const double range = expected_i_max - expected_i_min;
    return (range + (range * BOPublicConfig::EXPECTED_I_MARGIN_FACTOR * 2.0)) *
           BOPublicConfig::EXPECTED_I_MARGIN_RANGE_FACTOR;
}

// Divisor seguro: si |divisor| < epsilon, devuelve +-epsilon para evitar div/0
inline double safe_divisor(double divisor)
{

    return std::abs(divisor) < BOPublicConstants::SMALL_DIVISOR
               ? (divisor < 0.0 ? BOPublicConstants::NEGATIVE_SMALL_DIVISOR : BOPublicConstants::SMALL_DIVISOR)
               : divisor;
}

// Rangos de búsqueda para cada parámetro sináptico (igual que offline).
// Limbo trabaja en [0,1]^dim; estos rangos mapean [0,1] al espacio real.
struct BOParamRanges
{
    struct ParamRange
    {

        double min;   // Valor real mínimo
        double range; // max - min

        ParamRange() = default;

        ParamRange(double min, double max)
            : min(min), range(max - min)
        {
        }
    };

    ParamRange s_fast;
    ParamRange s_slow;
    ParamRange e_syn;
    ParamRange log_k1;     // En log-space
    ParamRange log_R;      // En log-space (R = k2/k1)
    ParamRange log_g_fast; // En log-space
    ParamRange log_g_slow; // En log-space
    ParamRange v_fast;
    ParamRange v_slow;

    // Calcula los rangos reales de cada parámetro sináptico (igual que offline + V_MARGIN_FACTOR en v_fast/v_slow)
    void init(double v_pre_min,
              double v_pre_max,
              double v_post_min,
              double v_post_max,
              double expected_i_min,
              double expected_i_max,
              unsigned int use_i_fast,
              unsigned int use_i_slow,
              unsigned int search_phase,
              double fc,
              double period_t,
              double dt)
    {
        constexpr double G_MIN_FACTOR = BOPublicConfig::G_MIN_FACTOR;
        constexpr double R_MIN = BOPublicConfig::R_MIN;
        constexpr double SMALL_LOG = BOPublicConstants::SMALL_LOG;
        constexpr double R_MAX = BOPublicConfig::R_MAX;

        const double v_pre_range = v_pre_max - v_pre_min;
        const double v_post_range = v_post_max - v_post_min;

        // e_syn: excitadora -> por encima de v_post_max; inhibidora -> por debajo de v_post_min
        const double e_syn_far_final_term = v_post_range * BOPublicConfig::E_SYN_FAR_TERM;
        const double e_syn_near_final_term = v_post_range * BOPublicConfig::E_SYN_NEAR_TERM;
        double e_syn_max = 0.0, e_syn_min = 0.0;
        if (search_phase)
        {
            // Excitadora: e_syn > v_post -> (v_post - e_syn) < 0 -> corriente entrante
            e_syn_min = v_post_max + e_syn_near_final_term;
            e_syn_max = v_post_max + e_syn_far_final_term;
        }
        else
        {
            // Inhibidora: e_syn < v_post -> (v_post - e_syn) > 0 -> corriente saliente
            e_syn_min = v_post_min - e_syn_far_final_term;
            e_syn_max = v_post_min - e_syn_near_final_term;
        }
        e_syn = ParamRange(e_syn_min, e_syn_max);

        // Margen de +-V_MARGIN_FACTOR sobre v_pre para compensar variabilidad de la señal biológica real
        const double v_pre_margin = v_pre_range * BOPublicConfig::V_MARGIN_FACTOR;
        const double v_pre_margin_min = v_pre_min - v_pre_margin;
        const double v_pre_margin_max = v_pre_max + v_pre_margin;

        // Margen sobre la corriente esperada (+-EXPECTED_I_MARGIN_FACTOR del rango) para ampliar la búsqueda de g
        const double expected_i_margin = (expected_i_max - expected_i_min) * BOPublicConfig::EXPECTED_I_MARGIN_FACTOR;
        const double expected_i_margin_min = expected_i_min - expected_i_margin;
        const double expected_i_margin_max = expected_i_max + expected_i_margin;
        const double safe_v_pre_range = safe_divisor(v_pre_range);

        if (use_i_fast)
        {
            // v_fast: umbral de la sigmoide rápida, cuarto superior de v_pre + margen online en el máximo
            v_fast = ParamRange(v_pre_min + (v_pre_range * BOPublicConfig::V_FAST_MIN_FACTOR),
                                v_pre_margin_max);
            // s_fast: pendiente de la sigmoide rápida, inversamente proporcional al rango de v_pre (1/V)
            const double s_fast_max = BOPublicConfig::S_FAST_MAX_FACTOR / safe_v_pre_range;
            s_fast = ParamRange(BOPublicConfig::S_FAST_MIN_FACTOR / safe_v_pre_range,
                                s_fast_max);

            // Sigmoide máxima: mayor pendiente (s_fast_max), umbral más bajo (v_fast.min), voltaje más alto (v_pre_max)
            const double sigmoid_fast_max = chemical_sigmoid(s_fast_max, v_fast.min, v_pre_max);
            // g_fast_max: invierte g = I / (sigma * (v_post - e_syn)) en los dos extremos de driving force (peor caso)
            const double g_fast_max = std::max(std::abs(expected_i_margin_max / safe_divisor((v_post_max - e_syn_min) * sigmoid_fast_max)),
                                               std::abs(expected_i_margin_min / safe_divisor((v_post_min - e_syn_max) * sigmoid_fast_max)));
            // g_fast_min = g_fast_max * G_MIN_FACTOR -> cubre ~3 órdenes de magnitud
            const double g_fast_min = g_fast_max * G_MIN_FACTOR;
            // Log-space para exploración uniforme en varias órdenes de magnitud
            log_g_fast = ParamRange(std::log(g_fast_min == 0.0 ? SMALL_LOG : g_fast_min),
                                    std::log(g_fast_max == 0.0 ? SMALL_LOG : g_fast_max));
        }

        if (use_i_slow)
        {
            // k1: tasa de apertura del canal lento, proporcional a fc (frecuencia de corte fast/slow)
            // Se usa fc en la escala del dt de la sinapsis (adimensional)
            const double fc_adim = fc * period_t / dt;
            const double k1_max = BOPublicConfig::K1_MAX_FACTOR * fc_adim;
            const double k1_min = BOPublicConfig::K1_MIN_FACTOR * fc_adim;

            // v_slow: umbral de la sigmoide lenta, todo el rango de v_pre + margen en ambos extremos (online)
            v_slow = ParamRange(v_pre_margin_min,
                                v_pre_margin_max);
            // s_slow: pendiente de la sigmoide lenta, más suave que s_fast (factores menores)
            const double s_slow_max = BOPublicConfig::S_SLOW_MAX_FACTOR / safe_v_pre_range;
            s_slow = ParamRange(BOPublicConfig::S_SLOW_MIN_FACTOR / safe_v_pre_range,
                                s_slow_max);

            // log_k1: en log-space porque abarca ~6 órdenes de magnitud
            log_k1 = ParamRange(std::log(k1_min == 0.0 ? SMALL_LOG : k1_min),
                                std::log(k1_max == 0.0 ? SMALL_LOG : k1_max));
            // log_R: R = k2/k1 en log-space; R grande -> m_slow pequeño, R pequeño -> m_slow grande
            log_R = ParamRange(std::log(R_MIN == 0.0 ? SMALL_LOG : R_MIN),
                               std::log(R_MAX == 0.0 ? SMALL_LOG : R_MAX));

            // k2_min para calcular m_max (peor caso: k2 mínimo -> m_slow máximo)
            const double k2_min = k1_min * R_MIN;
            // m_max: estado estacionario máximo de m_slow -> m_ss = (k1*sigma) / (k1*sigma + k2) con k1 max, sigma max, k2 min
            const double m_max_term = k1_max * chemical_sigmoid(s_slow_max, v_slow.min, v_pre_max);
            const double m_max = m_max_term / safe_divisor(m_max_term + k2_min);

            // g_slow_max: invierte g = I / (m_max * (v_post - e_syn)) en los dos extremos de driving force
            const double g_slow_max = std::max(std::abs(expected_i_margin_max / safe_divisor((v_post_max - e_syn_min) * m_max)),
                                               std::abs(expected_i_margin_min / safe_divisor((v_post_min - e_syn_max) * m_max)));
            // g_slow_min = g_slow_max * G_MIN_FACTOR
            const double g_slow_min = g_slow_max * G_MIN_FACTOR;
            // Log-space para exploración uniforme
            log_g_slow = ParamRange(std::log(g_slow_min == 0.0 ? SMALL_LOG : g_slow_min),
                                    std::log(g_slow_max == 0.0 ? SMALL_LOG : g_slow_max));
        }
    }

    BOParamRanges() = default;
};

// Buffers de padding preasignados para el filtro Butterworth bidireccional (uno por dirección)
struct EvaluationPadBuffers
{

    kfr::univector<double> padded_buff_12;
    kfr::univector<double> padded_buff_21;

    EvaluationPadBuffers(size_t padded_buff_size_12,
                         size_t padded_buff_size_21)
        : padded_buff_12(padded_buff_size_12),
          padded_buff_21(padded_buff_size_21)
    {
    }
};

// Functor de parada adicional para la BO online:
// comprueba el flag atómico stop_BO para permitir cancelación interactiva
struct StopFunctor
{

    static inline std::atomic<bool> *stop_BO_ptr = nullptr;

    template <typename BO, typename AggregatorFunction>
    bool operator()(const BO &bo, const AggregatorFunction &)
    {
        // Devuelve true si se ha solicitado parar -> Limbo detiene la optimización
        return stop_BO_ptr->load(std::memory_order_relaxed);
    }
};

// Parámetros de una sinapsis química (ambas componentes fast y slow)
struct ChemicalSynapseParams
{

    double e_syn;  // Potencial de reversa sináptico
    double g_fast; // Conductancia de la componente rápida
    double s_fast; // Pendiente de la sigmoide rápida
    double v_fast; // Umbral de voltaje de la sigmoide rápida
    double g_slow; // Conductancia de la componente lenta
    double k1;     // Tasa de apertura del canal lento
    double k2;     // Tasa de cierre del canal lento (k2 = k1 * R)
    double s_slow; // Pendiente de la sigmoide lenta
    double v_slow; // Umbral de voltaje de la sigmoide lenta
};

// Candidato bidireccional: parámetros para las dos sinapsis (1->2 y 2->1)
struct Candidate
{

    ChemicalSynapseParams params_12;
    ChemicalSynapseParams params_21;
};

// Resultado de evaluar un candidato: puntuaciones de rango y forma
struct ChemicalSynapseEvaluation
{

    double i_range_score;
    double i_shape_score;

    // Protección: valores no finitos se sustituyen por penalizaciones
    ChemicalSynapseEvaluation(double i_range_score_ = 0.0,
                              double i_shape_score_ = 0.0)
        : i_range_score(std::isfinite(i_range_score_) ? i_range_score_ : BOPublicConfig::VERY_BAD_RANGE_SCORE),
          i_shape_score(std::isfinite(i_shape_score_) ? i_shape_score_ : 0.0)
    {
    }
};

// Copia selectiva de parámetros sinápticos al array runtime, según las componentes activas.
// Se usa en el double-buffering para actualizar los parámetros del hilo RT de forma atómica.
inline void copy_selected_synapse_params(ChemicalSynapseParams &runtime_params,
                                         const ChemicalSynapseParams &params,
                                         unsigned int use_i_fast,
                                         unsigned int use_i_slow)
{

    if (use_i_fast || use_i_slow)
    {
        runtime_params.e_syn = params.e_syn;

        if (use_i_fast)
        {
            runtime_params.g_fast = params.g_fast;
            runtime_params.s_fast = params.s_fast;
            runtime_params.v_fast = params.v_fast;
        }

        if (use_i_slow)
        {
            runtime_params.g_slow = params.g_slow;
            runtime_params.v_slow = params.v_slow;
            runtime_params.k1 = params.k1;
            runtime_params.k2 = params.k2;
            runtime_params.s_slow = params.s_slow;
        }
    }
}

#endif
