#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <kfr/all.hpp>
#include "scaling.hpp"
#include <ChemicalSynapsis.h>

using namespace limbo;
using namespace nlohmann;

namespace BOPrivateConfig
{
    // Período mínimo de re-optimización de hiperparámetros del GP (en iteraciones)
    static constexpr int HP_PERIOD_MIN = 1;
    // Divisor: hp_period = max(HP_PERIOD_MIN, num_iters / HP_PERIOD_DIVISOR)
    static constexpr int HP_PERIOD_DIVISOR = 25;

    // Factores para escalar las iteraciones de los optimizadores internos con la dimensión
    static constexpr int RPROP_ITER_FACTOR = 50;  // iters_rprop = 50 * dim
    static constexpr int NLOPT_ITER_FACTOR = 150; // iters_nlopt = 150 * dim^2

    // Pesos relativos de la puntuación de rango y forma en el fitness total
    static constexpr double I_SHAPE_WEIGHT = 0.5;
    static constexpr double I_RANGE_WEIGH = 0.5;

    // --- Factores para calcular los rangos de búsqueda de los parámetros sinápticos ---

    // v_fast_min = v_pre_min + rango * V_FAST_MIN_FACTOR (la sigmoide rápida actúa en la parte alta)
    static constexpr double V_FAST_MIN_FACTOR = 0.25;
    // s_fast in [S_FAST_MIN_FACTOR, S_FAST_MAX_FACTOR] / rango_v_pre (inversamente proporcional)
    static constexpr double S_FAST_MIN_FACTOR = 4.12;
    static constexpr double S_FAST_MAX_FACTOR = 13.72;
    // s_slow in [S_SLOW_MIN_FACTOR, S_SLOW_MAX_FACTOR] / rango_v_pre
    static constexpr double S_SLOW_MIN_FACTOR = 1.72;
    static constexpr double S_SLOW_MAX_FACTOR = 3.43;
    // k1 in [K1_MIN_FACTOR * fc_adim, K1_MAX_FACTOR * fc_adim] (proporcional a la frecuencia de corte adimensional)
    static constexpr double K1_MAX_FACTOR = 29.845125;
    static constexpr double K1_MIN_FACTOR = 29.845125e-4;
    // e_syn se coloca a una distancia proporcional al rango_v_post:
    //   phase: [v_post_max + near*rango, v_post_max + far*rango] (por encima -> despolarizante)
    //   antiphase: [v_post_min - far*rango, v_post_min - near*rango] (por debajo -> hiperpolarizante)
    static constexpr double E_SYN_FAR_TERM = 3.86;
    static constexpr double E_SYN_NEAR_TERM = 0.2;
    // g_min = g_max * G_MIN_FACTOR (explora varias órdenes de magnitud de conductancia)
    static constexpr double G_MIN_FACTOR = 0.001;
    // R = k2/k1 in [R_MIN, R_MAX], en log-space
    static constexpr double R_MAX = 30.0;
    static constexpr double R_MIN = 0.001;
}

namespace BOPrivateConstants
{
    // Suma de pesos para normalizar el fitness
    static constexpr double TOTAL_WEIGHT = BOPrivateConfig::I_SHAPE_WEIGHT + BOPrivateConfig::I_RANGE_WEIGH;

    // Valor mínimo positivo (para log de valores cercanos a cero)
    static constexpr double SMALL_LOG = std::numeric_limits<double>::min();
}

// ==========================================
//  Configuración de Limbo (Bayesian Optimization)
// ==========================================
// Cada struct hereda de defaults y sobreescribe los parámetros relevantes.
// BO_PARAM = parámetro compilado en tiempo de compilación
// BO_DYN_PARAM = parámetro dinámico (se establece en runtime con set_xxx)
struct Params
{
    struct bayes_opt_bobase : public defaults::bayes_opt_bobase
    {
        // Acota las muestras al hipercubo [0,1]^dim
        BO_PARAM(bool, bounded, true);
        // Deshabilita la recolección de estadísticas (muestras, observaciones, etc.)
        BO_PARAM(bool, stats_enabled, false);
    };

    struct bayes_opt_boptimizer : public defaults::bayes_opt_boptimizer
    {
        // Cada cuántas iteraciones se re-optimizan los hiperparámetros del GP
        BO_DYN_PARAM(int, hp_period);
    };

    // Muestreo inicial Latin Hypercube Sampling
    struct init_lhs : public defaults::init_lhs
    {
        // Número de muestras iniciales del LHS
        BO_DYN_PARAM(int, samples);
    };

    // Configuración del kernel del GP
    struct kernel : public defaults::kernel
    {
        // Ruido observacional fijo (offline: señal determinista, poco ruido)
        BO_PARAM(double, noise, 1e-6);
        // No optimizar el ruido (es fijo)
        BO_PARAM(bool, optimize_noise, false);
    };

    // Kernel Squared Exponential con Automatic Relevance Determination (una longscale por dimensión)
    struct kernel_squared_exp_ard : public defaults::kernel_squared_exp_ard
    {
        // k=0: no se fuerza ningún número de componentes ARD
        BO_PARAM(int, k, 0);
        // Varianza de señal inicial del kernel
        BO_PARAM(double, sigma_sq, 1.0);
    };

    // Función de adquisición Expected Improvement
    struct acqui_ei : public defaults::acqui_ei
    {
        // Jitter (xi): controla la exploración vs explotación en EI
        BO_PARAM(double, jitter, 0.003);
    };

    // Optimizador sin gradiente NLOpt (LN_SBPLX = Subplex) para maximizar la adquisición
    struct opt_nloptnograd : public defaults::opt_nloptnograd
    {
        BO_DYN_PARAM(int, iterations);
        // -1 = deshabilitado (no parar por tolerancia)
        BO_PARAM(double, fun_tolerance, -1);
        BO_PARAM(double, xrel_tolerance, -1);
    };

    // Optimizador Rprop para los hiperparámetros del GP (gradiente)
    struct opt_rprop : public defaults::opt_rprop
    {
        BO_DYN_PARAM(int, iterations);
        BO_PARAM(double, eps_stop, 1e-6);
    };

    // Criterio de parada por número máximo de iteraciones
    struct stop_maxiterations : public defaults::stop_maxiterations
    {
        BO_DYN_PARAM(int, iterations);
    };

    // Media constante del GP (prior): fitness inicial esperado = 0.5
    struct mean_constant : public defaults::mean_constant
    {
        BO_PARAM(double, constant, 0.5);
    };
};

// Sigmoide sináptica: 1 / (1 + exp(s * (v_threshold - v_pre)))
static double chemical_sigmoid(double s,
                               double v_threshold,
                               double v_pre)
{

    return 1.0 / (1.0 + std::exp(s * (v_threshold - v_pre)));
}

// Rangos de búsqueda para cada parámetro sináptico.
// Limbo trabaja en [0,1]^dim; estos rangos mapean [0,1] al espacio real del parámetro.
// Algunos parámetros se optimizan en log-space (g_fast, g_slow, k1, R) para cubrir
// varias órdenes de magnitud.
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

    // Calcula los rangos reales de cada parámetro sináptico a partir de v_pre, v_post y corrientes esperadas
    void init(double v_pre_min,
              double v_pre_max,
              double v_post_min,
              double v_post_max,
              double expected_i_min,
              double expected_i_max,
              bool use_i_fast,
              bool use_i_slow,
              bool search_phase,
              double fc,
              double csv_step,
              double dt,
              size_t points_factor)
    {

        constexpr double G_MIN_FACTOR = BOPrivateConfig::G_MIN_FACTOR;
        constexpr double R_MIN = BOPrivateConfig::R_MIN;
        constexpr double SMALL_LOG = BOPrivateConstants::SMALL_LOG;
        constexpr double R_MAX = BOPrivateConfig::R_MAX;

        const double v_pre_range = v_pre_max - v_pre_min;
        const double v_post_range = v_post_max - v_post_min;

        // e_syn: excitadora -> por encima de v_post_max; inhibidora -> por debajo de v_post_min
        const double e_syn_far_final_term = v_post_range * BOPrivateConfig::E_SYN_FAR_TERM;
        const double e_syn_near_final_term = v_post_range * BOPrivateConfig::E_SYN_NEAR_TERM;
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

        // Margen sobre la corriente esperada (+-EXPECTED_I_MARGIN_FACTOR del rango) para ampliar la búsqueda de g
        const double expected_i_margin = (expected_i_max - expected_i_min) * BOPublicConfig::EXPECTED_I_MARGIN_FACTOR;
        const double expected_i_margin_min = expected_i_min - expected_i_margin;
        const double expected_i_margin_max = expected_i_max + expected_i_margin;
        const double safe_v_pre_range = safe_divisor(v_pre_range);

        if (use_i_fast)
        {
            // v_fast: umbral de la sigmoide rápida, restringido al cuarto superior de v_pre (solo responde a spikes)
            v_fast = ParamRange(v_pre_min + (v_pre_range * BOPrivateConfig::V_FAST_MIN_FACTOR),
                                v_pre_max);
            // s_fast: pendiente de la sigmoide rápida, inversamente proporcional al rango de v_pre (1/V)
            const double s_fast_max = BOPrivateConfig::S_FAST_MAX_FACTOR / safe_v_pre_range;
            s_fast = ParamRange(BOPrivateConfig::S_FAST_MIN_FACTOR / safe_v_pre_range,
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
            const double fc_adim = fc * csv_step / (dt * points_factor);
            const double k1_max = BOPrivateConfig::K1_MAX_FACTOR * fc_adim;
            const double k1_min = BOPrivateConfig::K1_MIN_FACTOR * fc_adim;

            // v_slow: umbral de la sigmoide lenta, cubre todo el rango de v_pre (puede activarse a cualquier nivel)
            v_slow = ParamRange(v_pre_min,
                                v_pre_max);
            // s_slow: pendiente de la sigmoide lenta, más suave que s_fast (factores menores)
            const double s_slow_max = BOPrivateConfig::S_SLOW_MAX_FACTOR / safe_v_pre_range;
            s_slow = ParamRange(BOPrivateConfig::S_SLOW_MIN_FACTOR / safe_v_pre_range,
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

// ==========================================
//  Tipos de Limbo: composición del pipeline de BO
// ==========================================

// Modelo GP: kernel SE-ARD con hiperparámetros optimizados vía Rprop,
// media constante como prior
using Model_t = model::GP<
    Params,
    kernel::SquaredExpARD<Params>,
    mean::Constant<Params>,
    model::gp::KernelLFOpt<Params, opt::Rprop<Params>>>;

// Función de adquisición: Expected Improvement
using Acqui_t = acqui::EI<Params, Model_t>;

// Optimizador de la adquisición: NLOpt Subplex (sin gradiente)
using AcquiOpt_t = opt::NLOptNoGrad<Params, nlopt::LN_SBPLX>;

// Estadísticas recogidas durante la optimización
using Stat_t = boost::fusion::vector<>;

// Inicialización: Latin Hypercube Sampling
using Init_t = init::LHS<Params>;

// Criterio de parada: número máximo de iteraciones
using Stop_t = stop::MaxIterations<Params>;

// BOptimizer configurado con todos los componentes
using BO_t = bayes_opt::BOptimizer<
    Params,
    stopcrit<Stop_t>,
    modelfun<Model_t>,
    acquifun<Acqui_t>,
    statsfun<Stat_t>,
    initfun<Init_t>,
    acquiopt<AcquiOpt_t>>;

// Decodifica un valor normalizado [0,1] al espacio real del parámetro usando su rango
static double decode_param(double x_val, const BOParamRanges::ParamRange &range)
{

    return range.min + (std::clamp(x_val, 0.0, 1.0) * range.range);
}

// Decodifica un vector de Limbo (en [0,1]^dim) a parámetros sinápticos reales.
// Los parámetros en log-space (g_fast, g_slow, k1, R) se des-logaritmizan con exp().
// k2 = k1 * R (no es un parámetro independiente).
static ChemicalSynapseParams decode_to_candidate(const Eigen::VectorXd &x,
                                                 const BOParamRanges &ranges,
                                                 bool use_i_fast,
                                                 bool use_i_slow)
{

    ChemicalSynapseParams candidate{};
    size_t idx = 0;

    candidate.e_syn = decode_param(x(idx++), ranges.e_syn);

    if (use_i_fast)
    {

        candidate.g_fast = std::exp(decode_param(x(idx++), ranges.log_g_fast));
        candidate.s_fast = decode_param(x(idx++), ranges.s_fast);
        candidate.v_fast = decode_param(x(idx++), ranges.v_fast);
    }

    if (use_i_slow)
    {
        candidate.g_slow = std::exp(decode_param(x(idx++), ranges.log_g_slow));
        candidate.v_slow = decode_param(x(idx++), ranges.v_slow);
        candidate.k1 = std::exp(decode_param(x(idx++), ranges.log_k1));

        // R = k2/k1 (ratio), se parametriza en log-space independientemente
        const double R = std::exp(decode_param(x(idx++), ranges.log_R));
        candidate.k2 = candidate.k1 * R;
        candidate.s_slow = decode_param(x(idx++), ranges.s_slow);
    }

    return candidate;
}

// Base para el functor de evaluación; proporciona dim_in como parámetro dinámico de Limbo
struct EvaluationFunctorBase
{

    BO_DYN_PARAM(int, dim_in);
};

// Functor que Limbo llama para evaluar cada candidato.
// operator() recibe un vector [0,1]^dim de Limbo, decodifica a parámetros sinápticos,
// simula la sinapsis, y devuelve la puntuación combinada de rango + forma.
template <typename Integrator, typename NeuronType,
          ResetStateFunc<NeuronType> ResetStateFuncType,
          GetVFunc<NeuronType> GetVFuncType>
struct EvaluationFunctor : public EvaluationFunctorBase
{
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    EvaluationFunctor(ChemicalSynapsisType &synapse_,
                      NeuronType &model_neur_,
                      const ScaledSigResult &scaled_result_,
                      EvaluationISigBuffers &buffers_,
                      const ConstantEvaluationVals &constant_evaluation_vals_,
                      const BOParamRanges &ranges_,
                      bool use_i_fast_,
                      bool use_i_slow_,
                      bool search_phase_,
                      size_t v_pre_sig_start_idx_,
                      double i_min_,
                      double i_max_,
                      bool verbose_,
                      ResetStateFuncType reset_state_neur_,
                      GetVFuncType get_v_neur_,
                      json *score_history_ptr_)
        : synapse(synapse_),
          model_neur(model_neur_),
          scaled_result(scaled_result_),
          buffers(buffers_),
          constant_evaluation_vals(constant_evaluation_vals_),
          ranges(ranges_),
          use_i_fast(use_i_fast_),
          use_i_slow(use_i_slow_),
          search_phase(search_phase_),
          v_pre_sig_start_idx(v_pre_sig_start_idx_),
          i_min(i_min_),
          i_max(i_max_),
          verbose(verbose_),
          reset_state_neur(reset_state_neur_),
          get_v_neur(get_v_neur_),
          score_history_ptr(score_history_ptr_)
    {
    }

    ChemicalSynapsisType &synapse;
    NeuronType &model_neur;
    const ScaledSigResult &scaled_result;
    EvaluationISigBuffers &buffers;
    const ConstantEvaluationVals &constant_evaluation_vals;
    const BOParamRanges &ranges;
    bool use_i_fast;
    bool use_i_slow;
    bool search_phase;
    size_t v_pre_sig_start_idx;
    double i_min;
    double i_max;
    bool verbose;
    ResetStateFuncType reset_state_neur;
    GetVFuncType get_v_neur;
    json *score_history_ptr;
    mutable size_t evaluation_count = 0;

    // La BO trabaja con dim_out = 1 (fitness escalar)
    BO_PARAM(int, dim_out, 1);

    Eigen::VectorXd operator()(const Eigen::VectorXd &x) const
    {

        ChemicalSynapseParams candidate = decode_to_candidate(x, ranges, use_i_fast, use_i_slow);

        ChemicalSynapseEvaluation evaluations = evaluate_candidate(
            candidate,
            synapse,
            model_neur,
            scaled_result,
            buffers,
            v_pre_sig_start_idx,
            constant_evaluation_vals,
            use_i_fast,
            use_i_slow,
            search_phase,
            i_min,
            i_max,
            reset_state_neur,
            get_v_neur);

        // Fitness = media ponderada de puntuación de rango y forma
        const double y = ((evaluations.i_range_score * BOPrivateConfig::I_RANGE_WEIGH) +
                          (evaluations.i_shape_score * BOPrivateConfig::I_SHAPE_WEIGHT)) /
                         BOPrivateConstants::TOTAL_WEIGHT;

        if (verbose)
        {
            std::cout << "Evaluation " << evaluation_count + 1 << ": " << y << " (range: " << evaluations.i_range_score << ", shape: " << evaluations.i_shape_score << ")" << std::endl;
        }

        if (score_history_ptr)
        {
            json &score_history = *score_history_ptr;
            score_history["scores"][evaluation_count] = y;
            score_history["range_scores"][evaluation_count] = evaluations.i_range_score;
            score_history["shape_scores"][evaluation_count] = evaluations.i_shape_score;
        }
        evaluation_count++;

        return limbo::tools::make_vector(y);
    }
};

// Declaraciones de parámetros dinámicos (requeridas por Limbo para la definición estática)
BO_DECLARE_DYN_PARAM(int, EvaluationFunctorBase, dim_in);
BO_DECLARE_DYN_PARAM(int, Params::bayes_opt_boptimizer, hp_period);
BO_DECLARE_DYN_PARAM(int, Params::init_lhs, samples);
BO_DECLARE_DYN_PARAM(int, Params::stop_maxiterations, iterations);
BO_DECLARE_DYN_PARAM(int, Params::opt_rprop, iterations);
BO_DECLARE_DYN_PARAM(int, Params::opt_nloptnograd, iterations);

// Cuenta el número de parámetros según las componentes activas:
// e_syn (1) + i_fast (3: g_fast, s_fast, v_fast) + i_slow (5: g_slow, v_slow, k1, R, s_slow)
static inline size_t count_params(bool use_i_fast,
                                  bool use_i_slow)
{

    size_t num_params = 1; // e_syn siempre
    if (use_i_fast)
        num_params += 3;
    if (use_i_slow)
        num_params += 5;
    return num_params;
}

// Función principal de BO offline
template <typename Integrator, typename NeuronType,
          CreateFunc<NeuronType> CreateFuncType,
          ResetStateFunc<NeuronType> ResetStateFuncType,
          GetVFunc<NeuronType> GetVFuncType>
std::optional<ChemicalSynapseParams> BO(const std::string &csv_path,
                                        size_t column_idx,
                                        double csv_step,
                                        double start_time,
                                        double stabilization_time,
                                        double evaluation_time,
                                        double observation_time,
                                        size_t initial_samples,
                                        size_t iterations,
                                        NumericIntegrator integrator,
                                        NeuronModel model,
                                        bool search_phase,
                                        bool check_drift,
                                        SynComponent syn_component,
                                        CreateFuncType create_neur,
                                        ResetStateFuncType reset_state_neur,
                                        GetVFuncType get_v_neur,
                                        typename NeuronType::variable neur_v_var,
                                        int syn_model_step_factor,
                                        double fc,
                                        double expected_i_min,
                                        double expected_i_max,
                                        double i_min,
                                        double i_max,
                                        bool verbose,
                                        const std::optional<std::string> &jsonl_history_file_path)
{
    const std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();

    if (csv_step <= 0.0 || evaluation_time <= 0.0 || observation_time <= 0.0 || stabilization_time < 0.0 || fc <= 0.0 || expected_i_min >= expected_i_max || i_min > i_max)
    {
        throw std::invalid_argument("Invalid arguments: csv_step, evaluation_time, observation_time, and fc must be positive; stabilization_time non-negative; expected_i_min must be less than expected_i_max; i_min must be less or equal to i_max");
    }

    // Determina qué componentes sinápticas optimizar
    const bool use_i_fast = (syn_component != SynComponent::ISLOW);
    const bool use_i_slow = (syn_component != SynComponent::IFAST);

    if (!use_i_fast && !use_i_slow)
    {
        throw std::invalid_argument("At least one of use_i_fast or use_i_slow must be true");
    }

    // Escala la señal del CSV al espacio del modelo neuronal
    const std::optional<ScaledSigResult> scaled_result_opt = scale_sig(
        csv_path, column_idx, csv_step, start_time, evaluation_time + stabilization_time,
        observation_time, integrator, model, check_drift);

    if (!scaled_result_opt)
    {
        return std::nullopt;
    }

    const ScaledSigResult &scaled_result = *scaled_result_opt;

    // Neurona post-sináptica modelo (recibe la corriente generada)
    NeuronType model_neur = create_neur(false);
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    // Inicializa la sinapsis, desactivando las componentes no usadas (a 0)
    typename ChemicalSynapsisType::ConstructorArgs syn_args{};
    if (use_i_slow && !use_i_fast)
    {

        syn_args.params[ChemicalSynapsisType::gfast] = 0.0;
        syn_args.params[ChemicalSynapsisType::sfast] = 0.0;
        syn_args.params[ChemicalSynapsisType::Vfast] = 0.0;
    }
    else if (use_i_fast && !use_i_slow)
    {

        syn_args.params[ChemicalSynapsisType::gslow] = 0.0;
        syn_args.params[ChemicalSynapsisType::Vslow] = 0.0;
        syn_args.params[ChemicalSynapsisType::k1] = 0.0;
        syn_args.params[ChemicalSynapsisType::k2] = 0.0;
        syn_args.params[ChemicalSynapsisType::sslow] = 0.0;
    }

    // Crea la sinapsis: neurona pre (vacía, su voltaje viene de la señal) y neurona post
    ChemicalSynapsisType synapse(create_neur(true), neur_v_var, model_neur, neur_v_var, syn_args, syn_model_step_factor);

    // Puntos de estabilización y evaluación
    const size_t stabilization_points = static_cast<size_t>(stabilization_time / csv_step);
    const size_t evaluation_points = scaled_result.sig.size() - stabilization_points;

    // Buffers preasignados para las corrientes (se reutilizan en cada evaluación)
    EvaluationISigBuffers buffers(evaluation_points, use_i_fast, use_i_slow);

    // Segmento de señal presináptica para la evaluación (sin estabilización)
    const kfr::univector_ref<const double> evaluation_v_pre_sig = scaled_result.sig.slice(stabilization_points, evaluation_points);

    // Rango dinámico de la neurona postsináptica (según el modelo)
    double v_post_min, v_post_max;
    if (model == NeuronModel::HINDMARSH_ROSE)
    {
        v_post_min = HindmarshRose::MIN;
        v_post_max = HindmarshRose::MAX;
    }
    else
    {
        throw std::invalid_argument("Unsupported neuron model");
    }

    // Min/max de la señal presináptica escalada (para calcular rangos de búsqueda)
    const double evaluation_v_pre_min = kfr::minof(evaluation_v_pre_sig);
    const double evaluation_v_pre_max = kfr::maxof(evaluation_v_pre_sig);

    // Calcula los rangos de búsqueda para cada parámetro sináptico
    BOParamRanges ranges;

    ranges.init(evaluation_v_pre_min, evaluation_v_pre_max, v_post_min, v_post_max,
                expected_i_min, expected_i_max,
                use_i_fast, use_i_slow,
                search_phase, fc,
                csv_step, scaled_result.dt, scaled_result.points_factor);

    // Precalcula las señales de referencia y rangos esperados de corriente
    const ConstantEvaluationVals constant_evaluation_vals = calc_constant_evaluation_vals(
        evaluation_v_pre_sig,
        evaluation_v_pre_min,
        evaluation_v_pre_max,
        csv_step,
        fc,
        use_i_fast,
        use_i_slow,
        expected_i_min,
        expected_i_max,
        search_phase);

    // --- Configuración de parámetros dinámicos de Limbo ---
    const int dim_in = static_cast<int>(count_params(use_i_fast, use_i_slow));
    const int iters = static_cast<int>(iterations);

    Params::init_lhs::set_samples(static_cast<int>(initial_samples));
    Params::stop_maxiterations::set_iterations(iters);
    // Re-optimizar hiperparámetros del GP cada iters/25 iteraciones (mínimo 1)
    Params::bayes_opt_boptimizer::set_hp_period(std::max(BOPrivateConfig::HP_PERIOD_MIN, iters / BOPrivateConfig::HP_PERIOD_DIVISOR));
    // Iteraciones de los optimizadores internos escalan con la dimensión
    Params::opt_rprop::set_iterations(BOPrivateConfig::RPROP_ITER_FACTOR * dim_in);
    Params::opt_nloptnograd::set_iterations(BOPrivateConfig::NLOPT_ITER_FACTOR * dim_in * dim_in);

    EvaluationFunctor<Integrator, NeuronType, ResetStateFuncType, GetVFuncType>::set_dim_in(dim_in);

    json history;
    json *score_history_ptr = nullptr;
    if (jsonl_history_file_path)
    {
        const size_t total_evaluations = initial_samples + iterations;
        json &score_history = history["score_history"];
        score_history["scores"] = std::vector<double>(total_evaluations);
        score_history["range_scores"] = std::vector<double>(total_evaluations);
        score_history["shape_scores"] = std::vector<double>(total_evaluations);
        score_history_ptr = &score_history;
    }

    EvaluationFunctor<Integrator, NeuronType, ResetStateFuncType, GetVFuncType> functor(
        synapse,
        model_neur,
        scaled_result,
        buffers,
        constant_evaluation_vals,
        ranges,
        use_i_fast,
        use_i_slow,
        search_phase,
        stabilization_points,
        i_min,
        i_max,
        verbose,
        reset_state_neur,
        get_v_neur,
        score_history_ptr);

    // Ejecuta la optimización bayesiana
    BO_t opt;

    opt.optimize(functor);

    const std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed = t_end - t_start;

    if (verbose)
    {
        std::cout << "Best: " << opt.best_observation()(0) << std::endl;
        std::cout << "Optimization time: " << elapsed.count() << " s" << std::endl;
    }

    // Decodifica la mejor muestra encontrada a parámetros sinápticos reales
    const ChemicalSynapseParams best_candidate = decode_to_candidate(opt.best_sample(), ranges, use_i_fast, use_i_slow);

    if (jsonl_history_file_path)
    {
        const Eigen::VectorXd &lengthscales = opt.model().kernel_function().ell();
        json &ls_json = history["ARD_lss"];
        json &best_params = history["best_params"];
        size_t idx = 0;

        history["optimization_time"] = elapsed.count();

        ls_json["e_syn"] = lengthscales(idx++);
        best_params["e_syn"] = best_candidate.e_syn;

        if (use_i_fast)
        {
            ls_json["g_fast"] = lengthscales(idx++);
            ls_json["s_fast"] = lengthscales(idx++);
            ls_json["v_fast"] = lengthscales(idx++);

            best_params["g_fast"] = best_candidate.g_fast;
            best_params["s_fast"] = best_candidate.s_fast;
            best_params["v_fast"] = best_candidate.v_fast;
        }
        if (use_i_slow)
        {
            ls_json["g_slow"] = lengthscales(idx++);
            ls_json["v_slow"] = lengthscales(idx++);
            ls_json["k1"] = lengthscales(idx++);
            ls_json["R"] = lengthscales(idx++);
            ls_json["s_slow"] = lengthscales(idx++);

            best_params["g_slow"] = best_candidate.g_slow;
            best_params["v_slow"] = best_candidate.v_slow;
            best_params["k1"] = best_candidate.k1;
            best_params["k2"] = best_candidate.k2;
            best_params["s_slow"] = best_candidate.s_slow;
        }

        std::ofstream fout(*jsonl_history_file_path, std::ios::app);
        if (fout.is_open())
        {
            fout << history.dump() << "\n";
        }
        else
        {
            std::cerr << "Warning: Could not open jsonl file " << *jsonl_history_file_path << " for appending.\n";
        }
    }

    return best_candidate;
}
