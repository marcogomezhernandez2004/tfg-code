#include <limbo/limbo.hpp>
#include "bidirectional_chemical_synapse_BO.h"
#include <nlohmann/json.hpp>
#include <chrono>
#include <fstream>
#include <cstdlib>

using namespace limbo;
using namespace nlohmann;

namespace BOPrivateConfig
{
    // Factor de longitud de padding para filtfilt
    static constexpr double PAD_LEN_FACTOR = 1.5;

    // Período mínimo de re-optimización de hiperparámetros del GP
    static constexpr int HP_PERIOD_MIN = 1;
    // hp_period = max(HP_PERIOD_MIN, iters / HP_PERIOD_DIVISOR)
    static constexpr int HP_PERIOD_DIVISOR = 25;

    // Factores para escalar las iteraciones de los optimizadores internos con la dimensión
    static constexpr int RPROP_ITER_FACTOR = 50;  // iters_rprop = 50 * dim
    static constexpr int NLOPT_ITER_FACTOR = 150; // iters_nlopt = 150 * dim^2

    // Pesos relativos de la puntuación de rango y forma en el fitness total
    static constexpr double I_SHAPE_WEIGHT = 0.5;
    static constexpr double I_RANGE_WEIGH = 0.5;
}

namespace BOPrivateConstants
{
    static constexpr double TOTAL_WEIGHT = BOPrivateConfig::I_SHAPE_WEIGHT + BOPrivateConfig::I_RANGE_WEIGH;
}

// ==========================================
//  Configuración de Limbo (Bayesian Optimization)
// ==========================================
struct Params
{
    struct bayes_opt_bobase : public defaults::bayes_opt_bobase
    {
        // Acota las muestras al hipercubo [0,1]^dim
        BO_PARAM(bool, bounded, true);
        // Deshabilita la recolección de estadísticas
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
        BO_DYN_PARAM(int, samples);
    };

    // Configuración del kernel del GP (online: más ruido porque la señal es real)
    struct kernel : public defaults::kernel
    {
        // Ruido observacional inicial más alto que offline (señal biológica ruidosa)
        BO_PARAM(double, noise, 0.05);
        // Permite optimizar el ruido (adapta a la variabilidad real de la señal)
        BO_PARAM(bool, optimize_noise, true);
    };

    // Kernel Squared Exponential ARD
    struct kernel_squared_exp_ard : public defaults::kernel_squared_exp_ard
    {
        BO_PARAM(int, k, 0);
        BO_PARAM(double, sigma_sq, 1.0);
    };

    // Función de adquisición Expected Improvement
    struct acqui_ei : public defaults::acqui_ei
    {
        // Jitter (xi) más bajo que offline: menos exploración, más explotación online
        BO_PARAM(double, jitter, 0.001);
    };

    // Optimizador sin gradiente NLOpt (Subplex) para maximizar la adquisición
    struct opt_nloptnograd : public defaults::opt_nloptnograd
    {
        BO_DYN_PARAM(int, iterations);
        BO_PARAM(double, fun_tolerance, -1);
        BO_PARAM(double, xrel_tolerance, -1);
    };

    // Optimizador Rprop para hiperparámetros del GP
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

    // Media constante del GP (prior): fitness esperado = 0.5
    struct mean_constant : public defaults::mean_constant
    {
        BO_PARAM(double, constant, 0.5);
    };
};

// Functor que Limbo llama para evaluar cada candidato.
// Decodifica el vector [0,1]^dim, configura las sinapsis vía double-buffering,
// espera a que el hilo RT recoja las señales, y calcula la puntuación.
struct EvaluationFunctor
{

    EvaluationFunctor(BidirectionalChemicalSynapseBO &module_,
                      const BOParamRanges &ranges_12_,
                      const BOParamRanges &ranges_21_,
                      double fs_,
                      size_t effective_pad_12_,
                      size_t effective_pad_21_,
                      EvaluationPadBuffers &pad_buffers_,
                      double i_dist_max_12_,
                      double i_dist_max_21_,
                      size_t curr_synapse_idx_,
                      json *score_history_ptr_)
        : module(module_),
          ranges_12(ranges_12_),
          ranges_21(ranges_21_),
          fs(fs_),
          effective_pad_12(effective_pad_12_),
          effective_pad_21(effective_pad_21_),
          pad_buffers(pad_buffers_),
          i_dist_max_12(i_dist_max_12_),
          i_dist_max_21(i_dist_max_21_),
          curr_synapse_idx(curr_synapse_idx_),
          score_history_ptr(score_history_ptr_)
    {
    }

    BidirectionalChemicalSynapseBO &module;
    const BOParamRanges &ranges_12;
    const BOParamRanges &ranges_21;
    double fs;
    size_t effective_pad_12, effective_pad_21;
    EvaluationPadBuffers &pad_buffers;
    double i_dist_max_12, i_dist_max_21;
    // Índice del buffer activo del double-buffering; se actualiza en cada evaluación
    mutable size_t curr_synapse_idx;
    json *score_history_ptr;
    mutable size_t evaluation_count = 0;

    // Dimensión de entrada (total de parámetros de ambas direcciones)
    BO_DYN_PARAM(int, dim_in);
    // La BO trabaja con fitness escalar
    BO_PARAM(int, dim_out, 1);

    Eigen::VectorXd operator()(const Eigen::VectorXd &x) const
    {
        // Decodifica el vector normalizado a parámetros sinápticos de ambas direcciones
        Candidate candidate = module.decode_to_candidate(x, ranges_12, ranges_21);

        // Evalúa el candidato: sincroniza con el hilo RT, recoge señales y calcula puntuaciones
        ChemicalSynapseEvaluation evaluations = module.evaluate_candidate(candidate,
                                                                          fs,
                                                                          effective_pad_12,
                                                                          effective_pad_21,
                                                                          pad_buffers,
                                                                          i_dist_max_12,
                                                                          i_dist_max_21,
                                                                          curr_synapse_idx);

        // Fitness = media ponderada de puntuación de rango y forma
        const double y = ((evaluations.i_range_score * BOPrivateConfig::I_RANGE_WEIGH) + (evaluations.i_shape_score * BOPrivateConfig::I_SHAPE_WEIGHT)) / BOPrivateConstants::TOTAL_WEIGHT;
        if (module.verbose.load(std::memory_order_relaxed))
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

// --- Tipos de Limbo ---

// Modelo GP: kernel SE-ARD, media constante, hiperparámetros optimizados vía Rprop
using Model_t = model::GP<
    Params,
    kernel::SquaredExpARD<Params>,
    mean::Constant<Params>,
    model::gp::KernelLFOpt<Params, opt::Rprop<Params>>>;

// Función de adquisición: Expected Improvement
using Acqui_t = acqui::EI<Params, Model_t>;

// Optimizador de la adquisición: NLOpt Subplex
using AcquiOpt_t = opt::NLOptNoGrad<Params, nlopt::LN_SBPLX>;

// Estadísticas recogidas durante la optimización
using Stat_t = boost::fusion::vector<>;

// Inicialización: Latin Hypercube Sampling
using Init_t = init::LHS<Params>;

// Criterio de parada: MaxIterations + StopFunctor (cancelación interactiva via stop_BO)
using Stop_t = boost::fusion::vector<
    stop::MaxIterations<Params>,
    StopFunctor>;

// BOptimizer con todos los componentes
using BO_t = bayes_opt::BOptimizer<
    Params,
    stopcrit<Stop_t>,
    modelfun<Model_t>,
    acquifun<Acqui_t>,
    statsfun<Stat_t>,
    initfun<Init_t>,
    acquiopt<AcquiOpt_t>>;

// Declaraciones de parámetros dinámicos (requeridas por Limbo)
BO_DECLARE_DYN_PARAM(int, EvaluationFunctor, dim_in);
BO_DECLARE_DYN_PARAM(int, Params::bayes_opt_boptimizer, hp_period);
BO_DECLARE_DYN_PARAM(int, Params::init_lhs, samples);
BO_DECLARE_DYN_PARAM(int, Params::stop_maxiterations, iterations);
BO_DECLARE_DYN_PARAM(int, Params::opt_rprop, iterations);
BO_DECLARE_DYN_PARAM(int, Params::opt_nloptnograd, iterations);

// Decodifica un valor normalizado [0,1] al espacio real del parámetro
static double decode_param(double x_val, const BOParamRanges::ParamRange &range)
{

    return range.min + (std::clamp(x_val, 0.0, 1.0) * range.range);
}

// Decodifica parámetros de una dirección sináptica (1->2 o 2->1).
// Los parámetros en log-space se des-logaritmizan con exp().
// idx se pasa por referencia para avanzar el puntero en el vector de Limbo.
static void decode_to_params(const Eigen::VectorXd &x,
                             size_t &idx,
                             ChemicalSynapseParams &params,
                             unsigned int use_i_fast,
                             unsigned int use_i_slow,
                             const BOParamRanges &ranges)
{

    if (use_i_fast || use_i_slow)
    {
        params.e_syn = decode_param(x(idx++), ranges.e_syn);

        if (use_i_fast)
        {

            params.g_fast = std::exp(decode_param(x(idx++), ranges.log_g_fast));
            params.s_fast = decode_param(x(idx++), ranges.s_fast);
            params.v_fast = decode_param(x(idx++), ranges.v_fast);
        }

        if (use_i_slow)
        {
            params.g_slow = std::exp(decode_param(x(idx++), ranges.log_g_slow));
            params.v_slow = decode_param(x(idx++), ranges.v_slow);
            params.k1 = std::exp(decode_param(x(idx++), ranges.log_k1));

            // R = k2/k1 (ratio), en log-space
            const double R = std::exp(decode_param(x(idx++), ranges.log_R));
            params.k2 = params.k1 * R;
            params.s_slow = decode_param(x(idx++), ranges.s_slow);
        }
    }
}

// Decodifica el vector completo de Limbo a un candidato bidireccional
Candidate BidirectionalChemicalSynapseBO::decode_to_candidate(const Eigen::VectorXd &x,
                                                              const BOParamRanges &ranges_12,
                                                              const BOParamRanges &ranges_21)
{
    Candidate candidate{};
    size_t idx = 0;
    // Primero los parámetros de la dirección 1->2, luego 2->1 (orden fijo en el vector)
    decode_to_params(x, idx, candidate.params_12,
                     use_i_fast_12, use_i_slow_12, ranges_12);
    decode_to_params(x, idx, candidate.params_21,
                     use_i_fast_21, use_i_slow_21, ranges_21);
    return candidate;
}

// Cuenta los parámetros de una dirección:
// e_syn (1) + i_fast (3: g, s, v) + i_slow (5: g, v, k1, R, s)
static size_t count_params_one_direction(unsigned int use_i_fast,
                                         unsigned int use_i_slow)
{

    size_t num_params = 0;
    if (use_i_fast || use_i_slow)
    {
        num_params += 1; // e_syn
        if (use_i_fast)
            num_params += 3;
        if (use_i_slow)
            num_params += 5;
    }
    return num_params;
}

// ==========================================
//  NRT_BO: bucle principal de BO en el hilo No-Real-Time
// ==========================================
// Este método se ejecuta en un thread separado (NRT).
// Configura la BO, lanza la optimización, y al finalizar publica los mejores
// parámetros al hilo RT mediante double-buffering atómico.
void BidirectionalChemicalSynapseBO::NRT_BO(double period_t)
{
    const std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();

    const double fs = 1.0 / safe_divisor(period_t);           // fs en kHz (period en ms)
    num_elements = static_cast<size_t>(evaluation_time * fs); // Puntos a almacenar por evaluación
    if (num_elements < 1)
    {
        QMetaObject::invokeMethod(this, "stop_BO_event_async", Qt::QueuedConnection);
        return;
    }
    // Factor de padding para filtfilt, proporcional a fs
    const size_t pad_len_factor_fs = static_cast<size_t>(BOPrivateConfig::PAD_LEN_FACTOR * fs);

    const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
    const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;

    if (!use_syn_12 && !use_syn_21)
    {
        QMetaObject::invokeMethod(this, "stop_BO_event_async", Qt::QueuedConnection);
        return;
    }

    // Calcula el padding del filtro para cada dirección (proporcional a fs/fc)
    size_t effective_pad_12 = 0;
    size_t padded_buff_size_12 = 0;
    if (use_syn_12)
    {

        effective_pad_12 = std::min(num_elements - 1,
                                    static_cast<size_t>(pad_len_factor_fs / safe_divisor(fc_1)));
        padded_buff_size_12 = num_elements + (2 * effective_pad_12);
    }

    size_t effective_pad_21 = 0;
    size_t padded_buff_size_21 = 0;
    if (use_syn_21)
    {
        effective_pad_21 = std::min(num_elements - 1,
                                    static_cast<size_t>(pad_len_factor_fs / safe_divisor(fc_2)));
        padded_buff_size_21 = num_elements + (2 * effective_pad_21);
    }

    // Buffers de padding preasignados para el filtro Butterworth de cada dirección
    EvaluationPadBuffers pad_buffers(padded_buff_size_12, padded_buff_size_21);

    // Inicializa los rangos de búsqueda para cada dirección sináptica.
    // Online usa v_pre/v_post por dirección (dinámicos o fijos según configuración)
    BOParamRanges ranges_12;
    if (use_syn_12)
    {
        // v_pre/v_post configurados explícitamente para la dirección 1->2
        ranges_12.init(v_pre_min_12, v_pre_max_12, v_post_min_12, v_post_max_12,
                       expected_i_min_12, expected_i_max_12,
                       use_i_fast_12, use_i_slow_12,
                       search_phase, fc_1, period_t, dt);

        // Redimensiona los buffers de señal del hilo RT para la fase de recogida
        v_pre_sig_12.resize(num_elements);
        if (use_i_fast_12)
            i_fast_sig_12.resize(num_elements);
        if (use_i_slow_12)
            i_slow_sig_12.resize(num_elements);
    }

    BOParamRanges ranges_21;
    if (use_syn_21)
    {
        // v_pre/v_post configurados explícitamente para la dirección 2->1
        ranges_21.init(v_pre_min_21, v_pre_max_21, v_post_min_21, v_post_max_21,
                       expected_i_min_21, expected_i_max_21,
                       use_i_fast_21, use_i_slow_21,
                       search_phase, fc_2, period_t, dt);
        v_pre_sig_21.resize(num_elements);
        if (use_i_fast_21)
            i_fast_sig_21.resize(num_elements);
        if (use_i_slow_21)
            i_slow_sig_21.resize(num_elements);
    }

    // Dimensión total = params(1->2) + params(2->1)
    const int dim_in = static_cast<int>(count_params_one_direction(use_i_fast_12, use_i_slow_12) + count_params_one_direction(use_i_fast_21, use_i_slow_21));
    const int iters = static_cast<int>(iterations);

    // --- Configuración de parámetros dinámicos de Limbo ---
    Params::init_lhs::set_samples(static_cast<int>(initial_samples));
    Params::stop_maxiterations::set_iterations(iters);
    Params::bayes_opt_boptimizer::set_hp_period(std::max(BOPrivateConfig::HP_PERIOD_MIN, iters / BOPrivateConfig::HP_PERIOD_DIVISOR));
    Params::opt_rprop::set_iterations(BOPrivateConfig::RPROP_ITER_FACTOR * dim_in);
    Params::opt_nloptnograd::set_iterations(BOPrivateConfig::NLOPT_ITER_FACTOR * dim_in * dim_in);

    EvaluationFunctor::set_dim_in(dim_in);

    // Precalcula i_dist_max solo cuando se usa una sola componente por dirección
    double i_dist_max_12 = 0.0;
    if (use_syn_12)
    {

        if (use_i_fast_12 != use_i_slow_12)
        {
            i_dist_max_12 = calculate_expected_i_dist_max(expected_i_min_12, expected_i_max_12);
        }
    }
    double i_dist_max_21 = 0.0;
    if (use_syn_21)
    {
        if (use_i_fast_21 != use_i_slow_21)
        {
            i_dist_max_21 = calculate_expected_i_dist_max(expected_i_min_21, expected_i_max_21);
        }
    }

    const char *jsonl_path_env = std::getenv("JSONL_HISTORY_FILE_PATH");
    const std::optional<std::string> jsonl_history_file_path = (jsonl_path_env && jsonl_path_env[0] != '\0') ? std::optional<std::string>(jsonl_path_env) : std::nullopt;

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

    EvaluationFunctor functor(*this,
                              ranges_12,
                              ranges_21,
                              fs,
                              effective_pad_12,
                              effective_pad_21,
                              pad_buffers,
                              i_dist_max_12,
                              i_dist_max_21,
                              // Índice actual del double-buffer
                              synapse_idx.load(std::memory_order_relaxed),
                              score_history_ptr);

    // Notifica a la GUI que empiezan las evaluaciones (reset contador)
    QMetaObject::invokeMethod(this, "set_evaluations_completed", Qt::QueuedConnection,
                              Q_ARG(double, 0));

    BO_t opt;
    try
    {
        // Ejecuta la optimización bayesiana
        opt.optimize(functor);
    }
    catch (const StopEvaluation &)
    {
    }

    const std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed = t_end - t_start;

    // Comprobación clave: asegurar que existan individuos evaluados
    const bool has_evaluations = !opt.observations().empty();

    if (verbose.load(std::memory_order_relaxed))
    {
        if (has_evaluations)
        {
            std::cout << "Best: " << opt.best_observation()(0) << std::endl;
        }
        std::cout << "Optimization time: " << elapsed.count() << " s" << std::endl;
    }

    Candidate best_candidate{};
    if (has_evaluations)
    {
        // Decodifica los mejores parámetros encontrados
        best_candidate = decode_to_candidate(opt.best_sample(), ranges_12, ranges_21);
    }

    if (jsonl_history_file_path)
    {
        history["optimization_time"] = elapsed.count();

        // Solo accedemos al modelo y guardamos resultados si hubo evaluaciones
        if (has_evaluations)
        {
            const Eigen::VectorXd &lengthscales = opt.model().kernel_function().ell();
            json &ls_json = history["ARD_lss"];
            json &best_params = history["best_params"];
            size_t ls_idx = 0;

            if (use_syn_12)
            {
                json &ls_12 = ls_json["1->2"];
                json &bp_12 = best_params["1->2"];
                const ChemicalSynapseParams &bc_params_12 = best_candidate.params_12;

                ls_12["e_syn"] = lengthscales(ls_idx++);
                bp_12["e_syn"] = bc_params_12.e_syn;

                if (use_i_fast_12)
                {
                    ls_12["g_fast"] = lengthscales(ls_idx++);
                    ls_12["s_fast"] = lengthscales(ls_idx++);
                    ls_12["v_fast"] = lengthscales(ls_idx++);

                    bp_12["g_fast"] = bc_params_12.g_fast;
                    bp_12["s_fast"] = bc_params_12.s_fast;
                    bp_12["v_fast"] = bc_params_12.v_fast;
                }
                if (use_i_slow_12)
                {
                    ls_12["g_slow"] = lengthscales(ls_idx++);
                    ls_12["v_slow"] = lengthscales(ls_idx++);
                    ls_12["k1"] = lengthscales(ls_idx++);
                    ls_12["R"] = lengthscales(ls_idx++);
                    ls_12["s_slow"] = lengthscales(ls_idx++);

                    bp_12["g_slow"] = bc_params_12.g_slow;
                    bp_12["v_slow"] = bc_params_12.v_slow;
                    bp_12["k1"] = bc_params_12.k1;
                    bp_12["k2"] = bc_params_12.k2;
                    bp_12["s_slow"] = bc_params_12.s_slow;
                }
            }

            if (use_syn_21)
            {
                json &ls_21 = ls_json["2->1"];
                json &bp_21 = best_params["2->1"];
                const ChemicalSynapseParams &bc_params_21 = best_candidate.params_21;

                ls_21["e_syn"] = lengthscales(ls_idx++);
                bp_21["e_syn"] = bc_params_21.e_syn;

                if (use_i_fast_21)
                {
                    ls_21["g_fast"] = lengthscales(ls_idx++);
                    ls_21["s_fast"] = lengthscales(ls_idx++);
                    ls_21["v_fast"] = lengthscales(ls_idx++);

                    bp_21["g_fast"] = bc_params_21.g_fast;
                    bp_21["s_fast"] = bc_params_21.s_fast;
                    bp_21["v_fast"] = bc_params_21.v_fast;
                }
                if (use_i_slow_21)
                {
                    ls_21["g_slow"] = lengthscales(ls_idx++);
                    ls_21["v_slow"] = lengthscales(ls_idx++);
                    ls_21["k1"] = lengthscales(ls_idx++);
                    ls_21["R"] = lengthscales(ls_idx++);
                    ls_21["s_slow"] = lengthscales(ls_idx++);

                    bp_21["g_slow"] = bc_params_21.g_slow;
                    bp_21["v_slow"] = bc_params_21.v_slow;
                    bp_21["k1"] = bc_params_21.k1;
                    bp_21["k2"] = bc_params_21.k2;
                    bp_21["s_slow"] = bc_params_21.s_slow;
                }
            }
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

    // Si se canceló durante la optimización o si no hay parámetros que actualizar
    if (stop_BO.load(std::memory_order_relaxed) || !has_evaluations)
    {
        QMetaObject::invokeMethod(this, "stop_BO_event_async", Qt::QueuedConnection);
        return;
    }

    // --- Publicación final de los mejores parámetros al hilo RT ---
    // Espera a que el hilo RT haya leído el índice actual antes de escribir el nuevo
    const size_t curr_synapse_idx = functor.curr_synapse_idx;

    if (!wait_until_RT_read_idx_or_stop(curr_synapse_idx))
    {
        QMetaObject::invokeMethod(this, "stop_BO_event_async", Qt::QueuedConnection);
        return;
    }

    // Escribe en el buffer alterno y publica con release
    const size_t new_synapse_idx = 1 - curr_synapse_idx;

    copy_selected_synapse_params(params_12[new_synapse_idx],
                                 best_candidate.params_12,
                                 use_i_fast_12,
                                 use_i_slow_12);
    copy_selected_synapse_params(params_21[new_synapse_idx],
                                 best_candidate.params_21,
                                 use_i_fast_21,
                                 use_i_slow_21);
    // Release: hace visibles las escrituras anteriores al hilo RT
    synapse_idx.store(new_synapse_idx, std::memory_order_release);

    QMetaObject::invokeMethod(this, "stop_BO_event_async", Qt::QueuedConnection);
}