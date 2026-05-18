#ifndef EVALUATION_H
#define EVALUATION_H

#include <cmath>
#include <cstddef>
#include <kfr/all.hpp>
#include <ChemicalSynapsis.h>
#include "scaling.hpp"
#include "utils.hpp"

namespace EvaluationPublicConfig
{
    // Puntuación de rango muy mala; se asigna si el resultado no es finito
    inline constexpr double VERY_BAD_RANGE_SCORE = -1e6;
}

// Resultado de evaluar un candidato de sinapsis: puntuaciones de rango y forma
struct ChemicalSynapseEvaluation
{
    // Puntuación de cuánto se acerca el rango [min,max] de la corriente al esperado
    double i_range_score;

    // Puntuación de la correlación de Pearson entre la forma de onda de la corriente y la referencia
    double i_shape_score;

    // Protección: si los valores no son finitos, se sustituyen por penalizaciones
    ChemicalSynapseEvaluation(double i_range_score_ = 0.0,
                              double i_shape_score_ = 0.0)
        : i_range_score(std::isfinite(i_range_score_) ? i_range_score_ : EvaluationPublicConfig::VERY_BAD_RANGE_SCORE),
          i_shape_score(std::isfinite(i_shape_score_) ? i_shape_score_ : 0.0)
    {
    }
};

// Valores constantes de evaluación precalculados una sola vez antes de la BO:
// señales de referencia (centradas), factores de Pearson, y rangos esperados de corriente
struct ConstantEvaluationVals
{
    // Señal de referencia i_fast centrada (media 0) y su factor sqrt(sum(x^2)) para Pearson
    kfr::univector<double> ref_i_fast_sig_centered;
    double ref_i_fast_sig_factor;
    // Señal de referencia i_slow centrada y su factor
    kfr::univector<double> ref_i_slow_sig_centered;
    double ref_i_slow_sig_factor;

    // Rangos esperados de corriente para el cálculo de range_score
    double ref_i_fast_min;
    double ref_i_fast_max;
    double i_fast_dist_max; // Distancia máxima admisible (para normalizar el error de rango)

    double ref_i_slow_min;
    double ref_i_slow_max;
    double i_slow_dist_max;
};

// Buffers preasignados para almacenar las señales de corriente simuladas de cada candidato
// Se reutilizan en cada evaluación para evitar realocaciones
struct EvaluationISigBuffers
{

    EvaluationISigBuffers(size_t size_to_reserve,
                          bool use_i_fast,
                          bool use_i_slow)
    {
        if (use_i_fast)
            i_fast_sig.resize(size_to_reserve);
        if (use_i_slow)
            i_slow_sig.resize(size_to_reserve);
    }

    kfr::univector<double> i_fast_sig;
    kfr::univector<double> i_slow_sig;
};

// Precalcula los valores constantes de evaluación (señales de referencia y rangos)
// a partir de la señal presináptica filtrada
ConstantEvaluationVals calc_constant_evaluation_vals(
    const kfr::univector_ref<double> &v_pre_sig,
    double v_pre_min,
    double v_pre_max,
    double csv_step,
    double fc,
    bool use_i_fast,
    bool use_i_slow,
    double expected_i_min,
    double expected_i_max,
    bool search_phase);

// Evalúa un candidato sináptico: simula la sinapsis con los parámetros dados,
// recoge las corrientes i_fast/i_slow, y calcula las puntuaciones de rango y forma
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
    GetVFuncType get_v_neur);

#include "evaluation.tpp"

#endif
