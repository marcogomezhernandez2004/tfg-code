#ifndef BO_H
#define BO_H

#include <limbo/limbo.hpp>
#include <Eigen/Core>
#include <cstddef>
#include <optional>
#include <string>
#include "utils.hpp"
#include "evaluation.hpp"
#include <nlohmann/json.hpp>
#include <fstream>

// Función principal de Optimización Bayesiana offline.
// Plantilla parametrizada por integrador numérico, tipo de neurona y funciones auxiliares.
// Lee una señal de CSV, la escala, y busca los parámetros sinápticos que maximicen
// la puntuación combinada de rango + forma de las corrientes generadas.
// Devuelve nullopt si no se puede escalar la señal (frecuencia de burst incompatible).
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
                                        typename NeuronType::VarType neur_v_var,
                                        int syn_model_step_factor,
                                        double fc,
                                        double expected_i_min,
                                        double expected_i_max,
                                        double i_min,
                                        double i_max,
                                        bool verbose,
                                        const std::optional<std::string> &jsonl_history_file_path);

#include "BO.tpp"

#endif
