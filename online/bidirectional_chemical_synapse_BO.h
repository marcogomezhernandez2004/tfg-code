/*
 * Copyright (C) 2011 Georgia Institute of Technology, University of Utah,
 * Weill Cornell Medical College
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BIDIRECTIONAL_CHEMICAL_SYNAPSE_BO_H
#define BIDIRECTIONAL_CHEMICAL_SYNAPSE_BO_H

#include <default_gui_model.h>
#include <kfr/all.hpp>
#include <thread>
#include <atomic>
#include "utils.hpp"

// Módulo RTXI para sinapsis química bidireccional con Optimización Bayesiana.
// Implementa dos hilos:
//   - RT (execute): calcula corrientes sinápticas en tiempo real y recoge señales
//   - NRT (NRT_BO): ejecuta la BO en un hilo separado, sincronizándose con RT
//     mediante double-buffering atómico y flags RT_storing/stop_BO
class BidirectionalChemicalSynapseBO : public DefaultGUIModel
{

  Q_OBJECT

public:
  BidirectionalChemicalSynapseBO(void);
  virtual ~BidirectionalChemicalSynapseBO(void);

  // Función RT: se llama cada período de RTXI
  void execute(void);

  void createGUI(DefaultGUIModel::variable_t *, int);
  void customizeGUI(void);

  // EvaluationFunctor necesita acceso a los miembros para la evaluación
  friend struct EvaluationFunctor;

protected:
  virtual void update(DefaultGUIModel::update_flags_t);

private:
  double evaluations_completed;               // Estado: número de evaluaciones completadas
  double i_fast_12;                           // Estado: I_fast 1->2
  double i_slow_12;                           // Estado: I_slow 1->2
  double i_fast_21;                           // Estado: I_fast 2->1
  double i_slow_21;                           // Estado: I_slow 2->1
  unsigned int initial_samples;               // Muestras iniciales LHS de la BO
  unsigned int iterations;                    // Iteraciones de BO tras el muestreo inicial
  double evaluation_time, stabilization_time; // Tiempos en ms
  unsigned int search_phase;                  // 1 = fase (excitadora), 0 = antifase

  // Rangos esperados de corriente por dirección (objetivo de la BO)
  double expected_i_min_12, expected_i_max_12, expected_i_min_21, expected_i_max_21;

  // Frecuencias de corte Butterworth (kHz) por dirección
  double fc_1, fc_2;

  // Modos de adquisición de v_pre/v_post por dirección: 0 = manual, 1 = dinámico (desde inputs)
  unsigned int dynamic_v_pre_min_max_12, dynamic_v_post_min_max_12;
  unsigned int dynamic_v_pre_min_max_21, dynamic_v_post_min_max_21;
  // Rangos de voltaje por dirección (para calcular rangos de búsqueda)
  double v_pre_min_12, v_pre_max_12, v_post_min_12, v_post_max_12;
  double v_pre_min_21, v_pre_max_21, v_post_min_21, v_post_max_21;

  // Clamps de corriente de salida por dirección
  double i_min_12, i_max_12, i_min_21, i_max_21;

  // Verbose atómico (se puede cambiar desde GUI incluso durante la BO)
  std::atomic<unsigned int> verbose;

  double dt; // dt = integración paso del integrador RK6(5) en ms

  // Selección de componentes sinápticas por dirección
  unsigned int use_i_fast_12, use_i_slow_12, use_i_fast_21, use_i_slow_21;

  // Double-buffer de parámetros sinápticos: [0] y [1] para cada dirección.
  // El hilo NRT escribe en el alterno y publica con synapse_idx atómico.
  ChemicalSynapseParams params_21[2];
  ChemicalSynapseParams params_12[2];

  // Variables de estado del canal lento (m_slow) por dirección
  double m_slow_21, m_slow_12;
  double period; // Período RT en ms

  // Hilo NRT para la optimización bayesiana
  std::thread BO_NRT_thread;
  std::atomic<bool> stop_BO; // Flag de cancelación
  bool BO_running;           // true mientras la BO está activa

  // --- Variables de sincronización RT <-> NRT ---
  // synapse_idx: índice del buffer activo (0 o 1); escrito por NRT con release,
  //              leído por RT con acquire
  std::atomic<size_t> synapse_idx;
  // last_synapse_idx_read_RT: último índice que el RT ha leído;
  //              escrito por RT, leído por NRT para saber cuándo puede escribir el alterno
  std::atomic<size_t> last_synapse_idx_read_RT;

  // RT_storing: true cuando el NRT quiere que el RT recoja señales.
  //             NRT lo pone a true (release), RT lo pone a false cuando termina (release).
  std::atomic<bool> RT_storing;
  size_t storing_idx;  // Índice de escritura actual en los buffers de señal
  size_t num_elements; // Número total de muestras a recoger por evaluación

  // Buffers de señal compartidos: RT escribe, NRT lee después de RT_storing=false
  kfr::univector<double> i_fast_sig_12;
  kfr::univector<double> i_fast_sig_21;
  kfr::univector<double> i_slow_sig_12;
  kfr::univector<double> i_slow_sig_21;
  kfr::univector<double> v_pre_sig_12;
  kfr::univector<double> v_pre_sig_21;

  QPushButton *BO_button;

  void initParameters();

  // Integra m_slow un paso con Runge-Kutta 6(5)
  void runge_kutta_65(double (*f)(double, double, const ChemicalSynapseParams &), double &m_slow, double v_pre, double dt, const ChemicalSynapseParams &params);

  // Calcula las corrientes sinápticas
  double compute_i_slow(double &m_slow, double v_pre, double v_post, const ChemicalSynapseParams &params);
  double compute_i_fast(double v_pre, double v_post, const ChemicalSynapseParams &params);

  // Derivada de m_slow: dm/dt = k1*(1-m)*sigmoid - k2*m
  static double sm_chemical_synapse_m(double m_slow, double v_pre, const ChemicalSynapseParams &params);

  // Evalúa un candidato bidireccional (sincroniza con RT, recoge señales, calcula puntuaciones)
  ChemicalSynapseEvaluation evaluate_candidate(
      const Candidate &candidate,
      double fs,
      size_t effective_pad_12,
      size_t effective_pad_21,
      EvaluationPadBuffers &pad_buffers,
      double i_dist_max_12,
      double i_dist_max_21,
      size_t &curr_synapse_idx);

  // Decodifica un vector de Limbo a un candidato bidireccional
  Candidate decode_to_candidate(const Eigen::VectorXd &x,
                                const BOParamRanges &ranges_12,
                                const BOParamRanges &ranges_21);

  // Bucle principal de BO en el hilo NRT
  void NRT_BO(double period_t);

  void set_params_read_only(bool read_only);

  void init_syn_params_and_vars(ChemicalSynapseParams &params);
  // Espera activa hasta que RT lea el índice dado, o hasta cancelación
  bool wait_until_RT_read_idx_or_stop(size_t idx_to_achieve);

private slots:
  void toggle_BO_event(void);
  void stop_BO_event_async(void);
  void update_params_gui(void);
  void set_evaluations_completed(double evals);
};

#endif
