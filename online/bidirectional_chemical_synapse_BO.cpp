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

#include "bidirectional_chemical_synapse_BO.h"
#include <chrono>
#include <algorithm>
#include <main_window.h>

namespace ModulePrivateConfig
{
  // Frecuencia de corte por defecto (kHz) para separar i_fast de i_slow
  static constexpr double FILTER_FC = 0.3;

  // Margen mínimo permitido para m_slow (evita bloqueo numérico fuera de [0,1])
  static constexpr double M_SLOW_MARGIN = 1e-6;
}

namespace ModuleConstants
{
  // Rango con margen para clampear m_slow (variable de gating lento)
  static constexpr double M_SLOW_MIN = -ModulePrivateConfig::M_SLOW_MARGIN;
  static constexpr double M_SLOW_MAX = 1.0 + ModulePrivateConfig::M_SLOW_MARGIN;
}

extern "C" Plugin::Object *
createRTXIPlugin(void)
{

  return new BidirectionalChemicalSynapseBO();
}

// Definición de todos los parámetros, estados, entradas y salidas del módulo RTXI
static DefaultGUIModel::variable_t vars[] = {
    {"BO evaluations completed", "Finishes when this is initial samples + iterations", DefaultGUIModel::STATE},
    {"I_fast 1->2 (nA)", "", DefaultGUIModel::STATE},
    {"I_slow 1->2 (nA)", "", DefaultGUIModel::STATE},
    {"I_fast 2->1 (nA)", "", DefaultGUIModel::STATE},
    {"I_slow 2->1 (nA)", "", DefaultGUIModel::STATE},

    {"BO initial samples", "Number of initialization samples for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"BO iterations", "Number of BO iterations after initial sampling", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"BO evaluation time (ms)", "Time to record signals per evaluation", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO stabilization time (ms)", "Wait time after setting params before recording", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO search phase (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"BO current min to achieve 1->2 (nA)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO current max to achieve 1->2 (nA)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO current min to achieve 2->1 (nA)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO current max to achieve 2->1 (nA)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO cutoff frequency 1 (kHz)", "To separate the I_fast and I_slow for BO in synapse 1->2", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO cutoff frequency 2 (kHz)", "To separate the I_fast and I_slow for BO in synapse 2->1", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Dynamic V_pre min and max 1->2 (1/0)", "1 = Enable, 0 = Disable; necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"V_pre min 1->2 (mV)", "Necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_pre max 1->2 (mV)", "Necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Dynamic V_post min and max 1->2 (1/0)", "1 = Enable, 0 = Disable; necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"V_post min 1->2 (mV)", "Necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_post max 1->2 (mV)", "Necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Dynamic V_pre min and max 2->1 (1/0)", "1 = Enable, 0 = Disable; necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"V_pre min 2->1 (mV)", "Necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_pre max 2->1 (mV)", "Necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Dynamic V_post min and max 2->1 (1/0)", "1 = Enable, 0 = Disable; necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"V_post min 2->1 (mV)", "Necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_post max 2->1 (mV)", "Necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Current min 1->2 (nA)", "Fixed output clamp min for current 1->2", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Current max 1->2 (nA)", "Fixed output clamp max for current 1->2", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Current min 2->1 (nA)", "Fixed output clamp min for current 2->1", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Current max 2->1 (nA)", "Fixed output clamp max for current 2->1", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Verbose (1/0)", "Enable/disable BO candidate evaluation logging", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},

    {"dt (ms)", "Integration step dt in ms", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Use I_fast 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Use I_slow 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Use I_fast 2->1 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Use I_slow 2->1 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},

    {"E_syn 1->2 (mV)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"g_fast 1->2 (uS)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"s_fast 1->2 (1/mV)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_fast 1->2 (mV)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"g_slow 1->2 (uS)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"k1 1->2 (1/ms)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"k2 1->2 (1/ms)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"s_slow 1->2 (1/mV)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_slow 1->2 (mV)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"E_syn 2->1 (mV)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"g_fast 2->1 (uS)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"s_fast 2->1 (1/mV)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_fast 2->1 (mV)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"g_slow 2->1 (uS)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"k1 2->1 (1/ms)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"k2 2->1 (1/ms)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"s_slow 2->1 (1/mV)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_slow 2->1 (mV)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Current 1->2 (nA)", "Total synaptic current 1->2", DefaultGUIModel::OUTPUT},
    {"Current 2->1 (nA)", "Total synaptic current 2->1", DefaultGUIModel::OUTPUT},

    {"V_pre 1->2 (mV)", "Presynaptic membrane potential 1->2", DefaultGUIModel::INPUT},
    {"V_post 1->2 (mV)", "Postsynaptic membrane potential 1->2", DefaultGUIModel::INPUT},
    {"V_pre 2->1 (mV)", "Presynaptic membrane potential 2->1", DefaultGUIModel::INPUT},
    {"V_post 2->1 (mV)", "Postsynaptic membrane potential 2->1", DefaultGUIModel::INPUT},
    {"V_pre min 1->2 (mV)", "Dynamic V_pre min 1->2", DefaultGUIModel::INPUT},
    {"V_pre max 1->2 (mV)", "Dynamic V_pre max 1->2", DefaultGUIModel::INPUT},
    {"V_post min 1->2 (mV)", "Dynamic V_post min 1->2", DefaultGUIModel::INPUT},
    {"V_post max 1->2 (mV)", "Dynamic V_post max 1->2", DefaultGUIModel::INPUT},
    {"V_pre min 2->1 (mV)", "Dynamic V_pre min 2->1", DefaultGUIModel::INPUT},
    {"V_pre max 2->1 (mV)", "Dynamic V_pre max 2->1", DefaultGUIModel::INPUT},
    {"V_post min 2->1 (mV)", "Dynamic V_post min 2->1", DefaultGUIModel::INPUT},
    {"V_post max 2->1 (mV)", "Dynamic V_post max 2->1", DefaultGUIModel::INPUT},
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

BidirectionalChemicalSynapseBO::BidirectionalChemicalSynapseBO(void)
    : DefaultGUIModel("Bidirectional Chemical Synapse BO", ::vars, ::num_vars)
{

  setWhatsThis("<p><b>Bidirectional Chemical Synapse BO</b></p>");
  DefaultGUIModel::createGUI(vars, num_vars);
  initParameters();
  customizeGUI();
  update(INIT);
  refresh();
  QTimer::singleShot(0, this, SLOT(resizeMe()));
}

BidirectionalChemicalSynapseBO::~BidirectionalChemicalSynapseBO(void)
{
  // Si la BO está corriendo, señaliza parada y espera a que termine el hilo NRT
  if (BO_NRT_thread.joinable())
  {

    stop_BO.store(true, std::memory_order_relaxed);
    BO_NRT_thread.join();
  }
}

// Integrador Runge-Kutta de orden 5 con 6 etapas para la variable m_slow del canal lento.
// Usa la misma función de derivada f(m_slow, v_pre, params) en cada etapa.
void BidirectionalChemicalSynapseBO::runge_kutta_65(double (*f)(double, double, const ChemicalSynapseParams &), double &m_slow, double v_pre, double dt, const ChemicalSynapseParams &params)
{

  double apoyo, retorno;
  double k[6];

  // 6 etapas del integrador RK5 con coeficientes de Butcher específicos
  retorno = (*f)(m_slow, v_pre, params);
  k[0] = dt * retorno;
  apoyo = m_slow + k[0] * 0.2;

  retorno = (*f)(apoyo, v_pre, params);
  k[1] = dt * retorno;
  apoyo = m_slow + k[0] * 0.075 + k[1] * 0.225;

  retorno = (*f)(apoyo, v_pre, params);
  k[2] = dt * retorno;
  apoyo = m_slow + k[0] * 0.3 - k[1] * 0.9 + k[2] * 1.2;

  retorno = (*f)(apoyo, v_pre, params);
  k[3] = dt * retorno;
  apoyo = m_slow + k[0] * 0.075 + k[1] * 0.675 - k[2] * 0.6 + k[3] * 0.75;

  retorno = (*f)(apoyo, v_pre, params);
  k[4] = dt * retorno;
  apoyo = m_slow + k[0] * 0.660493827160493 + k[1] * 2.5 - k[2] * 5.185185185185185 + k[3] * 3.888888888888889 - k[4] * 0.864197530864197;

  retorno = (*f)(apoyo, v_pre, params);
  k[5] = dt * retorno;

  // Actualiza m_slow con la combinación ponderada de las etapas
  m_slow += k[0] * 0.098765432098765 +
            k[2] * 0.396825396825396 +
            k[3] * 0.231481481481481 +
            k[4] * 0.308641975308641 -
            k[5] * 0.035714285714285;
}

// Derivada de m_slow: dm/dt = k1*(1-m)*sigmoid(v_pre) - k2*m
double BidirectionalChemicalSynapseBO::sm_chemical_synapse_m(double m_slow, double v_pre, const ChemicalSynapseParams &params)
{

  return (params.k1 * (1.0 - m_slow) * chemical_sigmoid(params.s_slow, params.v_slow, v_pre)) -
         (params.k2 * m_slow);
}

// Calcula i_slow: integra m_slow un paso con RK5 de 6 etapas y devuelve g_slow * m * (v_post - e_syn)
double BidirectionalChemicalSynapseBO::compute_i_slow(double &m_slow, double v_pre, double v_post, const ChemicalSynapseParams &params)
{
  runge_kutta_65(sm_chemical_synapse_m, m_slow, v_pre, dt, params);
  // Clampea m_slow para evitar que se salga de [0,1] por errores numéricos
  m_slow = std::clamp(m_slow, ModuleConstants::M_SLOW_MIN, ModuleConstants::M_SLOW_MAX);
  return params.g_slow * m_slow * (v_post - params.e_syn);
}

// Calcula i_fast: corriente instantánea sin dinámica temporal
// i_fast = g_fast * sigmoid(v_pre) * (v_post - e_syn)
double BidirectionalChemicalSynapseBO::compute_i_fast(double v_pre, double v_post, const ChemicalSynapseParams &params)
{

  return params.g_fast * (v_post - params.e_syn) * chemical_sigmoid(params.s_fast, params.v_fast, v_pre);
}

// ==========================================
//  execute(): función de tiempo real (hilo RT)
// ==========================================
// Se llama cada período de RTXI. Calcula las corrientes sinápticas y,
// si la BO está recogiendo datos (RT_storing), almacena las señales.
void BidirectionalChemicalSynapseBO::execute(void)
{

  // Solo actualiza v_pre/v_post dinámicos cuando la BO no está corriendo
  // (para evitar cambios de rango durante la optimización)
  if (!BO_running)
  {
    if (dynamic_v_pre_min_max_12)
    {
      v_pre_min_12 = input(4);
      v_pre_max_12 = input(5);
    }

    if (dynamic_v_post_min_max_12)
    {
      v_post_min_12 = input(6);
      v_post_max_12 = input(7);
    }

    if (dynamic_v_pre_min_max_21)
    {
      v_pre_min_21 = input(8);
      v_pre_max_21 = input(9);
    }

    if (dynamic_v_post_min_max_21)
    {
      v_post_min_21 = input(10);
      v_post_max_21 = input(11);
    }
  }

  const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
  const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;

  double v_pre_12 = 0.0;
  double v_post_12 = 0.0;
  double v_pre_21 = 0.0;
  double v_post_21 = 0.0;

  // --- Lectura del double-buffer atómico ---
  // Comprueba si el hilo NRT está pidiendo datos (RT_storing)
  const bool aux_RT_storing = RT_storing.load(std::memory_order_acquire);

  // Lee el índice del buffer activo (publicado por el hilo NRT con release)
  const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_acquire);

  // Confirma al hilo NRT que ya leyó este índice (para que pueda escribir el alterno)
  last_synapse_idx_read_RT.store(curr_synapse_idx, std::memory_order_relaxed);

  // --- Cálculo de corrientes sinápticas dirección 1->2 ---
  if (use_syn_12)
  {
    v_pre_12 = input(0);
    v_post_12 = input(1);

    const ChemicalSynapseParams &curr_params_12 = params_12[curr_synapse_idx];
    if (use_i_slow_12)
      i_slow_12 = compute_i_slow(m_slow_12, v_pre_12, v_post_12, curr_params_12);
    if (use_i_fast_12)
      i_fast_12 = compute_i_fast(v_pre_12, v_post_12, curr_params_12);
  }

  const double val_i_12 = i_fast_12 + i_slow_12;
  output(0) = std::clamp(val_i_12, i_min_12, i_max_12);

  // --- Cálculo de corrientes sinápticas dirección 2->1 ---
  if (use_syn_21)
  {
    v_pre_21 = input(2);
    v_post_21 = input(3);

    const ChemicalSynapseParams &curr_params_21 = params_21[curr_synapse_idx];
    if (use_i_slow_21)
      i_slow_21 = compute_i_slow(m_slow_21, v_pre_21, v_post_21, curr_params_21);
    if (use_i_fast_21)
      i_fast_21 = compute_i_fast(v_pre_21, v_post_21, curr_params_21);
  }

  const double val_i_21 = i_fast_21 + i_slow_21;
  output(1) = std::clamp(val_i_21, i_min_21, i_max_21);

  // --- Almacenamiento de señales para la BO (cuando RT_storing está activo) ---
  // El hilo NRT activa RT_storing = true y el RT llena los buffers muestra a muestra
  if (aux_RT_storing)
  {

    if (storing_idx < num_elements)
    {
      if (use_syn_12)
      {
        v_pre_sig_12[storing_idx] = v_pre_12;
        if (use_i_fast_12)
          i_fast_sig_12[storing_idx] = i_fast_12;
        if (use_i_slow_12)
          i_slow_sig_12[storing_idx] = i_slow_12;
      }
      if (use_syn_21)
      {
        v_pre_sig_21[storing_idx] = v_pre_21;
        if (use_i_fast_21)
          i_fast_sig_21[storing_idx] = i_fast_21;
        if (use_i_slow_21)
          i_slow_sig_21[storing_idx] = i_slow_21;
      }
      storing_idx++;
    }
    else
    {
      // Buffer lleno: señaliza al hilo NRT que los datos están listos
      RT_storing.store(false, std::memory_order_release);
    }
  }
}

void BidirectionalChemicalSynapseBO::initParameters(void)
{

  evaluations_completed = 0.0;
  initial_samples = 40u;
  iterations = 200u;
  evaluation_time = 2000.0;
  stabilization_time = 1000.0;
  search_phase = 1u;
  expected_i_min_12 = 0.0;
  expected_i_max_12 = 0.0;
  expected_i_min_21 = 0.0;
  expected_i_max_21 = 0.0;
  constexpr double FILTER_FC = ModulePrivateConfig::FILTER_FC;
  fc_1 = FILTER_FC;
  fc_2 = FILTER_FC;

  dynamic_v_pre_min_max_12 = 0u;
  v_pre_min_12 = 0.0;
  v_pre_max_12 = 0.0;

  dynamic_v_post_min_max_12 = 0u;
  v_post_min_12 = 0.0;
  v_post_max_12 = 0.0;

  dynamic_v_pre_min_max_21 = 0u;
  v_pre_min_21 = 0.0;
  v_pre_max_21 = 0.0;

  dynamic_v_post_min_max_21 = 0u;
  v_post_min_21 = 0.0;
  v_post_max_21 = 0.0;

  i_min_12 = 0.0;
  i_max_12 = 0.0;
  i_min_21 = 0.0;
  i_max_21 = 0.0;

  i_fast_12 = 0.0;
  i_slow_12 = 0.0;
  i_fast_21 = 0.0;
  i_slow_21 = 0.0;

  verbose.store(0u, std::memory_order_relaxed);

  dt = 0.001;

  use_i_fast_12 = 1u;
  use_i_slow_12 = 1u;
  use_i_fast_21 = 1u;
  use_i_slow_21 = 1u;

  init_syn_params_and_vars(params_12[0]);
  init_syn_params_and_vars(params_12[1]);
  init_syn_params_and_vars(params_21[0]);
  init_syn_params_and_vars(params_21[1]);

  m_slow_12 = 0.0;
  m_slow_21 = 0.0;

  // Inicializa el double-buffer con índice 0
  synapse_idx.store(0, std::memory_order_relaxed);
  last_synapse_idx_read_RT.store(0, std::memory_order_relaxed);

  stop_BO.store(false, std::memory_order_relaxed);
  RT_storing.store(false, std::memory_order_relaxed);
  BO_running = false;

  // Conecta el puntero estático del StopFunctor al flag atómico de esta instancia
  StopFunctor::stop_BO_ptr = &stop_BO;
}

void BidirectionalChemicalSynapseBO::init_syn_params_and_vars(ChemicalSynapseParams &params)
{
  params.e_syn = 0.0;
  params.g_fast = 0.0;
  params.s_fast = 0.0;
  params.v_fast = 0.0;
  params.g_slow = 0.0;
  params.k1 = 0.0;
  params.k2 = 0.0;
  params.s_slow = 0.0;
  params.v_slow = 0.0;
}

// Espera activa (polling) hasta que el hilo RT haya leído el índice dado,
// o hasta que se solicite parar. Devuelve false si hubo cancelación.
bool BidirectionalChemicalSynapseBO::wait_until_RT_read_idx_or_stop(size_t idx_to_achieve)
{
  const std::chrono::duration<double, std::milli> active_wait_duration(BOPublicConfig::ACTIVE_WAIT_MS);

  while (last_synapse_idx_read_RT.load(std::memory_order_relaxed) != idx_to_achieve)
  {
    if (stop_BO.load(std::memory_order_relaxed))
      return false;

    std::this_thread::sleep_for(active_wait_duration);
  }

  return true;
}

// Gestiona los eventos de la GUI de RTXI (init, modify, period, pause, unpause)
void BidirectionalChemicalSynapseBO::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag)
  {
  case INIT:
  {

    period = RT::System::getInstance()->getPeriod() * 1e-6; // ns -> ms

    setState("BO evaluations completed", evaluations_completed);
    setState("I_fast 1->2 (nA)", i_fast_12);
    setState("I_slow 1->2 (nA)", i_slow_12);
    setState("I_fast 2->1 (nA)", i_fast_21);
    setState("I_slow 2->1 (nA)", i_slow_21);

    setParameter("BO initial samples", initial_samples);
    setParameter("BO iterations", iterations);
    setParameter("BO evaluation time (ms)", evaluation_time);
    setParameter("BO stabilization time (ms)", stabilization_time);
    setParameter("BO search phase (1/0)", search_phase);
    setParameter("BO current min to achieve 1->2 (nA)", expected_i_min_12);
    setParameter("BO current max to achieve 1->2 (nA)", expected_i_max_12);
    setParameter("BO current min to achieve 2->1 (nA)", expected_i_min_21);
    setParameter("BO current max to achieve 2->1 (nA)", expected_i_max_21);
    setParameter("BO cutoff frequency 1 (kHz)", fc_1);
    setParameter("BO cutoff frequency 2 (kHz)", fc_2);

    setParameter("Dynamic V_pre min and max 1->2 (1/0)", dynamic_v_pre_min_max_12);
    setParameter("V_pre min 1->2 (mV)", v_pre_min_12);
    setParameter("V_pre max 1->2 (mV)", v_pre_max_12);
    setParameter("Dynamic V_post min and max 1->2 (1/0)", dynamic_v_post_min_max_12);
    setParameter("V_post min 1->2 (mV)", v_post_min_12);
    setParameter("V_post max 1->2 (mV)", v_post_max_12);
    setParameter("Dynamic V_pre min and max 2->1 (1/0)", dynamic_v_pre_min_max_21);
    setParameter("V_pre min 2->1 (mV)", v_pre_min_21);
    setParameter("V_pre max 2->1 (mV)", v_pre_max_21);
    setParameter("Dynamic V_post min and max 2->1 (1/0)", dynamic_v_post_min_max_21);
    setParameter("V_post min 2->1 (mV)", v_post_min_21);
    setParameter("V_post max 2->1 (mV)", v_post_max_21);

    setParameter("Current min 1->2 (nA)", i_min_12);
    setParameter("Current max 1->2 (nA)", i_max_12);
    setParameter("Current min 2->1 (nA)", i_min_21);
    setParameter("Current max 2->1 (nA)", i_max_21);

    setParameter("Verbose (1/0)", verbose.load(std::memory_order_relaxed));

    setParameter("dt (ms)", dt);

    setParameter("Use I_fast 1->2 (1/0)", use_i_fast_12);
    setParameter("Use I_slow 1->2 (1/0)", use_i_slow_12);
    setParameter("Use I_fast 2->1 (1/0)", use_i_fast_21);
    setParameter("Use I_slow 2->1 (1/0)", use_i_slow_21);

    update_params_gui();

    break;
  }
  case MODIFY:
  {
    // Los clamps y verbose se pueden cambiar durante la BO
    i_min_12 = getParameter("Current min 1->2 (nA)").toDouble();
    i_max_12 = getParameter("Current max 1->2 (nA)").toDouble();
    i_min_21 = getParameter("Current min 2->1 (nA)").toDouble();
    i_max_21 = getParameter("Current max 2->1 (nA)").toDouble();

    verbose.store(getParameter("Verbose (1/0)").toUInt(), std::memory_order_relaxed);

    // El resto de parámetros solo se puede cambiar cuando la BO no está corriendo
    if (!BO_running)
    {
      initial_samples = getParameter("BO initial samples").toUInt();
      iterations = getParameter("BO iterations").toUInt();
      evaluation_time = getParameter("BO evaluation time (ms)").toDouble();
      stabilization_time = getParameter("BO stabilization time (ms)").toDouble();
      search_phase = getParameter("BO search phase (1/0)").toUInt();
      expected_i_min_12 = getParameter("BO current min to achieve 1->2 (nA)").toDouble();
      expected_i_max_12 = getParameter("BO current max to achieve 1->2 (nA)").toDouble();
      expected_i_min_21 = getParameter("BO current min to achieve 2->1 (nA)").toDouble();
      expected_i_max_21 = getParameter("BO current max to achieve 2->1 (nA)").toDouble();
      fc_1 = getParameter("BO cutoff frequency 1 (kHz)").toDouble();
      fc_2 = getParameter("BO cutoff frequency 2 (kHz)").toDouble();

      dynamic_v_pre_min_max_12 = getParameter("Dynamic V_pre min and max 1->2 (1/0)").toUInt();
      const double new_v_pre_min_12 = getParameter("V_pre min 1->2 (mV)").toDouble();
      const double new_v_pre_max_12 = getParameter("V_pre max 1->2 (mV)").toDouble();
      if (!dynamic_v_pre_min_max_12)
      {
        v_pre_min_12 = new_v_pre_min_12;
        v_pre_max_12 = new_v_pre_max_12;
      }

      dynamic_v_post_min_max_12 = getParameter("Dynamic V_post min and max 1->2 (1/0)").toUInt();
      const double new_v_post_min_12 = getParameter("V_post min 1->2 (mV)").toDouble();
      const double new_v_post_max_12 = getParameter("V_post max 1->2 (mV)").toDouble();
      if (!dynamic_v_post_min_max_12)
      {
        v_post_min_12 = new_v_post_min_12;
        v_post_max_12 = new_v_post_max_12;
      }

      dynamic_v_pre_min_max_21 = getParameter("Dynamic V_pre min and max 2->1 (1/0)").toUInt();
      const double new_v_pre_min_21 = getParameter("V_pre min 2->1 (mV)").toDouble();
      const double new_v_pre_max_21 = getParameter("V_pre max 2->1 (mV)").toDouble();
      if (!dynamic_v_pre_min_max_21)
      {
        v_pre_min_21 = new_v_pre_min_21;
        v_pre_max_21 = new_v_pre_max_21;
      }

      dynamic_v_post_min_max_21 = getParameter("Dynamic V_post min and max 2->1 (1/0)").toUInt();
      const double new_v_post_min_21 = getParameter("V_post min 2->1 (mV)").toDouble();
      const double new_v_post_max_21 = getParameter("V_post max 2->1 (mV)").toDouble();
      if (!dynamic_v_post_min_max_21)
      {
        v_post_min_21 = new_v_post_min_21;
        v_post_max_21 = new_v_post_max_21;
      }

      dt = getParameter("dt (ms)").toDouble();

      use_i_fast_12 = getParameter("Use I_fast 1->2 (1/0)").toUInt();
      use_i_slow_12 = getParameter("Use I_slow 1->2 (1/0)").toUInt();
      use_i_fast_21 = getParameter("Use I_fast 2->1 (1/0)").toUInt();
      use_i_slow_21 = getParameter("Use I_slow 2->1 (1/0)").toUInt();
      if (!use_i_fast_12)
        i_fast_12 = 0.0;
      if (!use_i_slow_12)
        i_slow_12 = 0.0;
      if (!use_i_fast_21)
        i_fast_21 = 0.0;
      if (!use_i_slow_21)
        i_slow_21 = 0.0;

      // Actualiza los parámetros sinápticos del buffer activo actual
      const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_relaxed);

      params_12[curr_synapse_idx].e_syn = getParameter("E_syn 1->2 (mV)").toDouble();
      params_12[curr_synapse_idx].g_fast = getParameter("g_fast 1->2 (uS)").toDouble();
      params_12[curr_synapse_idx].s_fast = getParameter("s_fast 1->2 (1/mV)").toDouble();
      params_12[curr_synapse_idx].v_fast = getParameter("V_fast 1->2 (mV)").toDouble();
      params_12[curr_synapse_idx].g_slow = getParameter("g_slow 1->2 (uS)").toDouble();
      params_12[curr_synapse_idx].k1 = getParameter("k1 1->2 (1/ms)").toDouble();
      params_12[curr_synapse_idx].k2 = getParameter("k2 1->2 (1/ms)").toDouble();
      params_12[curr_synapse_idx].s_slow = getParameter("s_slow 1->2 (1/mV)").toDouble();
      params_12[curr_synapse_idx].v_slow = getParameter("V_slow 1->2 (mV)").toDouble();

      params_21[curr_synapse_idx].e_syn = getParameter("E_syn 2->1 (mV)").toDouble();
      params_21[curr_synapse_idx].g_fast = getParameter("g_fast 2->1 (uS)").toDouble();
      params_21[curr_synapse_idx].s_fast = getParameter("s_fast 2->1 (1/mV)").toDouble();
      params_21[curr_synapse_idx].v_fast = getParameter("V_fast 2->1 (mV)").toDouble();
      params_21[curr_synapse_idx].g_slow = getParameter("g_slow 2->1 (uS)").toDouble();
      params_21[curr_synapse_idx].k1 = getParameter("k1 2->1 (1/ms)").toDouble();
      params_21[curr_synapse_idx].k2 = getParameter("k2 2->1 (1/ms)").toDouble();
      params_21[curr_synapse_idx].s_slow = getParameter("s_slow 2->1 (1/mV)").toDouble();
      params_21[curr_synapse_idx].v_slow = getParameter("V_slow 2->1 (mV)").toDouble();
    }

    break;
  }

  case UNPAUSE:
    break;

  case PERIOD:
  {
    // Si cambia el período RT, se para la BO
    // (porque cambiaría el número de puntos a almacenar)
    const double new_period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
    if (new_period != period)
    {
      period = new_period;
      if (BO_running)
      {

        stop_BO.store(true, std::memory_order_relaxed); // Porque cambiaría el número de puntos a almacenar
      }
    }

    break;
  }

  case PAUSE:
    output(0) = 0;
    output(1) = 0;
    break;

  default:
    break;
  }
}

void BidirectionalChemicalSynapseBO::customizeGUI(void)
{
  QGridLayout *customlayout = DefaultGUIModel::getLayout();

  BO_button = new QPushButton("Start BO");
  QObject::connect(BO_button, SIGNAL(clicked()), this, SLOT(toggle_BO_event()));
  customlayout->addWidget(BO_button, 0, 0);

  parameter["I_fast 1->2 (nA)"].label->hide();
  parameter["I_fast 1->2 (nA)"].edit->hide();
  parameter["I_slow 1->2 (nA)"].label->hide();
  parameter["I_slow 1->2 (nA)"].edit->hide();
  parameter["I_fast 2->1 (nA)"].label->hide();
  parameter["I_fast 2->1 (nA)"].edit->hide();
  parameter["I_slow 2->1 (nA)"].label->hide();
  parameter["I_slow 2->1 (nA)"].edit->hide();

  setLayout(customlayout);
}

// Alterna entre iniciar y detener la BO desde la GUI.
// Lanza el hilo NRT para la optimización bayesiana.
void BidirectionalChemicalSynapseBO::toggle_BO_event(void)
{
  if (!BO_running)
  {
    // Si hay un hilo previo (ya terminado), hacer join antes de lanzar uno nuevo
    if (BO_NRT_thread.joinable())
    {
      BO_NRT_thread.join();
      stop_BO.store(false, std::memory_order_relaxed);
    }
    BO_running = true;
    BO_button->setText("Stop BO");
    set_params_read_only(true);

    // Lanza la BO en un hilo NRT separado, pasando el período actual
    BO_NRT_thread = std::thread(&BidirectionalChemicalSynapseBO::NRT_BO, this, period);
  }
  else
  {
    // Solicita parada (el hilo NRT la detectará en StopFunctor o evaluate_candidate)
    stop_BO.store(true, std::memory_order_relaxed);
  }
}

// Callback invocado desde el hilo NRT vía QMetaObject::invokeMethod cuando la BO termina
void BidirectionalChemicalSynapseBO::stop_BO_event_async(void)
{
  update_params_gui();
  set_params_read_only(false);
  BO_button->setText("Start BO");
  BO_running = false;
}

// Bloquea/desbloquea los campos de la GUI según si la BO está corriendo
void BidirectionalChemicalSynapseBO::set_params_read_only(bool read_only)
{

  parameter["BO initial samples"].edit->setReadOnly(read_only);
  parameter["BO iterations"].edit->setReadOnly(read_only);
  parameter["BO evaluation time (ms)"].edit->setReadOnly(read_only);
  parameter["BO stabilization time (ms)"].edit->setReadOnly(read_only);
  parameter["BO search phase (1/0)"].edit->setReadOnly(read_only);
  parameter["BO current min to achieve 1->2 (nA)"].edit->setReadOnly(read_only);
  parameter["BO current max to achieve 1->2 (nA)"].edit->setReadOnly(read_only);
  parameter["BO current min to achieve 2->1 (nA)"].edit->setReadOnly(read_only);
  parameter["BO current max to achieve 2->1 (nA)"].edit->setReadOnly(read_only);
  parameter["BO cutoff frequency 1 (kHz)"].edit->setReadOnly(read_only);
  parameter["BO cutoff frequency 2 (kHz)"].edit->setReadOnly(read_only);

  parameter["Dynamic V_pre min and max 1->2 (1/0)"].edit->setReadOnly(read_only);
  parameter["V_pre min 1->2 (mV)"].edit->setReadOnly(read_only);
  parameter["V_pre max 1->2 (mV)"].edit->setReadOnly(read_only);
  parameter["Dynamic V_post min and max 1->2 (1/0)"].edit->setReadOnly(read_only);
  parameter["V_post min 1->2 (mV)"].edit->setReadOnly(read_only);
  parameter["V_post max 1->2 (mV)"].edit->setReadOnly(read_only);
  parameter["Dynamic V_pre min and max 2->1 (1/0)"].edit->setReadOnly(read_only);
  parameter["V_pre min 2->1 (mV)"].edit->setReadOnly(read_only);
  parameter["V_pre max 2->1 (mV)"].edit->setReadOnly(read_only);
  parameter["Dynamic V_post min and max 2->1 (1/0)"].edit->setReadOnly(read_only);
  parameter["V_post min 2->1 (mV)"].edit->setReadOnly(read_only);
  parameter["V_post max 2->1 (mV)"].edit->setReadOnly(read_only);

  parameter["dt (ms)"].edit->setReadOnly(read_only);

  parameter["Use I_fast 1->2 (1/0)"].edit->setReadOnly(read_only);
  parameter["Use I_slow 1->2 (1/0)"].edit->setReadOnly(read_only);
  parameter["Use I_fast 2->1 (1/0)"].edit->setReadOnly(read_only);
  parameter["Use I_slow 2->1 (1/0)"].edit->setReadOnly(read_only);

  if (use_i_fast_12 || use_i_slow_12)
  {
    parameter["E_syn 1->2 (mV)"].edit->setReadOnly(read_only);
    if (use_i_fast_12)
    {
      parameter["g_fast 1->2 (uS)"].edit->setReadOnly(read_only);
      parameter["s_fast 1->2 (1/mV)"].edit->setReadOnly(read_only);
      parameter["V_fast 1->2 (mV)"].edit->setReadOnly(read_only);
    }
    if (use_i_slow_12)
    {
      parameter["g_slow 1->2 (uS)"].edit->setReadOnly(read_only);
      parameter["k1 1->2 (1/ms)"].edit->setReadOnly(read_only);
      parameter["k2 1->2 (1/ms)"].edit->setReadOnly(read_only);
      parameter["s_slow 1->2 (1/mV)"].edit->setReadOnly(read_only);
      parameter["V_slow 1->2 (mV)"].edit->setReadOnly(read_only);
    }
  }

  if (use_i_fast_21 || use_i_slow_21)
  {
    parameter["E_syn 2->1 (mV)"].edit->setReadOnly(read_only);
    if (use_i_fast_21)
    {
      parameter["g_fast 2->1 (uS)"].edit->setReadOnly(read_only);
      parameter["s_fast 2->1 (1/mV)"].edit->setReadOnly(read_only);
      parameter["V_fast 2->1 (mV)"].edit->setReadOnly(read_only);
    }
    if (use_i_slow_21)
    {
      parameter["g_slow 2->1 (uS)"].edit->setReadOnly(read_only);
      parameter["k1 2->1 (1/ms)"].edit->setReadOnly(read_only);
      parameter["k2 2->1 (1/ms)"].edit->setReadOnly(read_only);
      parameter["s_slow 2->1 (1/mV)"].edit->setReadOnly(read_only);
      parameter["V_slow 2->1 (mV)"].edit->setReadOnly(read_only);
    }
  }
}

// Actualiza los campos de la GUI con los parámetros sinápticos actuales del double-buffer
void BidirectionalChemicalSynapseBO::update_params_gui(void)
{

  const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_acquire);

  // Confirma lectura del índice para liberar el buffer alterno
  last_synapse_idx_read_RT.store(curr_synapse_idx, std::memory_order_relaxed);

  if (use_i_fast_12 || use_i_slow_12)
  {
    setParameter("E_syn 1->2 (mV)", params_12[curr_synapse_idx].e_syn);
    if (use_i_fast_12)
    {
      setParameter("g_fast 1->2 (uS)", params_12[curr_synapse_idx].g_fast);
      setParameter("s_fast 1->2 (1/mV)", params_12[curr_synapse_idx].s_fast);
      setParameter("V_fast 1->2 (mV)", params_12[curr_synapse_idx].v_fast);
    }
    if (use_i_slow_12)
    {
      setParameter("g_slow 1->2 (uS)", params_12[curr_synapse_idx].g_slow);
      setParameter("k1 1->2 (1/ms)", params_12[curr_synapse_idx].k1);
      setParameter("k2 1->2 (1/ms)", params_12[curr_synapse_idx].k2);
      setParameter("s_slow 1->2 (1/mV)", params_12[curr_synapse_idx].s_slow);
      setParameter("V_slow 1->2 (mV)", params_12[curr_synapse_idx].v_slow);
    }
  }

  if (use_i_fast_21 || use_i_slow_21)
  {
    setParameter("E_syn 2->1 (mV)", params_21[curr_synapse_idx].e_syn);
    if (use_i_fast_21)
    {
      setParameter("g_fast 2->1 (uS)", params_21[curr_synapse_idx].g_fast);
      setParameter("s_fast 2->1 (1/mV)", params_21[curr_synapse_idx].s_fast);
      setParameter("V_fast 2->1 (mV)", params_21[curr_synapse_idx].v_fast);
    }
    if (use_i_slow_21)
    {
      setParameter("g_slow 2->1 (uS)", params_21[curr_synapse_idx].g_slow);
      setParameter("k1 2->1 (1/ms)", params_21[curr_synapse_idx].k1);
      setParameter("k2 2->1 (1/ms)", params_21[curr_synapse_idx].k2);
      setParameter("s_slow 2->1 (1/mV)", params_21[curr_synapse_idx].s_slow);
      setParameter("V_slow 2->1 (mV)", params_21[curr_synapse_idx].v_slow);
    }
  }
}

void BidirectionalChemicalSynapseBO::set_evaluations_completed(double evals)
{
  evaluations_completed = evals;
}
