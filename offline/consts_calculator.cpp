#include <iostream>
#include <vector>
#include <array>
#include <iomanip>
#include <string>
#include <stdexcept>
#include "scaling.hpp"
#include "utils.hpp"

// Programa auxiliar que precalcula las constantes del modelo neuronal
// (min, max, y puntos promedio por burst para cada dt) que se usan en utils.hpp.
// Se ejecuta una sola vez y su salida se copia al código fuente como tablas constantes.

namespace ConstCalculatorConstants
{
    // Número de bursts sobre los que promediar para calcular puntos/burst
    static constexpr int BURSTS_TO_AVERAGE = 20;
}

namespace ConstCalculatorConfig
{
    // Array de dts a evaluar (de fino a grueso)
    static constexpr std::array<double, 144> DTS = {
        0.000500, 0.000600, 0.000700, 0.000800, 0.000900, 0.001000, 0.001100, 0.001200,
        0.001300, 0.001400, 0.001500, 0.001600, 0.001800, 0.002000, 0.002200, 0.002500,
        0.002800, 0.002900, 0.003000, 0.003100, 0.003200, 0.003300, 0.003400, 0.003500,
        0.003600, 0.003700, 0.003800, 0.003900, 0.004000, 0.004100, 0.004200, 0.004300,
        0.004400, 0.004500, 0.004600, 0.004700, 0.004800, 0.004900, 0.005000, 0.005100,
        0.005200, 0.005400, 0.005600, 0.005800, 0.006000, 0.006200, 0.006400, 0.006600,
        0.006800, 0.007000, 0.007200, 0.007400, 0.007700, 0.008000, 0.008300, 0.008600,
        0.008900, 0.009200, 0.009600, 0.010000, 0.010400, 0.010900, 0.011400, 0.011900,
        0.012500, 0.013100, 0.013800, 0.014600, 0.015400, 0.016300, 0.017300, 0.018500,
        0.019900, 0.021500, 0.023300, 0.025500, 0.028100, 0.028400, 0.028700, 0.029000,
        0.029400, 0.029800, 0.030200, 0.030600, 0.031000, 0.031400, 0.031800, 0.032200,
        0.032600, 0.033000, 0.033400, 0.033900, 0.034400, 0.034900, 0.035400, 0.035900,
        0.036400, 0.036900, 0.037400, 0.038000, 0.038600, 0.039200, 0.039800, 0.040400,
        0.041000, 0.041700, 0.042400, 0.043100, 0.043800, 0.044500, 0.045300, 0.046100,
        0.046900, 0.047700, 0.048600, 0.049500, 0.050400, 0.051400, 0.052400, 0.053400,
        0.054500, 0.055600, 0.056800, 0.058000, 0.059300, 0.060600, 0.062000, 0.063400,
        0.064900, 0.066500, 0.068200, 0.069900, 0.071700, 0.073600, 0.075600, 0.077700,
        0.079900, 0.082300, 0.084800, 0.087500, 0.090300, 0.093300, 0.096500, 0.100000};

    static constexpr double OBSERVATION_TIME = 2000.0;
    static constexpr double MINMAX_DT = DTS[0]; // El dt más fino para calcular min/max precisos
    static constexpr double STABILIZATION_TIME = 2000.0;
}

static constexpr size_t DTS_SIZE = ConstCalculatorConfig::DTS.size();

struct MinMaxResult
{

    double min;
    double max;
};

struct PtsResult
{

    std::array<double, DTS_SIZE> pts; // Puntos promedio por burst para cada dt
    std::vector<double> invalid_dts;  // dts donde no se detectaron bursts
};

// Simula la neurona con el dt más fino y registra su rango dinámico [min, max].
// Se usa el dt más pequeño disponible para obtener la máxima precisión en los extremos.
// La estabilización previa evita incluir el transitorio inicial del modelo.
template <typename NeuronType, CreateFunc<NeuronType> CreateFuncType, ResetStateFunc<NeuronType> ResetFuncType, GetVFunc<NeuronType> GetVFuncType>
static inline MinMaxResult calculate_min_max(
    CreateFuncType create_neur,
    ResetFuncType reset_state_neur,
    GetVFuncType get_v_neur,
    double observation_time,
    double dt,
    double stabilization_time)
{

    if (observation_time <= 0 || dt <= 0 || stabilization_time < 0)
    {
        throw std::invalid_argument("observation_time and dt must be positive, stabilization_time non-negative");
    }

    NeuronType neuron = create_neur(false);

    reset_state_neur(neuron);

    // Fase de estabilización: avanza la neurona hasta que el transitorio inicial
    // se disipe y la dinámica alcance el ciclo límite (atractor bursting)
    const size_t stabilization_steps = static_cast<size_t>(stabilization_time / dt);

    for (size_t i = 0; i < stabilization_steps; i++)
        neuron.step(dt);

    // Inicializa min/max con valores extremos opuestos (sentinelas):
    // min = +inf para que cualquier valor real sea menor, max = -inf analogamente
    double min = GeneralConstants::DOUBLE_MAX;
    double max = GeneralConstants::DOUBLE_MIN;
    const size_t obs_steps = static_cast<size_t>(observation_time / dt);

    // Fase de observación: registra los extremos del voltaje en el ciclo límite
    for (size_t step = 0; step < obs_steps; step++)
    {
        neuron.step(dt);

        double val = get_v_neur(neuron);

        if (val > max)
            max = val;
        if (val < min)
            min = val;
    }

    MinMaxResult result;
    result.min = min;
    result.max = max;

    return result;
}

// Para cada dt, simula la neurona y mide el número promedio de pasos de integración
// que componen un burst completo. Esto produce la tabla PTS que se usa en select_dt.
//
// Método:
//   1. Detecta bursts con histéresis (th_on = 10% del rango, th_up = 90%).
//   2. Un burst se define como el intervalo entre dos cruces ascendentes de th_up.
//   3. Se promedian BURSTS_TO_AVERAGE bursts para reducir variabilidad.
//   4. bursts_seen empieza en -1: el primer cruce ascendente marca el inicio del
//      primer burst pero NO contar como burst completo (el conteo de pasos empieza ahí).
template <typename NeuronType, size_t N, CreateFunc<NeuronType> CreateFuncType, ResetStateFunc<NeuronType> ResetFuncType, GetVFunc<NeuronType> GetVFuncType>
static inline PtsResult calculate_pts(
    CreateFuncType create_neur,
    ResetFuncType reset_state_neur,
    GetVFuncType get_v_neur,
    const std::array<double, N> &dts,
    double observation_time,
    double stabilization_time,
    double min_val,
    double max_val)
{

    NeuronType neuron = create_neur(false);

    PtsResult result;

    const double range = max_val - min_val;

    // Umbrales de detección de burst por histéresis (mismas constantes que scaling.cpp)
    // th_on: señal "activa" (baja, para permitir nuevo cruce); th_up: señal "arriba" (alta)
    const double th_on = SigPublicConfig::SIG_PERCENTAGE_MIN * range + min_val;
    const double th_up = SigPublicConfig::SIG_PERCENTAGE_MAX * range + min_val;

    std::array<double, N> &pts = result.pts;
    std::vector<double> &invalid_dts = result.invalid_dts;

    for (size_t i = 0; i < N; i++)
    {
        const double dt = dts[i];

        // Resetea la neurona para cada dt (independizar las mediciones)
        reset_state_neur(neuron);

        // Estabilización: deja que el transitorio se disipe antes de medir
        const size_t stabilization_steps = static_cast<size_t>(stabilization_time / dts[i]);

        for (size_t j = 0; j < stabilization_steps; j++)
            neuron.step(dt);

        bool up = (get_v_neur(neuron) > th_up);
        double total_steps = 0;

        // bursts_seen = -1: el primer cruce ascendente (burst 0) solo sirve para marcar
        // el inicio del conteo. Los pasos acumulados antes de ese primer cruce se descartan
        // porque corresponden a un burst parcial (empezamos a contar en un punto arbitrario).
        int bursts_seen = -1;
        // Contador de pasos entre dos cruces ascendentes consecutivos = 1 burst completo
        size_t steps_in_current_burst = 0;
        double act_time = 0.0;

        while (bursts_seen < ConstCalculatorConstants::BURSTS_TO_AVERAGE && act_time < observation_time)
        {
            neuron.step(dt);
            act_time += dt;

            double val = get_v_neur(neuron);

            // Cruce ascendente del umbral superior -> inicio de un nuevo burst
            if (!up && val > th_up)
            {

                up = true;
                bursts_seen++;

                // Acumula los pasos del burst que acaba de terminar.
                // Cuando bursts_seen pasa de -1 a 0, total_steps += 0 (descarta el parcial).
                // A partir de bursts_seen >= 1, total_steps acumula bursts completos.
                total_steps += steps_in_current_burst;
                // Reinicia el contador para el siguiente burst
                steps_in_current_burst = 0;
            }
            else if (up && val < th_on)
            {
                // Cruce descendente: habilita futura detección de nuevo cruce ascendente
                up = false;
            }

            // Se incrementa siempre, incluyendo paso actual (cuenta todos los pasos del burst)
            steps_in_current_burst++;
        }

        if (bursts_seen <= 0)
        {
            // No se detectaron bursts completos con este dt -> marca como inválido
            pts[i] = GeneralConstants::DOUBLE_MAX;
            invalid_dts.push_back(dts[i]);
        }
        else
        {
            // Promedio: total_steps contiene la suma de pasos de bursts_seen bursts completos
            pts[i] = total_steps / static_cast<double>(bursts_seen);
        }
    }

    return result;
}

static inline std::string integrator_to_string(NumericIntegrator integrator)
{

    switch (integrator)
    {
    case RK4:
        return "RK4";
    default:
        return "UNKNOWN";
    }
}

// Imprime las tablas en formato C++ para copiar a utils.hpp
static inline void print_tables(const PtsResult &pr, NumericIntegrator integrator)
{

    const std::string integrator_str = integrator_to_string(integrator);

    const size_t ds_size_minus_1 = DTS_SIZE - 1;
    std::cout << "inline constexpr std::array<double, " << DTS_SIZE << "> DTS_" << integrator_str << " = {";
    for (size_t i = 0; i < DTS_SIZE; i++)
    {

        if (i % 8 == 0)
            std::cout << "\n    ";
        std::cout << ConstCalculatorConfig::DTS[i];
        if (i < ds_size_minus_1)
            std::cout << ", ";
    }
    std::cout << "};\n";

    const std::array<double, DTS_SIZE> &pts = pr.pts;
    std::cout << "inline constexpr std::array<double, " << DTS_SIZE << "> PTS_" << integrator_str << " = {";
    for (size_t i = 0; i < DTS_SIZE; i++)
    {

        if (i % 8 == 0)
            std::cout << "\n    ";
        std::cout << pts[i];
        if (i < ds_size_minus_1)
            std::cout << ", ";
    }
    std::cout << "};\n";

    const std::vector<double> &invalid_dts = pr.invalid_dts;
    if (!invalid_dts.empty())
    {
        std::cout << "Invalid dts (no bursts detected): ";
        for (double dt : invalid_dts)
        {
            std::cout << dt << " ";
        }
        std::cout << "\n";
    }
}

typedef RungeKutta4 Integrator;
typedef HindmarshRoseNeuron<Integrator> NeuronType;

int main()
{

    NumericIntegrator integrator = RK4;

    CreateFunc<NeuronType> auto create_func = &create_hindmarsh_rose<Integrator>;
    ResetStateFunc<NeuronType> auto reset_func = &reset_state_hindmarsh_rose<Integrator>;
    GetVFunc<NeuronType> auto get_v_func = &get_v_hindmarsh_rose<Integrator>;

    constexpr double OBSERVATION_TIME = ConstCalculatorConfig::OBSERVATION_TIME;
    constexpr double STABILIZATION_TIME = ConstCalculatorConfig::STABILIZATION_TIME;

    // 1. Calcula el rango dinámico [min, max] de la neurona con el dt más fino
    MinMaxResult mmr = calculate_min_max<NeuronType>(
        create_func,
        reset_func,
        get_v_func,
        OBSERVATION_TIME,
        ConstCalculatorConfig::MINMAX_DT,
        STABILIZATION_TIME);

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "inline constexpr double MIN = " << mmr.min << ";\n";
    std::cout << "inline constexpr double MAX = " << mmr.max << ";\n";

    // 2. Calcula los puntos promedio por burst para cada dt
    PtsResult pr = calculate_pts<NeuronType>(
        create_func,
        reset_func,
        get_v_func,
        ConstCalculatorConfig::DTS,
        OBSERVATION_TIME,
        STABILIZATION_TIME,
        mmr.min,
        mmr.max);

    print_tables(pr, integrator);

    return 0;
}
