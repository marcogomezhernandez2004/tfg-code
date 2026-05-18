#ifndef SCALING_UTILS_H
#define SCALING_UTILS_H

#include <limits>
#include <array>
#include <cstddef>
#include <concepts>
#include <type_traits>
#include <cmath>

#include <DifferentialNeuronWrapper.h>
#include <HindmarshRoseModel.h>
#include <SystemWrapper.h>
#include <RungeKutta4.h>

namespace BOPublicConfig
{
    // Porcentaje de margen sobre el rango esperado de corriente para ampliar la zona de búsqueda
    inline constexpr double EXPECTED_I_MARGIN_FACTOR = 0.5;
}

namespace BOPublicConstants
{
    // Umbral mínimo para divisores; evita divisiones por cero
    inline constexpr double SMALL_DIVISOR = std::numeric_limits<double>::epsilon();
    inline constexpr double NEGATIVE_SMALL_DIVISOR = -SMALL_DIVISOR;
}

// Parámetros del modelo de neurona Hindmarsh-Rose
namespace HindmarshRoseParams
{

    inline constexpr double e = 3.281;
    inline constexpr double mu = 0.0021;
    inline constexpr double S = 1.0;
    inline constexpr double a = 1.0;
    inline constexpr double b = 3.0;
    inline constexpr double c = 1.0;
    inline constexpr double d = 5.0;
    inline constexpr double xr = -1.6;
    inline constexpr double vh = 0.1;
}

// Estado inicial de la neurona Hindmarsh-Rose (valores de x, y, z ya estabilizados)
namespace HindmarshRoseInitialState
{

    inline constexpr double x = -0.712841;
    inline constexpr double y = -1.93688;
    inline constexpr double z = 3.16568;
}

// Alias: neurona HR envuelta con un integrador numérico genérico
template <typename Integrator>
using HindmarshRoseNeuron = DifferentialNeuronWrapper<SystemWrapper<HindmarshRoseModel<double>>, Integrator>;

// Componente sináptica a optimizar: solo ifast, solo islow, o ambas
enum SynComponent
{
    IFAST = 0,
    ISLOW = 1,
    BOTH = 2
};

namespace GeneralConstants
{

    inline constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();
    inline constexpr double DOUBLE_MIN = std::numeric_limits<double>::lowest();
}

namespace SigConstants
{
    // Valores centinela para dt y pts inválidos
    inline constexpr double INVALID_DT = -1.0;
    inline constexpr double INVALID_PTS = -1.0;
}

namespace SigPublicConfig
{
    // Porcentajes del rango de la señal para definir los umbrales de detección de burst
    // th_on = SIG_PERCENTAGE_MIN * rango + min,  th_up = SIG_PERCENTAGE_MAX * rango + min
    inline constexpr double SIG_PERCENTAGE_MIN = 0.10;
    inline constexpr double SIG_PERCENTAGE_MAX = 0.90;
}

// Concepto: función que crea una neurona; recibe bool (true = vacía, false = con params)
template <typename F, typename NeuronType>
concept CreateFunc = std::invocable<F, bool> && std::convertible_to<std::invoke_result_t<F, bool>, NeuronType>;

// Concepto: función que resetea el estado de una neurona
template <typename F, typename NeuronType>
concept ResetStateFunc = std::invocable<F, NeuronType &>;

// Concepto: función que extrae el voltaje (variable x en HR) de una neurona
template <typename F, typename NeuronType>
concept GetVFunc = std::invocable<F, const NeuronType &> && std::convertible_to<std::invoke_result_t<F, const NeuronType &>, double>;

// Resetea el estado de la neurona HR a los valores iniciales predefinidos
template <typename Integrator>
inline void reset_state_hindmarsh_rose(HindmarshRoseNeuron<Integrator> &neuron)
{

    using NeuronType = HindmarshRoseNeuron<Integrator>;
    neuron.set(NeuronType::x, HindmarshRoseInitialState::x);
    neuron.set(NeuronType::y, HindmarshRoseInitialState::y);
    neuron.set(NeuronType::z, HindmarshRoseInitialState::z);
    neuron.reset_synaptic_input();
}

// Devuelve el potencial de membrana (variable x) de la neurona HR
template <typename Integrator>
double get_v_hindmarsh_rose(const HindmarshRoseNeuron<Integrator> &neuron)
{

    return neuron.get(HindmarshRoseNeuron<Integrator>::x);
}

// Crea una neurona HR; si empty=true no asigna parámetros (neurona "vacía" para la post-sináptica)
template <typename Integrator>
HindmarshRoseNeuron<Integrator> create_hindmarsh_rose(bool empty)
{

    using NeuronType = HindmarshRoseNeuron<Integrator>;
    typename NeuronType::ConstructorArgs args;

    if (!empty)
    {
        args.params[NeuronType::e] = HindmarshRoseParams::e;
        args.params[NeuronType::mu] = HindmarshRoseParams::mu;
        args.params[NeuronType::S] = HindmarshRoseParams::S;
        args.params[NeuronType::a] = HindmarshRoseParams::a;
        args.params[NeuronType::b] = HindmarshRoseParams::b;
        args.params[NeuronType::c] = HindmarshRoseParams::c;
        args.params[NeuronType::d] = HindmarshRoseParams::d;
        args.params[NeuronType::xr] = HindmarshRoseParams::xr;
        args.params[NeuronType::vh] = HindmarshRoseParams::vh;
    }

    return NeuronType(args);
}

// Divisor seguro: si |divisor| < epsilon, devuelve +-epsilon para evitar div/0
inline double safe_divisor(double divisor)
{

    return std::abs(divisor) < BOPublicConstants::SMALL_DIVISOR
               ? (divisor < 0.0 ? BOPublicConstants::NEGATIVE_SMALL_DIVISOR : BOPublicConstants::SMALL_DIVISOR)
               : divisor;
}

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

namespace HindmarshRose
{
    // Rango dinámico del modelo HR (mín y máx de la variable x) precalculado por consts_calculator
    inline constexpr double MIN = -1.668473;
    inline constexpr double MAX = 1.764310;
    // Tablas precalculadas: DTS_RK4 = pasos de tiempo posibles,
    // PTS_RK4 = puntos promedio por burst para cada dt (usadas para elegir el dt del modelo)
    inline constexpr std::array<double, 144> DTS_RK4 = {
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
    inline constexpr std::array<double, 144> PTS_RK4 = {
        1689645.000000, 1408038.000000, 1206890.000000, 1056028.000000, 938692.000000, 844822.500000, 768021.000000, 704019.000000,
        649863.500000, 603445.000000, 563215.000000, 528014.000000, 469346.000000, 422411.000000, 384010.500000, 337929.000000,
        301722.500000, 291318.000000, 281607.500000, 272523.500000, 264007.000000, 256007.000000, 248477.000000, 241378.000000,
        234673.000000, 228330.500000, 222322.000000, 216621.000000, 211205.500000, 206054.500000, 201148.500000, 196470.500000,
        192005.000000, 187738.500000, 183657.000000, 179749.500000, 176005.000000, 172412.500000, 168964.500000, 165651.500000,
        162466.000000, 156448.500000, 150861.500000, 145659.000000, 140803.500000, 136262.000000, 132003.500000, 128003.500000,
        124238.500000, 120689.000000, 117336.500000, 114165.000000, 109717.000000, 105602.500000, 101786.000000, 98235.000000,
        94923.500000, 91828.500000, 88002.500000, 84482.000000, 81233.000000, 77506.500000, 74107.500000, 70993.500000,
        67585.500000, 64490.500000, 61219.000000, 57864.500000, 54858.500000, 51829.500000, 48834.000000, 45666.000000,
        42453.500000, 39294.000000, 36258.500000, 33130.000000, 30064.500000, 29747.000000, 29436.000000, 29132.000000,
        28735.500000, 28349.500000, 27974.000000, 27608.500000, 27252.000000, 26905.000000, 26566.500000, 26237.000000,
        25915.000000, 25600.500000, 25294.000000, 24920.500000, 24558.500000, 24207.000000, 23865.000000, 23532.500000,
        23209.000000, 22895.000000, 22589.000000, 22232.000000, 21886.500000, 21551.500000, 21226.500000, 20911.500000,
        20605.000000, 20259.500000, 19925.000000, 19601.500000, 19287.500000, 18984.500000, 18649.000000, 18326.000000,
        18013.000000, 17711.000000, 17382.500000, 17066.500000, 16762.000000, 16436.000000, 16122.500000, 15820.000000,
        15500.500000, 15194.000000, 14873.000000, 14565.500000, 14246.000000, 13940.000000, 13625.000000, 13324.500000,
        13016.500000, 12703.000000, 12386.500000, 12085.000000, 11781.500000, 11477.500000, 11173.000000, 10871.000000,
        10571.500000, 10263.000000, 9960.500000, 9652.500000, 9353.000000, 9052.000000, 8751.500000, 8444.500000};
};

#endif
