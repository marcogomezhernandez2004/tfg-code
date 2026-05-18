#include "scaling.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "utils.hpp"

using namespace kfr;

namespace SigPrivateConfig
{
    // Tolerancia al seleccionar dt: se permite que la parte fraccionaria de
    // (pts_burst_modelo / pts_real) sea hasta este factor * parte entera
    static constexpr double DT_SELECTION_TOLERANCE = 0.1;

    // Cada cuántos bursts (en puntos) se recalcula el escalado por drift
    static constexpr size_t DRIFT_N_BURST = 2;
}

// Factores de escalado lineal: valor_virtual = valor_real * scale + offset
struct ScalingFactors
{

    double scale_real_to_virtual;
    double offset_real_to_virtual;
};

// dt seleccionado del modelo y sus puntos promedio por burst asociados
struct DTSelection
{

    double dt;
    double pts_burst;
};

// Estadísticas de la señal real (del CSV)
struct SigStats
{

    double real_abs_min; // Mínimo absoluto de la señal
    double real_abs_max; // Máximo absoluto de la señal
    double real_rel_min; // Umbral inferior de detección de burst (10% del rango)
    double real_rel_max; // Umbral superior de detección de burst (90% del rango)
    double sig_period;   // Período medio de la señal (tiempo entre bursts)
};

// Lee num_points valores de la columna column_idx de un CSV, empezando en start_idx
static inline void read_csv_column(univector<double> &data, const std::string &csv_path, size_t column_idx,
                                   size_t start_idx, size_t num_points)
{

    std::ifstream file(csv_path);

    if (!file.is_open())
    {
        throw std::runtime_error("Unable to open CSV file: " + csv_path);
    }

    data.reserve(num_points);

    const size_t end_idx = start_idx + num_points;

    std::string line;
    std::string val;
    std::stringstream ss;
    size_t current_line = 0;

    while (current_line < end_idx && std::getline(file, line))
    {
        if (current_line >= start_idx)
        {

            ss.str(line);
            ss.clear();
            size_t current_col = 0;

            while (std::getline(ss, val, ','))
            {
                if (current_col == column_idx)
                {

                    data.push_back(std::stod(val));
                    break;
                }
                current_col++;
            }
        }

        current_line++;
    }

    const size_t data_size = data.size();

    if (data_size < num_points)
    {
        data.shrink_to_fit();
        std::cout << "Warning: Fewer data points read (" << data_size << ") than expected (" << num_points << ") from CSV file: " << csv_path << std::endl;
    }

    file.close();
}

// Calcula el período medio de la señal (en tiempo) contando los cruces ascendentes del umbral th_up.
// Usa histéresis con dos umbrales para evitar detecciones espurias por ruido:
//   - th_up (90% del rango): umbral de cruce ascendente -> marca inicio de burst
//   - th_on (10% del rango): umbral de cruce descendente -> permite un nuevo cruce ascendente
// Sin histéresis, oscilaciones de ruido cerca del umbral generarían falsos cruces.
static inline double sig_period(double observation_time, const univector_ref<const double> &sig,
                                double th_up, double th_on)
{

    // Estado inicial: determina si la señal empieza "arriba" (ya dentro de un burst)
    bool up = (sig.front() > th_up);

    double changes = 0.0;

    for (const double val : sig)
    {
        if (!up && val > th_up)
        {
            // Cruce ascendente: la señal sube por encima de th_up -> nuevo burst detectado
            changes++;
            up = true;
        }
        else if (up && val < th_on)
        {
            // Cruce descendente: la señal cae por debajo de th_on -> sale del burst.
            // Solo entonces se permite detectar otro cruce ascendente (histéresis).
            up = false;
        }
    }

    // frecuencia_bursts = changes / observation_time  (bursts por unidad de tiempo)
    // Período = 1 / frecuencia = observation_time / changes
    return 1.0 / (changes / observation_time);
}

// Calcula factores de escalado lineal del espacio real al virtual (modelo neuronal).
// Mapea [real_min, real_max] -> [virtual_min, virtual_max] con: v = r * scale + offset
// Derivación del offset:
//   virtual_min = real_min * scale + offset  ->  offset = virtual_min - real_min * scale
static ScalingFactors calculate_scaling(double virtual_min, double virtual_max,
                                        double real_min, double real_max)
{

    const double virtual_range = virtual_max - virtual_min;
    const double real_range = real_max - real_min;

    ScalingFactors factors;
    // scale = rango_virtual / rango_real (factor de proporcionalidad)
    const double scale_real_to_virtual = virtual_range / real_range;
    factors.scale_real_to_virtual = scale_real_to_virtual;
    // offset para que real_min se mapee exactamente a virtual_min
    factors.offset_real_to_virtual = virtual_min - (real_min * scale_real_to_virtual);

    return factors;
}

// Selecciona el dt del modelo neuronal cuyo burst (en puntos) sea múltiplo cercano de pts_real.
//
// Contexto: la señal del CSV tiene un cierto período de burst (medido en muestras = pts_real).
// El modelo neuronal, simulado con un determinado dt, produce bursts de pts[i] pasos.
// Necesitamos un dt tal que un burst del modelo contenga un número entero (o casi) de
// bursts de la señal real -> así se alinean temporalmente.
//
// Algoritmo:
//   1. Multiplica pts_real por factor creciente (1, 2, 3, ...) hasta alcanzar el rango de pts[]
//   2. Para cada factor, busca de dt grande (computacionalmente barato) a dt pequeño (caro)
//      el primer dt cuyo pts[i] > pts_real*factor
//   3. Comprueba divisibilidad: pts[i] / pts_real ~ entero (parte fraccional <= tolerancia * parte entera)
//   4. Si encuentra uno con buena divisibilidad, lo acepta. Si no, sigue con el siguiente factor.
//   5. Fallback: si ninguno tuvo buena divisibilidad, toma el último candidato encontrado.
template <size_t N>
static inline std::optional<DTSelection> select_dt_neuron_model(const std::array<double, N> &dts,
                                                                const std::array<double, N> &pts,
                                                                double pts_real)
{

    double aux = pts_real;

    double factor = 1.0;
    double intpart, fractpart;
    bool flag = false;

    constexpr double INVALID_DT = SigConstants::INVALID_DT;

    double dt_candidate = INVALID_DT;
    double pts_burst_candidate = SigConstants::INVALID_PTS;

    // Itera con factor creciente: busca que pts_real*factor sea alcanzable por algún pts[i].
    // pts[0] corresponde al dt más fino (máximos pts por burst). Si aux >= pts[0],
    // no tiene sentido seguir porque ningún dt producirá tantos puntos.
    while (aux < pts[0])
    {

        aux = pts_real * factor;
        factor += 1.0;

        // Recorre las tablas en orden inverso (dt grande -> dt pequeño, es decir,
        // pts pequeño -> pts grande). Busca el primer dt cuyo pts[i] supere aux.
        // Esto prioriza dts grandes (simulación más barata) que aún cubran el burst.
        for (size_t i = N - 1; i >= 0; i--)
        {
            if (pts[i] > aux)
            {

                dt_candidate = dts[i];
                pts_burst_candidate = pts[i];

                // Test de divisibilidad: descomponemos pts_burst / pts_real en parte
                // entera y fraccionaria. Si la parte fraccionaria es pequeña respecto
                // a la entera, los bursts se alinean bien (ej: 3.02 -> intpart=3, fract=0.02).
                // Criterio: fract <= TOLERANCE * intpart (ej: 0.02 <= 0.1 * 3 = 0.3 -> OK)
                fractpart = std::modf(pts_burst_candidate / pts_real, &intpart);

                if (fractpart <= SigPrivateConfig::DT_SELECTION_TOLERANCE * intpart)
                {
                    flag = true;
                }

                break;
            }
        }

        if (flag)
        {
            break;
        }
    }

    // Fallback: si ningún factor produjo buena divisibilidad, acepta el último candidato.
    // Es posible que con todos los factores probados, el mejor dt sea "aceptable" aunque
    // no ideal. Repite la búsqueda con el último aux para no quedarse sin resultado.
    if (!flag)
    {
        for (size_t i = N - 1; i >= 0; i--)
        {
            if (pts[i] > aux)
            {

                dt_candidate = dts[i];
                pts_burst_candidate = pts[i];
                break;
            }
        }
    }

    if (dt_candidate == INVALID_DT)
    {
        return std::nullopt;
    }

    DTSelection selection;
    selection.dt = dt_candidate;
    selection.pts_burst = pts_burst_candidate;

    return selection;
}

// Despacha la selección de dt según modelo e integrador
static inline std::optional<DTSelection> set_pts_burst(NeuronModel model, NumericIntegrator integrator, double pts_real)
{

    if (model == NeuronModel::HINDMARSH_ROSE)
    {
        if (integrator == NumericIntegrator::RK4)
        {
            return select_dt_neuron_model(HindmarshRose::DTS_RK4, HindmarshRose::PTS_RK4, pts_real);
        }
        else
        {
            throw std::invalid_argument("Unsupported integrator");
        }
    }
    else
    {
        throw std::invalid_argument("Unsupported neuron model");
    }
}

// Inicializa las estadísticas de la señal: extremos, umbrales de burst y período
static inline SigStats init_stats(const univector<double> &sig, size_t obs_points, double csv_step)
{

    const double observation_time_to_use = obs_points * csv_step;

    SigStats stats;

    const univector_ref<const double> obs_sig = sig.slice(0, obs_points);
    const double abs_max = maxof(obs_sig);
    const double abs_min = minof(obs_sig);

    stats.real_abs_min = abs_min;
    stats.real_abs_max = abs_max;

    const double range = abs_max - abs_min;

    // Umbrales de detección de burst basados en porcentajes del rango
    const double real_rel_min = SigPublicConfig::SIG_PERCENTAGE_MIN * range + abs_min;
    const double real_rel_max = SigPublicConfig::SIG_PERCENTAGE_MAX * range + abs_min;
    stats.real_rel_min = real_rel_min;
    stats.real_rel_max = real_rel_max;

    stats.sig_period = sig_period(observation_time_to_use, obs_sig, real_rel_max, real_rel_min);

    return stats;
}

// Función principal de escalado: lee señal del CSV, la escala al rango del modelo neuronal,
// selecciona el dt adecuado y genera puntos interpolados para sub-stepping
std::optional<ScaledSigResult> scale_sig(
    const std::string &csv_path,
    size_t column_idx,
    double csv_step,
    double start_time,
    double use_time,
    double observation_time,
    NumericIntegrator integrator,
    NeuronModel model,
    bool check_drift)
{

    if (csv_step <= 0 || use_time <= 0 || observation_time <= 0 || start_time < 0 || column_idx < 0 || csv_path.empty())
    {
        throw std::invalid_argument("Invalid arguments: csv_step, use_time, observation_time must be positive, start_time and column_idx non-negative, csv_path non-empty");
    }

    ScaledSigResult result;

    univector<double> &sig = result.sig;
    univector<double> &interpolated_points = result.interpolated_points;

    const size_t start_idx = static_cast<size_t>(start_time / csv_step);
    const size_t use_points = static_cast<size_t>(use_time / csv_step);
    if (use_points == 0)
    {
        throw std::invalid_argument("use_time is too short to read any points with given csv_step");
    }

    const size_t obs_points = static_cast<size_t>(observation_time / csv_step);
    if (obs_points == 0)
    {
        throw std::invalid_argument("observation_time is too short to read any points with given csv_step");
    }

    // Lee el máximo entre puntos de uso y de observación (la observación se usa para estadísticas)
    const size_t read_points = std::max(use_points, obs_points);

    read_csv_column(sig, csv_path, column_idx, start_idx, read_points);

    size_t sig_size = sig.size();
    if (sig_size == 0)
    {
        throw std::runtime_error("No data read from CSV file");
    }

    SigStats stats = init_stats(sig, obs_points, csv_step);

    // Recorta la señal al tamaño de uso si se leyó más (por observación)
    if (sig_size > use_points)
    {
        sig.resize(use_points);
        sig_size = use_points;
    }

    // Puntos por burst en la señal externa (CSV)
    const double external_pts_per_burst = stats.sig_period / csv_step;

    // Busca el dt del modelo cuyo burst coincida en período con la señal externa
    std::optional<DTSelection> selection;
    double model_abs_min, model_abs_max;
    if (model == NeuronModel::HINDMARSH_ROSE)
    {

        model_abs_min = HindmarshRose::MIN;
        model_abs_max = HindmarshRose::MAX;
        selection = set_pts_burst(model, integrator, external_pts_per_burst);
    }
    else
    {
        throw std::invalid_argument("Unsupported neuron model");
    }

    if (!selection)
    {
        // No se encontró un dt compatible con la frecuencia de burst de la señal
        return std::nullopt;
    }

    result.dt = selection->dt;

    // Factor de sub-pasos: cuántos pasos del modelo se ejecutan por cada muestra del CSV.
    // Ej: si el modelo con este dt necesita 84000 pasos/burst y la señal tiene 42000 muestras/burst,
    // s_points = 84000/42000 = 2 -> por cada muestra del CSV se dan 2 pasos del modelo.
    // Esto permite que el modelo avance a su ritmo natural mientras procesa la señal externa.
    size_t s_points = static_cast<size_t>(selection->pts_burst / external_pts_per_burst);

    if (s_points == 0)
        s_points = 1;

    result.points_factor = s_points;

    double &real_abs_min = stats.real_abs_min;
    double &real_abs_max = stats.real_abs_max;

    // Escalado lineal: mapea [real_min, real_max] -> [model_min, model_max]
    ScalingFactors factors = calculate_scaling(model_abs_min, model_abs_max, real_abs_min, real_abs_max);
    double scale_real_to_virtual = factors.scale_real_to_virtual;
    double offset_real_to_virtual = factors.offset_real_to_virtual;

    if (check_drift)
    {
        // Corrección de drift no causal: la señal biológica puede derivar lentamente.
        // Dado que se dispone de toda la señal offline, se procesa en ventanas temporales.
        // En cada ventana se calculan los extremos locales (filtrando outliers) y
        // se aplican los factores de escalado derivados de ellos a los puntos de esa misma ventana.

        constexpr double DOUBLE_MAX = GeneralConstants::DOUBLE_MAX;
        constexpr double DOUBLE_MIN = GeneralConstants::DOUBLE_MIN;
        // Inicializa min/max invertidos para que la primera comparación siempre actualice
        double window_max, window_min;
        // Rango de la señal global: se usa para filtrar outliers (artefactos de medida)
        const double drift_aux_range = real_abs_max - real_abs_min;
        const size_t window_size = static_cast<size_t>(SigPrivateConfig::DRIFT_N_BURST * external_pts_per_burst);

        for (size_t start = 0; start < sig_size; start += window_size)
        {
            size_t end = std::min(start + window_size, sig_size);

            window_min = DOUBLE_MAX;
            window_max = DOUBLE_MIN;

            univector_ref<double> sig_window = sig.slice(start, end - start);

            for (double &val : sig_window)
            {
                // Actualiza window_min solo si val es un nuevo mínimo Y no es un outlier.
                // Filtro de outliers: descarta valores que caigan más de un rango completo
                // por debajo del mínimo global. Ej: si min=-2 y rango=3, descarta val < -5.
                // Esto evita que artefactos puntuales distorsionen el escalado.
                if ((window_min > val) && (val > (real_abs_min - drift_aux_range)))
                {
                    window_min = val;
                }
                // Análogo para window_max: descarta outliers por encima de max + rango
                if ((window_max < val) && (val < (real_abs_max + drift_aux_range)))
                {
                    window_max = val;
                }
            }

            if (window_max != DOUBLE_MIN && window_min != DOUBLE_MAX)
            {
                // Recalcula scale/offset adaptados al drift actual
                factors = calculate_scaling(model_abs_min, model_abs_max, window_min, window_max);
                scale_real_to_virtual = factors.scale_real_to_virtual;
                offset_real_to_virtual = factors.offset_real_to_virtual;
            }

            // Aplica el escalado actual (posiblemente actualizado por drift): real -> virtual
            process(sig_window, (sig_window * scale_real_to_virtual) + offset_real_to_virtual);
        }
    }
    else
    {
        // Sin drift: aplica escalado lineal uniforme a toda la señal
        sig = (sig * scale_real_to_virtual) + offset_real_to_virtual;
    }

    // Genera puntos interpolados linealmente entre muestras consecutivas de sig
    // para sub-stepping del modelo sináptico.
    // Ej con s_points=3: entre sig[i] y sig[i+1] se generan 2 puntos intermedios:
    //   interp_1 = sig[i] + (1/3) * (sig[i+1] - sig[i])  -> a 1/3 del camino
    //   interp_2 = sig[i] + (2/3) * (sig[i+1] - sig[i])  -> a 2/3 del camino
    // La secuencia total por muestra sería: sig[i], interp_1, interp_2, sig[i+1], ...
    // En evaluate_candidate: se usa sig[i] como primer paso y luego los interpolados.
    const size_t interpolated_size = (sig_size - 1) * (s_points - 1);
    interpolated_points.reserve(interpolated_size);
    const double *sig_ptr = sig.data();

    for (size_t i = 0; i < sig_size - 1; i++)
    {
        // j va de 1 a s_points-1 (j=0 sería sig[i], que ya existe; j=s_points sería sig[i+1])
        for (double j = 1.0; j < s_points; j++)
        {
            // alpha in (0, 1): fracción del intervalo entre sig[i] y sig[i+1]
            double alpha = j / s_points;
            // Interpolación lineal: sig[i] * (1 - alpha) + sig[i+1] * alpha
            // Reescrita como: sig[i] + alpha * (sig[i+1] - sig[i])
            double interp_val = sig_ptr[i] + (alpha * (sig_ptr[i + 1] - sig_ptr[i]));
            interpolated_points.push_back(interp_val);
        }
    }

    return result;
}
