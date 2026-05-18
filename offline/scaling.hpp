#ifndef SCALING_H
#define SCALING_H

#include <cstddef>
#include <optional>
#include <string>
#include <kfr/all.hpp>
using namespace kfr;

enum NumericIntegrator
{
    RK4
};

enum NeuronModel
{
    HINDMARSH_ROSE
};

// Resultado de escalar la señal del CSV al espacio del modelo neuronal
struct ScaledSigResult
{
    // Señal ya escalada al rango dinámico del modelo (ej: HR [-1.67, 1.76])
    univector<double> sig;

    // Puntos intermedios interpolados linealmente entre cada par de muestras de sig,
    // para sub-stepping del modelo sináptico dentro de cada paso de la señal
    univector<double> interpolated_points;

    // Número de sub-pasos del modelo por cada muestra de sig (factor de interpolación)
    size_t points_factor;

    // Paso temporal del modelo neuronal (ms) seleccionado para coincidir con el período de burst
    double dt;
};

// Lee, escala y sub-muestrea una señal de un CSV para hacerla compatible con el modelo neuronal.
// Devuelve nullopt si no se encuentra un dt adecuado para la frecuencia de burst.
std::optional<ScaledSigResult> scale_sig(
    const std::string &csv_path,
    size_t column_idx,
    double csv_step,
    double start_time,
    double use_time,
    double observation_time,
    NumericIntegrator integrator,
    NeuronModel model,
    bool check_drift);

#endif
