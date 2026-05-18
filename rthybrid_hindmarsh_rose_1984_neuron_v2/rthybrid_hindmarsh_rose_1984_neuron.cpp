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

/*
 * This is a template implementation file for a user module derived from
 * DefaultGUIModel with a custom GUI.
 */

#include "rthybrid_hindmarsh_rose_1984_neuron.h"
#include <iostream>
#include <main_window.h>

// Variables
#define NM_HINDMARSH_ROSE_1986_V 0
#define NM_HINDMARSH_ROSE_1986_Y 1
#define NM_HINDMARSH_ROSE_1986_Z 2

// Parameters
#define NM_HINDMARSH_ROSE_1986_DT 0
#define NM_HINDMARSH_ROSE_1986_B 1
#define NM_HINDMARSH_ROSE_1986_A 2
#define NM_HINDMARSH_ROSE_1986_I 3
#define NM_HINDMARSH_ROSE_1986_SYN 4
#define NM_HINDMARSH_ROSE_1986_C 5
#define NM_HINDMARSH_ROSE_1986_D 6
#define NM_HINDMARSH_ROSE_1986_R 7
#define NM_HINDMARSH_ROSE_1986_S 8
#define NM_HINDMARSH_ROSE_1986_XR 9

#define NM_HINDMARSH_ROSE_1986_VH 10

extern "C" Plugin::Object*
createRTXIPlugin(void)
{
  return new RTHybridHindmarshRose1984Neuron();
}

static DefaultGUIModel::variable_t vars[] = {
    {"Burst duration (s)", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    {"I", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

    { "a", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    { "b", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    { "c", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    { "d", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    { "r", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    { "s", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    { "xr", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

    { "Vh", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

    { "x0", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    { "y0", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    { "z0", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

    { "Vm (v)", "Membrane potential (in V)", DefaultGUIModel::OUTPUT, },
    { "Vm (mV)", "Membrane potential (in mV)", DefaultGUIModel::OUTPUT, },

    { "Isyn (nA)", "Synaptic input current (in nA)", DefaultGUIModel::INPUT,},
    { "Burst duration (s)", "", DefaultGUIModel::INPUT,},

    { "v", "", DefaultGUIModel::STATE,},
    { "s_points", "", DefaultGUIModel::STATE,},
    { "dt", "", DefaultGUIModel::STATE,},
    { "syn", "", DefaultGUIModel::STATE,},
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

RTHybridHindmarshRose1984Neuron::RTHybridHindmarshRose1984Neuron(void)
  : DefaultGUIModel("RTHybrid Hindmarsh-Rose (1984) neuron model", ::vars, ::num_vars)
{
  setWhatsThis("<p><b>RTHybrid Hindmarsh-Rose (1984) neuron model</b></p>");
  DefaultGUIModel::createGUI(vars,
                             num_vars); // this is required to create the GUI
  //customizeGUI();
  initParameters();
  update(INIT); // this is optional, you may place initialization code directly
                // into the constructor
  refresh();    // this is required to update the GUI with parameter and state
                // values
  QTimer::singleShot(0, this, SLOT(resizeMe()));
}

RTHybridHindmarshRose1984Neuron::~RTHybridHindmarshRose1984Neuron(void)
{
}

void RTHybridHindmarshRose1984Neuron::runge_kutta_65 (void (*f) (double *, double *, double *, double), int dim, double dt, double * vars, double * params, double aux) {
    double apoyo[dim], retorno[dim];
    double k[6][dim];
    int j;

    (*f)(vars, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[0][j] = dt * retorno[j];
        apoyo[j] = vars[j] + k[0][j] * 0.2;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[1][j] = dt * retorno[j];
        apoyo[j] = vars[j] + k[0][j] * 0.075 + k[1][j] * 0.225;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[2][j] = dt * retorno[j];
        apoyo[j] = vars[j] + k[0][j] * 0.3 - k[1][j] * 0.9 + k[2][j] * 1.2;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[3][j] = dt * retorno[j];
        apoyo[j] = vars[j] + k[0][j] * 0.075 + k[1][j] * 0.675 - k[2][j] * 0.6 + k[3][j] * 0.75;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[4][j] = dt * retorno[j];
        apoyo[j] = vars[j]
                + k[0][j] * 0.660493827160493
                + k[1][j] * 2.5
                - k[2][j] * 5.185185185185185
                + k[3][j] * 3.888888888888889
                - k[4][j] * 0.864197530864197;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[5][j] = dt * retorno[j];
        apoyo[j] = vars[j]
                + k[0][j]*0.1049382716049382
                + k[2][j]*0.3703703703703703
                + k[3][j]*0.2777777777777777
                + k[4][j]*0.2469135802469135;
    }


    for(j = 0; j < dim; ++j) {
        vars[j] += k[0][j]*0.098765432098765+
                   k[2][j]*0.396825396825396+
                   k[3][j]*0.231481481481481+
                   k[4][j]*0.308641975308641-
                   k[5][j]*0.035714285714285;
    }

    return;
}


void RTHybridHindmarshRose1984Neuron::select_dt_neuron_model (double * dts, double * pts, unsigned int length, double pts_live, double * dt, double * pts_burst) {
    double aux = pts_live;
    double factor = 1;
    double intpart, fractpart;
    int flag = 0;
    int i;

    *dt = -1;
    *pts_burst = -1;

    while (aux < pts[0]) {
        aux = pts_live * factor;
        factor += 1;

        for (i = length - 1; i >= 0; i--) {
          if (pts[i] > aux) {
            *dt = dts[i];
              *pts_burst = pts[i];

              fractpart = modf(*pts_burst / pts_live, &intpart);

              if (fractpart <= 0.1*intpart) flag = 1;

              break;
          }
            
        }

        if (flag == 1) break;
    }

    if (flag == 0) {
        for (i = length - 1; i >= 0; i--) {
          if (pts[i] > aux) {
            *dt = dts[i];
              *pts_burst = pts[i];

              break;  
          }
        }
    }

    return;
}


double RTHybridHindmarshRose1984Neuron::set_pts_burst (double sec_per_burst) {
	int length = 0;
	int method = 3;
	double pts_match = sec_per_burst * freq;
	double pts_burst, dt;

	length =  144;
	double dts[] = {0.000500, 0.000600, 0.000700, 0.000800, 0.000900, 0.001000, 0.001100, 0.001200, 0.001300, 0.001400, 0.001500, 0.001600, 0.001800, 0.002000, 0.002200, 0.002500, 0.002800, 0.002900, 0.003000, 0.003100, 0.003200, 0.003300, 0.003400, 0.003500, 0.003600, 0.003700, 0.003800, 0.003900, 0.004000, 0.004100, 0.004200, 0.004300, 0.004400, 0.004500, 0.004600, 0.004700, 0.004800, 0.004900, 0.005000, 0.005100, 0.005200, 0.005400, 0.005600, 0.005800, 0.006000, 0.006200, 0.006400, 0.006600, 0.006800, 0.007000, 0.007200, 0.007400, 0.007700, 0.008000, 0.008300, 0.008600, 0.008900, 0.009200, 0.009600, 0.010000, 0.010400, 0.010900, 0.011400, 0.011900, 0.012500, 0.013100, 0.013800, 0.014600, 0.015400, 0.016300, 0.017300, 0.018500, 0.019900, 0.021500, 0.023300, 0.025500, 0.028100, 0.028400, 0.028700, 0.029000, 0.029400, 0.029800, 0.030200, 0.030600, 0.031000, 0.031400, 0.031800, 0.032200, 0.032600, 0.033000, 0.033400, 0.033900, 0.034400, 0.034900, 0.035400, 0.035900, 0.036400, 0.036900, 0.037400, 0.038000, 0.038600, 0.039200, 0.039800, 0.040400, 0.041000, 0.041700, 0.042400, 0.043100, 0.043800, 0.044500, 0.045300, 0.046100, 0.046900, 0.047700, 0.048600, 0.049500, 0.050400, 0.051400, 0.052400, 0.053400, 0.054500, 0.055600, 0.056800, 0.058000, 0.059300, 0.060600, 0.062000, 0.063400, 0.064900, 0.066500, 0.068200, 0.069900, 0.071700, 0.073600, 0.075600, 0.077700, 0.079900, 0.082300, 0.084800, 0.087500, 0.090300, 0.093300, 0.096500, 0.100000};
	double pts[] = {577638.000000, 481366.000000, 412599.000000, 357615.500000, 317880.000000, 286092.500000, 259143.333333, 237548.000000, 218869.500000, 203236.000000, 189687.000000, 177634.000000, 157897.000000, 142001.833333, 129024.142857, 113496.125000, 101304.555556, 97811.222222, 94527.400000, 91478.200000, 88619.400000, 85916.636364, 83389.636364, 81007.090909, 78743.583333, 76615.416667, 74599.250000, 72676.000000, 70859.076923, 69130.846154, 67476.642857, 65907.357143, 64402.666667, 62971.466667, 61602.533333, 60286.187500, 59030.250000, 57825.562500, 56664.411765, 55553.294118, 54485.000000, 52463.222222, 50586.263158, 48841.842105, 47211.050000, 45685.666667, 44255.818182, 42914.772727, 41650.739130, 40459.083333, 39335.208333, 38270.680000, 36778.346154, 35398.000000, 34117.571429, 32926.517241, 31815.833333, 30777.612903, 29493.939394, 28313.588235, 27223.638889, 25974.405405, 24834.410256, 23790.268293, 22647.767442, 21609.977778, 20513.166667, 19388.627451, 18381.132075, 17365.719298, 16361.600000, 15299.937500, 14223.202899, 13164.400000, 12147.123457, 11098.876404, 10071.693878, 9965.282828, 9861.100000, 9759.059406, 9626.242718, 9497.009615, 9371.179245, 9248.672897, 9129.293578, 9012.981818, 8899.594595, 8789.000000, 8681.149123, 8575.896552, 8473.170940, 8348.176471, 8226.809917, 8108.934426, 7994.379032, 7883.015873, 7774.710938, 7669.348837, 7566.801527, 7447.308271, 7331.525926, 7219.291971, 7110.435714, 7004.816901, 6902.298611, 6786.417808, 6674.355705, 6565.940397, 6460.987013, 6359.339744, 6247.012579, 6138.592593, 6033.866667, 5932.660714, 5822.777778, 5716.896552, 5614.796610, 5505.541436, 5400.467391, 5299.324468, 5192.348958, 5089.615385, 4982.075000, 4878.985294, 4772.014354, 4669.633803, 4564.178899, 4463.381166, 4360.214912, 4255.294872, 4149.216667, 4048.296748, 3946.654762, 3844.764479, 3743.041353, 3641.872263, 3541.583630, 3438.300000, 3336.926421, 3233.951299, 3133.666667, 3032.899696, 2932.320588, 2829.684659};
	select_dt_neuron_model(dts, pts, length, pts_match, &(params_model[NM_HINDMARSH_ROSE_1986_DT]), &(pts_burst));

	return pts_burst;
}


double RTHybridHindmarshRose1984Neuron::nm_hindmarsh_rose_1986_v (double * vars, double * params) {
	return vars[NM_HINDMARSH_ROSE_1986_Y] + params[NM_HINDMARSH_ROSE_1986_B] * (vars[NM_HINDMARSH_ROSE_1986_V]*vars[NM_HINDMARSH_ROSE_1986_V]) - params[NM_HINDMARSH_ROSE_1986_A] * (vars[NM_HINDMARSH_ROSE_1986_V]*vars[NM_HINDMARSH_ROSE_1986_V]*vars[NM_HINDMARSH_ROSE_1986_V]) - vars[NM_HINDMARSH_ROSE_1986_Z] + params[NM_HINDMARSH_ROSE_1986_I] - params[NM_HINDMARSH_ROSE_1986_SYN];
}

double RTHybridHindmarshRose1984Neuron::nm_hindmarsh_rose_1986_y (double * vars, double * params) {
	return params[NM_HINDMARSH_ROSE_1986_C] - params[NM_HINDMARSH_ROSE_1986_D] * vars[NM_HINDMARSH_ROSE_1986_V]*vars[NM_HINDMARSH_ROSE_1986_V] - vars[NM_HINDMARSH_ROSE_1986_Y];
}

double RTHybridHindmarshRose1984Neuron::nm_hindmarsh_rose_1986_z (double * vars, double * params) {
	return params[NM_HINDMARSH_ROSE_1986_R] * (params[NM_HINDMARSH_ROSE_1986_S] * (vars[NM_HINDMARSH_ROSE_1986_V] - params[NM_HINDMARSH_ROSE_1986_XR]) - (vars[NM_HINDMARSH_ROSE_1986_Z]*params[NM_HINDMARSH_ROSE_1986_VH]) );
}

/**
 * @brief Hindmarsh_Rose_1986 neuron model differential equations.
 * @param[in] vars Neuron model variables
 * @param[out] ret Return values array
 * @param[in] params Neuron models parameters
 * @param[in] syn Synapse input current value
 */

void RTHybridHindmarshRose1984Neuron::nm_hindmarsh_rose_1986_f (double * vars, double * ret, double * params, double syn) {
	params[NM_HINDMARSH_ROSE_1986_SYN] = syn;

	ret[NM_HINDMARSH_ROSE_1986_V] = nm_hindmarsh_rose_1986_v(vars, params);
	ret[NM_HINDMARSH_ROSE_1986_Y] = nm_hindmarsh_rose_1986_y(vars, params);
	ret[NM_HINDMARSH_ROSE_1986_Z] = nm_hindmarsh_rose_1986_z(vars, params);
}

void
RTHybridHindmarshRose1984Neuron::execute(void)
{
	int i;

    if (burst_duration_value <= -1) {
        burst_duration = input(1);
        s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq));
        if (s_points < 1) s_points = 1;
    }
    

	for (i = 0; i < s_points; i++) {
	  runge_kutta_65(nm_hindmarsh_rose_1986_f, 3, params_model[NM_HINDMARSH_ROSE_1986_DT], vars_model, params_model, input(0));
	}

	output(0) = vars_model[NM_HINDMARSH_ROSE_1986_V] / 1000.0;
	output(1) = vars_model[NM_HINDMARSH_ROSE_1986_V];

	return;
}

void
RTHybridHindmarshRose1984Neuron::initParameters(void)
{
  period = RT::System::getInstance()->getPeriod() * 1e-6;

	burst_duration_value = 1.0;
    burst_duration = burst_duration_value;
    freq = 1.0 / (period * 1e-3);
    s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq)); 
    if (s_points == 0) s_points = 1;

    vars_model[NM_HINDMARSH_ROSE_1986_V] = -0.712841;
	vars_model[NM_HINDMARSH_ROSE_1986_Y] = -1.936880;
	vars_model[NM_HINDMARSH_ROSE_1986_Z] = 3.165680;

	params_model[NM_HINDMARSH_ROSE_1986_I] = 3.0;

	params_model[NM_HINDMARSH_ROSE_1986_R] = 0.0021;
	params_model[NM_HINDMARSH_ROSE_1986_S] = 4.0;
	params_model[NM_HINDMARSH_ROSE_1986_XR] = -1.6;

	params_model[NM_HINDMARSH_ROSE_1986_VH] = 0.2;

	params_model[NM_HINDMARSH_ROSE_1986_A] = 1.0;
    params_model[NM_HINDMARSH_ROSE_1986_B] = 3.0;
    params_model[NM_HINDMARSH_ROSE_1986_C] = 1.0;
    params_model[NM_HINDMARSH_ROSE_1986_D] = 5.0;
}

void
RTHybridHindmarshRose1984Neuron::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag) {
    case INIT:
  		period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
  		freq = 1.0 / (period * 1e-3);

  		setParameter("Burst duration (s)", burst_duration_value);
  		setParameter("I", params_model[NM_HINDMARSH_ROSE_1986_I]);

  		setParameter("a", params_model[NM_HINDMARSH_ROSE_1986_A]);
  		setParameter("b", params_model[NM_HINDMARSH_ROSE_1986_B]);
  		setParameter("c", params_model[NM_HINDMARSH_ROSE_1986_C]);
  		setParameter("d", params_model[NM_HINDMARSH_ROSE_1986_D]);

  		setParameter("r", params_model[NM_HINDMARSH_ROSE_1986_R]);
  		setParameter("s", params_model[NM_HINDMARSH_ROSE_1986_S]);
  		setParameter("xr", params_model[NM_HINDMARSH_ROSE_1986_XR]);

  		setParameter("Vh", params_model[NM_HINDMARSH_ROSE_1986_VH]);

  		setParameter("x0", vars_model[NM_HINDMARSH_ROSE_1986_V]);
  		setParameter("y0", vars_model[NM_HINDMARSH_ROSE_1986_Y]);
  		setParameter("z0", vars_model[NM_HINDMARSH_ROSE_1986_Z]);

  		setState("v", vars_model[NM_HINDMARSH_ROSE_1986_V]);
  		setState("s_points", s_points);
  		setState("dt", params_model[NM_HINDMARSH_ROSE_1986_DT]);
  		setState("syn", params_model[NM_HINDMARSH_ROSE_1986_SYN]);

      break;

    case MODIFY:
        burst_duration_value = getParameter("Burst duration (s)").toDouble();
        freq = 1.0 / (period * 1e-3);

        if (burst_duration_value <= -1) {
            burst_duration = input(1);
        } else {
            burst_duration = burst_duration_value;
        }

        s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq));
        if (s_points < 1) s_points = 1;

  		params_model[NM_HINDMARSH_ROSE_1986_A] = getParameter("a").toDouble();
  		params_model[NM_HINDMARSH_ROSE_1986_B] = getParameter("b").toDouble();
  		params_model[NM_HINDMARSH_ROSE_1986_C] = getParameter("c").toDouble();
  		params_model[NM_HINDMARSH_ROSE_1986_D] = getParameter("d").toDouble();

  		params_model[NM_HINDMARSH_ROSE_1986_R] = getParameter("r").toDouble();
  		params_model[NM_HINDMARSH_ROSE_1986_S] = getParameter("s").toDouble();
  		params_model[NM_HINDMARSH_ROSE_1986_XR] = getParameter("xr").toDouble();

  		params_model[NM_HINDMARSH_ROSE_1986_VH] = getParameter("Vh").toDouble();

  		params_model[NM_HINDMARSH_ROSE_1986_I] = getParameter("I").toDouble();

  		vars_model[NM_HINDMARSH_ROSE_1986_V] = getParameter("x0").toDouble();
  		vars_model[NM_HINDMARSH_ROSE_1986_Y] = getParameter("y0").toDouble();
  		vars_model[NM_HINDMARSH_ROSE_1986_Z] = getParameter("z0").toDouble();
      break;

    case UNPAUSE:
  		freq = 1.0 / (period * 1e-3);
  		s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq)); 
  		if (s_points == 0) s_points = 1;

      break;

    case PAUSE:
      output(0) = 0;
      output(1) = 0;
      break;

    case PERIOD:
		period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
		freq = 1.0 / (period * 1e-3);
		s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq)); 
		if (s_points == 0) s_points = 1;

      break;

    default:
      break;
  }
}

void
RTHybridHindmarshRose1984Neuron::customizeGUI(void)
{
  QGridLayout* customlayout = DefaultGUIModel::getLayout();

  QGroupBox* button_group = new QGroupBox;

  QPushButton* abutton = new QPushButton("Button A");
  QPushButton* bbutton = new QPushButton("Button B");
  QHBoxLayout* button_layout = new QHBoxLayout;
  button_group->setLayout(button_layout);
  button_layout->addWidget(abutton);
  button_layout->addWidget(bbutton);
  QObject::connect(abutton, SIGNAL(clicked()), this, SLOT(aBttn_event()));
  QObject::connect(bbutton, SIGNAL(clicked()), this, SLOT(bBttn_event()));

  customlayout->addWidget(button_group, 0, 0);
  setLayout(customlayout);
}

// functions designated as Qt slots are implemented as regular C++ functions
void
RTHybridHindmarshRose1984Neuron::aBttn_event(void)
{
}

void
RTHybridHindmarshRose1984Neuron::bBttn_event(void)
{
}
