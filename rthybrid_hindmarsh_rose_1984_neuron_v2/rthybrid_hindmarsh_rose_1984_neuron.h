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
 * This is a template header file for a user modules derived from
 * DefaultGUIModel with a custom GUI.
 */

#include <default_gui_model.h>

class RTHybridHindmarshRose1984Neuron : public DefaultGUIModel
{

  Q_OBJECT

public:
  RTHybridHindmarshRose1984Neuron(void);
  virtual ~RTHybridHindmarshRose1984Neuron(void);

  void execute(void);
  void createGUI(DefaultGUIModel::variable_t*, int);
  void customizeGUI(void);

protected:
  virtual void update(DefaultGUIModel::update_flags_t);

private:
  double vars_model[3];
  double params_model[11];
  double period, freq;
  double burst_duration, burst_duration_value;
  double s_points;

  void initParameters();

  void runge_kutta_65 (void (*f) (double *, double *, double *, double), int dim, double dt, double * vars, double * params, double aux);
  void select_dt_neuron_model (double * dts, double * pts, unsigned int length, double pts_live, double * dt, double * pts_burst);
  double set_pts_burst (double sec_per_burst);
  static double nm_hindmarsh_rose_1986_v (double * vars, double * params);
  static double nm_hindmarsh_rose_1986_y (double * vars, double * params);
  static double nm_hindmarsh_rose_1986_z (double * vars, double * params);
  static void nm_hindmarsh_rose_1986_f (double * vars, double * ret, double * params, double syn);

private slots:
  // these are custom functions that can also be connected to events
  // through the Qt API. they must be implemented in plugin_template.cpp

  void aBttn_event(void);
  void bBttn_event(void);
};
