/*
 * Developed by Rodrigo Amaducci (rodrigo.amaducci@uam.es) from the Plugin Template provided by RTXI
 * Grupo de Neurocomputación Biológica (GNB), Universidad Autónoma de Madrid, 2018
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

#include <default_gui_model.h>

class TsodyksMarkramSynapse : public DefaultGUIModel
{

  Q_OBJECT

public:
  TsodyksMarkramSynapse(void);
  virtual ~TsodyksMarkramSynapse(void);

  void execute(void);
  void createGUI(DefaultGUIModel::variable_t*, int);
  void customizeGUI(void);

protected:
  virtual void update(DefaultGUIModel::update_flags_t);

private:
  double period;


  double vars_syn[3];
  double params_syn[5];


  void runge_kutta_65 (void (*f) (double *, double *, double *, double), int dim, double dt, double * vars, double * params, double aux);
  static void function (double * vars, double * ret, double * params, double vpre);


  void initParameters();

private slots:
  // these are custom functions that can also be connected to events
  // through the Qt API. they must be implemented in plugin_template.cpp

  void aBttn_event(void);
  void bBttn_event(void);
};
