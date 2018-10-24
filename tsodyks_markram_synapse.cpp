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

#include "tsodyks_markram_synapse.h"
#include <iostream>
#include <main_window.h>

extern "C" Plugin::Object*
createRTXIPlugin(void)
{
  return new TsodyksMarkramSynapse();
}


#define X 0
#define Y 1
#define U 2

#define TAU_REC 0
#define TAU_INACT 1
#define TAU_REL 2
#define RELEASE_PROB 3
#define A 4

static DefaultGUIModel::variable_t vars[] = {
  { "Vpre (V)", "Presynaptic neuron membrane potential", DefaultGUIModel::INPUT,},
  { "I (nA)", "Synaptic current", DefaultGUIModel::OUTPUT,},
  { "Tau recovery", "Tau recovery parameter", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  { "Tau inactivation", "Tau inactivation parameter", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  { "Tau release", "Tau release parameter", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  { "A", "Absolute synaptic efficacy", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  { "Release probability", "Release probability", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

TsodyksMarkramSynapse::TsodyksMarkramSynapse(void)
  : DefaultGUIModel("Tsodyks Markram Synapse", ::vars, ::num_vars)
{
  setWhatsThis("<p><b>Tsodyks Markram Synapse:</b><br>synapse model from (Tsodyks and Markram, 1998).</p>");
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

TsodyksMarkramSynapse::~TsodyksMarkramSynapse(void)
{
}


void TsodyksMarkramSynapse::runge_kutta_65 (void (*f) (double *, double *, double *, double), int dim, double dt, double * vars, double * params, double aux) {
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


void TsodyksMarkramSynapse::function(double * vars, double * ret, double * params, double vpre) {
    ret[X] = (1 - vars[X] - vars[Y]) / params[TAU_REC];
    ret[Y] = (-vars[Y]) / params[TAU_INACT];
    ret[U] = (params[RELEASE_PROB] - vars[U]) / params[TAU_REL];

    if (vpre > 0) {
        ret[X] += -vars[U] * vars[X];
        ret[Y] += vars[U] * vars[X];
        ret[U] += params[RELEASE_PROB] * (1.0 - vars[U]);
    }
}

void
TsodyksMarkramSynapse::execute(void)
{
    double vpre = input(0);

    runge_kutta_65(&function, 3, 0.001, vars_syn, params_syn, vpre);

    output(0) = params_syn[A] * vars_syn[Y];

    return;
}

void
TsodyksMarkramSynapse::initParameters(void)
{
  vars_syn[X] = 0.0;
  vars_syn[Y] = 0.0;
  vars_syn[U] = 0.0;

  params_syn[TAU_REC] = 800.0;
  params_syn[TAU_INACT] = 3.0;
  params_syn[TAU_REL] = 50.0;
  params_syn[RELEASE_PROB] = 0.67;
  params_syn[A] = 250.0;
}

void
TsodyksMarkramSynapse::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag) {
    case INIT:
      period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
      setParameter("Tau recovery", params_syn[TAU_REC]);
      setParameter("Tau inactivation", params_syn[TAU_INACT]);
      setParameter("Tau release", params_syn[TAU_REL]);
      setParameter("A", params_syn[A]);
      setParameter("Release probability", params_syn[RELEASE_PROB]);

      break;

    case MODIFY:
      params_syn[TAU_REC] = getParameter("Tau recovery").toDouble();
      params_syn[TAU_INACT] = getParameter("Tau inactivation").toDouble();
      params_syn[TAU_REL] = getParameter("Tau release").toDouble();
      params_syn[A] = getParameter("A").toDouble();
      params_syn[RELEASE_PROB] = getParameter("Release probability").toDouble();

      break;

    case UNPAUSE:
      break;

    case PAUSE:
      break;

    case PERIOD:
      period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
      break;

    default:
      break;
  }
}

void
TsodyksMarkramSynapse::customizeGUI(void)
{
}

// functions designated as Qt slots are implemented as regular C++ functions
void
TsodyksMarkramSynapse::aBttn_event(void)
{
}

void
TsodyksMarkramSynapse::bBttn_event(void)
{
}
