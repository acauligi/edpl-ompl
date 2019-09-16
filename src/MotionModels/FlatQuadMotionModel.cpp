/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2014, Texas A&M University
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Texas A&M University nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Authors: Abhishek Cauligi */ 

#include <cassert>
#include <tinyxml.h>

#include "Utils/FIRMUtils.h"
#include "Spaces/FlatQuadBeliefSpace.h"
#include "MotionModels/FlatQuadMotionModel.h"

//Produce the next state, given the current state, a control and a noise
void FlatQuadMotionModel::Evolve(const ompl::base::State *state, const ompl::control::Control *control, const NoiseType& w, ompl::base::State *result) {
  using namespace arma;

  typedef typename MotionModelMethod::StateType StateType;

  arma::colvec u = OMPL2ARMA(control);

  const colvec& Un = w.subvec(0, this->controlDim_-1);
  const colvec& Wg = w.subvec(this->controlDim_, this->noiseDim_-1);

  // TODO(acauligi): what is appropriate noise model for triple integrator?
  colvec x = state->as<StateType>()->getArmaData();
  // x = this->Ak_*x + this->Bk_*u + this->Gk_*(Un+Wg); 
  x = this->Ak_*x + this->Bk_*u + this->Gk_ * Un; 

  result->as<StateType>()->setArmaData(x);
}


void FlatQuadMotionModel::generateOpenLoopControls(const ompl::base::State *startState,
                                                  const ompl::base::State *endState,
                                                  std::vector<ompl::control::Control*> &openLoopControls) {
  // Jerk control input given by cubic polynomial over Tf horizon 
  using namespace arma;
  typedef typename MotionModelMethod::StateType StateType;

  double dh = this->dt_;
  double Tf = this->Tf_;
  double translation_steps = Tf / dh;

  colvec start = startState->as<StateType>()->getArmaData(); // turn into colvec (in radian)
  colvec target = endState->as<StateType>()->getArmaData(); // turn into colvec (in radian)

  // Math from "A computationally efficient motion primitive for quadrocopter
  //  trajectory generation" by Mueller et al. (2015)
  // TODO(acauligi): add bisection search for time-optimality and check for constraint satisfaction
  double Tf5 = pow(Tf,5);
  mat inv_matrix = zeros<mat>(3,3);   // Eq. 25 from paper
  inv_matrix.row(0) = rowvec({720, -360*Tf, 60*pow(Tf,2)});
  inv_matrix.row(1) = rowvec({-360*Tf, 168*pow(Tf,2), -24*pow(Tf,3)});
  inv_matrix.row(2) = rowvec({60*pow(Tf,2), -24*pow(Tf,3), 3*pow(Tf,4)});
  inv_matrix *= 1/Tf5;

  // difference in pose, vel, and acc for each coordinate of flat state
  // each column of deltas is (dp,dv,da) for (x,y,z,yaw)
  mat deltas = zeros<mat>(3,4);
  for (int ii=0; ii<4; ii++) {
    deltas.col(ii) = colvec({target(ii)-start(ii), 
                        target(4+ii)-start(4+ii), 
                        target(8+ii)-start(8+ii)});
  }

  // Recover (alpha,beta,gamma) coefficients for polynomial
  // Returns a 3x4 matrix with each col as (alpha,beta,gamma)
  // for the four flat states
  mat poly_coeffs = zeros<mat>(4,3);
  poly_coeffs = inv_matrix * deltas;

  // Sample continuous time control input at discrete interval
  colvec u_const = zeros<colvec>(this->controlDim_);
  for(int ii=0; ii<translation_steps; ii++) {
    double t = ii*dh/Tf;
    // Iterate through to get (jx,jy,jz,jyaw)
    for (int jj=0; jj<4; jj++) {
      double alpha = poly_coeffs(0,jj);
      double beta = poly_coeffs(1,jj);
      double gamma = poly_coeffs(2,jj);
      u_const(jj) = 0.5*alpha*pow(t,2) + beta*t+ gamma;
    }

    ompl::control::Control *tempControl = si_->allocControl();
    ARMA2OMPL(u_const, tempControl);
    openLoopControls.push_back(tempControl);
  }
}

void FlatQuadMotionModel::generateOpenLoopControlsForPath(const ompl::geometric::PathGeometric path, std::vector<ompl::control::Control*> &openLoopControls) {
  for(int i=0;i<path.getStateCount()-1;i++) {
    std::vector<ompl::control::Control*> olc;
    this->generateOpenLoopControls(path.getState(i),path.getState(i+1),olc) ;
    openLoopControls.insert(openLoopControls.end(),olc.begin(),olc.end());
  }
}


typename FlatQuadMotionModel::NoiseType
FlatQuadMotionModel::generateNoise(const ompl::base::State *state, const ompl::control::Control* control) {
  using namespace arma;

  NoiseType noise(this->noiseDim_);

  colvec indepUn = randn(this->controlDim_,1);
  
  mat P_Un = controlNoiseCovariance(control);
  
  colvec Un = indepUn % sqrt((P_Un.diag()));

  colvec Wg = sqrt(P_Wg_) * randn(this->stateDim_,1);
  
  noise = join_cols(Un, Wg);

  return noise;
}

typename FlatQuadMotionModel::JacobianType
FlatQuadMotionModel::getStateJacobian(const ompl::base::State *state, const ompl::control::Control* control, const NoiseType& w) {
  using namespace arma;
  return this->Ak_; 
}

typename FlatQuadMotionModel::JacobianType
FlatQuadMotionModel::getControlJacobian(const ompl::base::State *state, const ompl::control::Control* control, const NoiseType& w) {
  using namespace arma;
  typedef typename MotionModelMethod::StateType StateType;
  return this->Bk_; 
}


typename FlatQuadMotionModel::JacobianType
FlatQuadMotionModel::getNoiseJacobian(const ompl::base::State *state, const ompl::control::Control* control, const NoiseType& w) {
  using namespace arma;
  typedef typename MotionModelMethod::StateType StateType;
  return this->Gk_; 
}

arma::mat FlatQuadMotionModel::processNoiseCovariance(const ompl::base::State *state, const ompl::control::Control* control) {
  using namespace arma;

  mat P_Un = controlNoiseCovariance(control);

  mat Q_processNoise = zeros<mat>(P_Un.n_rows + P_Wg_.n_rows, P_Un.n_cols + P_Wg_.n_cols);

  Q_processNoise.submat(0, 0, P_Un.n_rows-1, P_Un.n_cols-1) = P_Un;
  Q_processNoise.submat(P_Un.n_rows, P_Un.n_cols,
                P_Un.n_rows + P_Wg_.n_rows -1,
                P_Un.n_cols + P_Wg_.n_cols -1) = P_Wg_;

  return Q_processNoise;
}


arma::mat FlatQuadMotionModel::controlNoiseCovariance(const ompl::control::Control* control) {
  using namespace arma;

  arma::colvec u = OMPL2ARMA(control);

  colvec uStd = eta_ % u + sigma_;

  mat P_Un = diagmat(square(uStd));

  return P_Un;
}

void FlatQuadMotionModel::constructAB() {
  // LTI triple integrator model

  using namespace arma;
  // Set up continuous time dynamics matrices 
  this->A_ = zeros(this->stateDim_, this->stateDim_);
  this->B_ = zeros(this->stateDim_, this->controlDim_);
  for (int ii=0; ii<this->controlDim_; ii++) {
    this->A_[4+ii,4+ii] = 1;
    this->A_[8+ii,8+ii] = 1;
    this->B_[8+ii,ii] = 1;
  }

  // Set up discrete time update matrices
  this->Ak_ = eye(this->stateDim_, this->stateDim_);
  this->Bk_ = zeros(this->stateDim_, this->controlDim_);
  this->Gk_ = sqrt(this->dt_)*eye(this->stateDim_, this->controlDim_);
  for (int ii=0; ii<this->controlDim_; ii++) {
    this->Ak_[ii,4+ii] = this->dt_;
    this->Ak_[ii,8+ii] = 0.5*this->dt_ * this->dt_;
    this->Ak_[4+ii,8+ii] = this->dt_;

    this->Bk_[ii,ii] = 1/6*this->dt_*this->dt_*this->dt_;
    this->Bk_[4+ii,ii] = 1/2*this->dt_*this->dt_;
    this->Bk_[8+ii,ii] = this->dt_;
  }
}

void FlatQuadMotionModel::loadParameters(const char *pathToSetupFile) {
  using namespace arma;

  TiXmlDocument doc(pathToSetupFile);
  bool loadOkay = doc.LoadFile();

  if ( !loadOkay ) {
    printf( "Could not load setup file in motion model. Error='%s'. Exiting.\n", doc.ErrorDesc() );
    exit( 1 );
  }

  TiXmlNode* node = 0;
  TiXmlElement* itemElement = 0;

  node = doc.FirstChild( "MotionModels" );
  assert( node );

  TiXmlNode* child = 0;

  child = node->FirstChild("FlatQuadMotionModel");

  assert( child );
  itemElement = child->ToElement();
  assert( itemElement );

  double attribute_val = 0;

  // Bias standard deviation of the motion noise
  itemElement->QueryDoubleAttribute("sigmaV", &attribute_val) ;
  this->sigma_ = attribute_val * arma::ones<colvec>(this->controlDim_); 
  
  // Proportional standard deviation of the motion noise 
  itemElement->QueryDoubleAttribute("etaV", &attribute_val) ;
  this->eta_ = attribute_val * arma::ones<colvec>(this->controlDim_); 

  // Covariance of state additive noise
  itemElement->QueryDoubleAttribute("wind_noise_pos", &attribute_val) ;
  rowvec Wg_root_vec(this->stateDim_);
  Wg_root_vec = attribute_val * arma::ones<rowvec>(this->stateDim_);
  this->P_Wg_ = diagmat(square(Wg_root_vec));
  
  // MAV constraints in flat space
  itemElement->QueryDoubleAttribute("max_linear_velocity", &attribute_val) ;
  this->max_linear_velocity_      = attribute_val;
  itemElement->QueryDoubleAttribute("max_linear_acceleration", &attribute_val) ;
  this->max_linear_acceleration_  = attribute_val;
  itemElement->QueryDoubleAttribute("max_linear_jerk", &attribute_val) ;
  this->max_linear_jerk_          = attribute_val;
  itemElement->QueryDoubleAttribute("dt", &attribute_val) ;
  this->dt_                     = attribute_val;
  itemElement->QueryDoubleAttribute("Tf", &attribute_val) ;
  this->Tf_ = attribute_val;

  OMPL_INFORM("FlatQuadMotionModel: sigma_ = ");
  std::cout<<this->sigma_<<std::endl;

  OMPL_INFORM("FlatQuadMotionModel: eta_ = ");
  std::cout<<this->eta_<<std::endl;

  OMPL_INFORM("FlatQuadMotionModel: P_Wg_ = ");
  std::cout<<this->P_Wg_<<std::endl;

  OMPL_INFORM("FlatQuadMotionModel: max Linear Velocity (m/s)    = %f",       this->max_linear_velocity_);
  OMPL_INFORM("FlatQuadMotionModel: max Linear Acceleration (m/s^2)    = %f", this->max_linear_acceleration_);
  OMPL_INFORM("FlatQuadMotionModel: max Linear Jerk (m/s^3)    = %f",         this->max_linear_jerk_);

  OMPL_INFORM("FlatQuadMotionModel: Timestep (seconds) = %f", this->dt_);
  OMPL_INFORM("FlatQuadMotionModel: Maximum final time (seconds) = %f", this->Tf_);
}
