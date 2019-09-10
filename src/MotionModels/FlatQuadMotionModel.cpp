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


#include <tinyxml.h>
#include "Spaces/FlatQuadBeliefSpace.h"
#include "MotionModels/FlatQuadMotionModel.h"
#include "Utils/FIRMUtils.h"

#include<cassert>

//Produce the next state, given the current state, a control and a noise
void FlatQuadMotionModel::Evolve(const ompl::base::State *state, const ompl::control::Control *control, const NoiseType& w, ompl::base::State *result) {
  using namespace arma;

  typedef typename MotionModelMethod::StateType StateType;

  arma::colvec u = OMPL2ARMA(control);

  const colvec& Un = w.subvec(0, this->controlDim_-1);
  const colvec& Wg = w.subvec(this->controlDim_, this->noiseDim_-1);

  // TODO(acauligi): what is appropriate noise model for triple integrator?
  colvec x = state->as<StateType>()->getArmaData();
  x = this->Ak_*x + this->Bk_*u + this->Gk_*(Un+Wg); 

  result->as<StateType>()->setArmData(x);
}


void FlatQuadMotionModel::generateOpenLoopControls(const ompl::base::State *startState,
                                                  const ompl::base::State *endState,
                                                  std::vector<ompl::control::Control*> &openLoopControls) {
  using namespace arma;
  typedef typename MotionModelMethod::StateType StateType;

  colvec start = startState->as<StateType>()->getArmaData(); // turn into colvec (in radian)
  colvec target = endState->as<StateType>()->getArmaData(); // turn into colvec (in radian)

  double delta_disp = 0;

  double translation_steps = 0;

  translation_steps = floor( std::max( fabs((target[0]-start[0])/(maxLinearVelocity_*this->dt_)), fabs((target[1]-start[1])/(maxLinearVelocity_*this->dt_))) );

  colvec u_const;
  u_const << (target[0]-start[0]) / (translation_steps*this->dt_) << endr
          << (target[1]-start[1]) / (translation_steps*this->dt_) << endr;

  for(int i=0; i<translation_steps; i++) {
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

  colvec xData = state->as<StateType>()->getArmaData();

  assert (xData.n_rows == (size_t)this->stateDim_);

  mat G(this->stateDim_,this->noiseDim_);

  G   << 1 << 0 << 1 << 0 << endr
      << 0 << 1 << 0 << 1 << endr;

  G *= sqrt(this->dt_);
  return G;
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

  double sigmaV=0;
  double etaV = 0;
  double windNoisePos=0;
  double maxLinearVelocity=0;
  double maxLinearAcceleration=0;
  double maxLinearJerk=0;
  double dt = 0;

  itemElement->QueryDoubleAttribute("sigmaV", &sigmaV) ;
  itemElement->QueryDoubleAttribute("etaV", &etaV) ;
  itemElement->QueryDoubleAttribute("wind_noise_pos", &windNoisePos) ;
  itemElement->QueryDoubleAttribute("max_linear_velocity", &maxLinearVelocity) ;
  itemElement->QueryDoubleAttribute("max_linear_acceleration", &maxLinearAcceleration) ;
  itemElement->QueryDoubleAttribute("max_linear_jerk", &maxLinearJerk) ;
  itemElement->QueryDoubleAttribute("dt", &dt) ;

  // parameters used to compute process noise values
  this->sigma_ << sigmaV << sigmaV << sigmaV << sigmaV << endr;
  this->eta_  << etav << etav << etav << etav << endr;

  rowvec Wg_root_vec(4);
  Wg_root_vec << windNoisePos << windNoisePos << windNoisePos << windNoisePos << endr;
  this->P_Wg_ = diagmat(square(Wg_root_vec));

  this->maxLinearVelocity_      = maxLinearVelocity;
  this->maxLinearAcceleration_  = maxLinearAcceleration;
  this->maxLinearJerk_          = maxLinearJerk;
  this->dt_                     = dt;

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
  this->Gk_ = sqrt(this->dt_)*eye(this->stateDim_, this->stateDim_);
  for (int ii=0; ii<this->controlDim_; ii++) {
    this->Ak_[ii,4+ii] = this->dt_;
    this->Ak_[ii,8+ii] = 0.5*this->dt_ * this->dt_;
    this->Ak_[4+ii,8+ii] = this->dt_;

    this->Bk_[ii,ii] = 1/6*this->dt_*this->dt_*this->dt_;
    this->Bk_[4+ii,ii] = 1/2*this->dt_*this->dt_;
    this->Bk_[8+ii,ii] = this->dt_;
  }

  OMPL_INFORM("FlatQuadMotionModel: sigma_ = ");
  std::cout<<this->sigma_<<std::endl;

  OMPL_INFORM("FlatQuadMotionModel: eta_ = ");
  std::cout<<this->eta_<<std::endl;

  OMPL_INFORM("FlatQuadMotionModel: P_Wg_ = ");
  std::cout<<this->P_Wg_<<std::endl;

  OMPL_INFORM("FlatQuadMotionModel: max Linear Velocity (m/s)    = %f", maxLinearVelocity);
  OMPL_INFORM("FlatQuadMotionModel: max Linear Acceleration (m/s^2)    = %f", maxLinearAcceleration);
  OMPL_INFORM("FlatQuadMotionModel: max Linear Jerk (m/s^3)    = %f", maxLinearJerk);

  OMPL_INFORM("FlatQuadMotionModel: Timestep (seconds) = %f", dt);
}
