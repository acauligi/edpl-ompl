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

#include "Spaces/FlatQuadBeliefSpace.h"

double FlatQuadBeliefSpace::StateType::meanNormWeight_  = -1;
double FlatQuadBeliefSpace::StateType::covNormWeight_   = -1;
double FlatQuadBeliefSpace::StateType::reachDist_   = -1;
arma::colvec FlatQuadBeliefSpace::StateType::normWeights_ = arma::zeros<arma::colvec>(12);

bool FlatQuadBeliefSpace::StateType::isReached(ompl::base::State *state, bool relaxedConstraint) const {
  // subtract the two beliefs and get the norm
  arma::colvec stateDiff = this->getArmaData() - state->as<FlatQuadBeliefSpace::StateType>()->getArmaData();

  if(stateDiff[2] > boost::math::constants::pi<double>() ) {
    stateDiff[2] =  (stateDiff[2] - 2*boost::math::constants::pi<double>()) ;
  }
  if( stateDiff[2] < -boost::math::constants::pi<double>() ) {
    stateDiff[2] =  stateDiff[2] + 2*boost::math::constants::pi<double>() ;
  }

  arma::mat covDiff = this->getCovariance() -  state->as<FlatQuadBeliefSpace::StateType>()->getCovariance();

  arma::colvec covDiffDiag = covDiff.diag();

  // Need weighted supNorm of difference in means
  double meanNorm = arma::norm(stateDiff % normWeights_, "inf");

  double covDiagNorm = arma::norm(sqrt(abs(covDiffDiag)) % normWeights_, "inf");

  double norm2 =  std::max(meanNorm*meanNormWeight_, covDiagNorm*covNormWeight_) ;

  double reachConstraint  = reachDist_;

  if(relaxedConstraint) {
    reachConstraint *= 4;
  }

  if(norm2 <= reachConstraint) {
    return true;
  }

  return false;

}


ompl::base::State* FlatQuadBeliefSpace::allocState(void) const {
  StateType *state = new StateType();
  
  allocStateComponents(state);
  
  state->setYaw(0.0);
  
  return state;
}

void FlatQuadBeliefSpace::copyState(State *destination, const State *source) const {
  destination->as<StateType>()->setX(source->as<StateType>()->getX());
  destination->as<StateType>()->setY(source->as<StateType>()->getY());
  destination->as<StateType>()->setZ(source->as<StateType>()->getZ());
  destination->as<StateType>()->setYaw(source->as<StateType>()->getYaw());
  destination->as<StateType>()->setVelocity(source->as<StateType>()->getVelocity());
  destination->as<StateType>()->setAcceleration(source->as<StateType>()->getAcceleration());
  destination->as<StateType>()->setCovariance(source->as<StateType>()->getCovariance());
}

void FlatQuadBeliefSpace::freeState(State *state) const {
  CompoundStateSpace::freeState(state);
}

double FlatQuadBeliefSpace::distance(const State* state1, const State *state2) {
  double dx = state1->as<StateType>()->getX() - state2->as<StateType>()->getX();
  double dy = state1->as<StateType>()->getY() - state2->as<StateType>()->getY();
  double dy = state1->as<StateType>()->getZ() - state2->as<StateType>()->getZ();

  // TODO(acauligi): include rotation distance metric?
  return pow(dx*dx+dy*dy+dz*dz, 0.5);
}

void FlatQuadBeliefSpace::getRelativeState(const State *from, const State *to, State *state) {
	state->as<StateType>()->setX(to->as<StateType>()->getX() - from->as<StateType>()->getX());
	state->as<StateType>()->setY(to->as<StateType>()->getY() - from->as<StateType>()->getY());
	state->as<StateType>()->setZ(to->as<StateType>()->getZ() - from->as<StateType>()->getZ());
  // TODO(acauligi): what about relative velocity and acceleration?

	/*
    	Calculating relative angle is a bit tricky.
    	Refer to "interpolate" function of SO2StateSpace at line 122 of SO2StateSpace.h  in OMPL lib
        to see the original implementation in OMPL
	*/
	double diff = to->as<StateType>()->getYaw() - from->as<StateType>()->getYaw();
	if (fabs(diff) <= boost::math::constants::pi<double>()) {
        state->as<StateType>()->setYaw(diff);
  } else {
    double v;
    if (diff > 0.0) {
      diff = 2.0 * boost::math::constants::pi<double>() - diff;
    } else {
      diff = -2.0 * boost::math::constants::pi<double>() - diff;
    }

    v = - diff ;
    // input states are within bounds, so the following check is sufficient
    if (v > boost::math::constants::pi<double>()) {
      v -= 2.0 * boost::math::constants::pi<double>();
    } else {
      if (v < -boost::math::constants::pi<double>()) {
        v += 2.0 * boost::math::constants::pi<double>();
      }
    }
    state->as<StateType>()->setYaw(v);
  }

  arma::mat fcov = from->as<StateType>()->getCovariance();
  arma::mat tocov = to->as<StateType>()->getCovariance();

  if (fcov.n_rows != 0 && fcov.n_cols != 0 && tocov.n_rows != 0 && tocov.n_cols != 0 ) {
    state->as<StateType>()->setCovariance(tocov - fcov);
  }
}

arma::mat FlatQuadBeliefSpace::flatToDCM(const State *state) {
  arma::mat R(3,3);

  arma::colvec acc(4), t(3), xc(3), yc(3), zb(3), yb(3), xb(3);
  double psi = state->getYaw();

  acc = state->getAccel();
  t[0] = acc[0];
  t[1] = acc[1];
  t[2] = acc[2];    // TODO(acauligi): add gravity value
  zb = arma::normalise(t);

  xc[0] = cos(psi);   xc[1] = sin(psi); xc[2] = 0;
  yc[0] = -sin(psi);  yc[2] = cos(psi); yc[2] = 0;

  // check whether xc or yc is more normal with zb
  double xc_zb_dist = fabs(arma::cross(zb,xc));
  double yc_zb_dist = fabs(arma::cross(zb,yc));

  if (xc_zb_dixt > yc_zb_dist) {
    yb = arma::cross(zb,xc);
    yb = arma::normalise(yb);
    xb = arma::cross(yb,zb);
  } else {
    xb = arma::cross(yc,zb);
    xb = arma::normalise(xb);
    yb = arma::cross(zb,xb);
  }

  R.col(0) = xb;
  R.col(1) = yb;
  R.col(2) = zb;
  return R;
}

void FlatQuadBeliefSpace::printBeliefState(const State *state) {
  arma::colvec vel(4);
  arma::colvec acc(4);
  vel = state->as<FlatQuadBeliefSpace::StateType>()->getVelocity();
  acc = state->as<FlatQuadBeliefSpace::StateType>()->getAcceleration();

  std::cout<<"----Printing BeliefState----"<<std::endl;
  std::cout<<"State [X, Y, Z, Yaw, Velocity, Acceleration]: ";
  std::cout<<"["<<state->as<FlatQuadBeliefSpace::StateType>()->getX()<<", "<<state->as<FlatQuadBeliefSpace::StateType>()->getY()
      <<", "<<state->as<FlatQuadBeliefSpace::StateType>()->getYaw()
      <<", "<< vel[0] <<", "<< vel[1] <<", " << vel[2] <<", " << vel[3]
      <<", "<< acc[0] <<", "<< acc[1] <<", " << acc[2] <<", " << acc[3]
      <<"]"<<std::endl;
  std::cout<<"Covariance  is" <<std::endl;
  std::cout<<state->as<FlatQuadBeliefSpace::StateType>()->getCovariance()<<std::endl;
  std::cout<<"------End BeliefState-------"<<std::endl;
}
