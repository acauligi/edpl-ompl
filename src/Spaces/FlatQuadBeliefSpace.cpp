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
double FlatQuadBeliefSpace::StateType::reachDistPos_   = -1;
double FlatQuadBeliefSpace::StateType::reachDistOri_   = -1;
double FlatQuadBeliefSpace::StateType::reachDistCov_   = -1;
arma::colvec FlatQuadBeliefSpace::StateType::normWeights_ = arma::zeros<arma::colvec>(12);

bool FlatQuadBeliefSpace::StateType::isReached(ompl::base::State *state, bool relaxedConstraint) const {
  // check if position and orientation errors are less than thresholds
  if(!this->isReachedPose(state)) {
    return false;
  }

  // check if covariance error is less than a threshold
  if(!this->isReachedCov(state)) {
    return false;
  }

  // otherwise, the given state is considered to have reached this state
  return true;
}

bool FlatQuadBeliefSpace::StateType::isReachedWithinNEpsilon(const ompl::base::State *state, const double nEpsilon) const {
  // check if position and orientation errors are less than thresholds
  if(!this->isReachedPose(state, nEpsilon)) {
    return false;
  }

  // check if covariance error is less than a threshold
  if(!this->isReachedCov(state, nEpsilon)) {
    return false;
  }
  // otherwise, the given state is considered to have reached this state
  return true;
}

bool FlatQuadBeliefSpace::StateType::isReachedPose(const ompl::base::State *state, const double nEpsilon) const {
  // subtract the two beliefs and get the norm
  arma::colvec stateDiff = state->as<FlatQuadBeliefSpace::StateType>()->getArmaData() - this->getArmaData();

  if(stateDiff[3] > boost::math::constants::pi<double>()) {
    stateDiff[3] = (stateDiff[3] - 2*boost::math::constants::pi<double>());
  }
  if(stateDiff[3] < -boost::math::constants::pi<double>()) {
    stateDiff[3] = stateDiff[3] + 2*boost::math::constants::pi<double>();
  }

  // compute position and orientation errors
  double pos_distance_to_goal = arma::norm(stateDiff.subvec(0,2), 2);
  double ori_distance_to_goal = std::abs(stateDiff[3]);

  // check for position and orientation thresholds
  if(pos_distance_to_goal > nEpsilon * reachDistPos_) {
    return false;
  }
  if(ori_distance_to_goal > nEpsilon * reachDistOri_) {
    return false;
  }
  return true;
}

bool FlatQuadBeliefSpace::StateType::isReachedCov(const ompl::base::State *state, const double nEpsilon) const {
  // subtract the two covariances
  arma::mat covDiff = state->as<FlatQuadBeliefSpace::StateType>()->getCovariance() - this->getCovariance();
  arma::colvec covDiffDiag = covDiff.diag();

  // NOTE if the given state's covariance is already smaller than this goal state, set the difference to zero
  for (int i=0; i<covDiffDiag.size(); i++) {
    if(covDiffDiag[i] < 0.0) {
      covDiffDiag[i] = 0.0;
    }
  }

  // compute covariance error
  //double cov_distance_to_goal = arma::norm(sqrt(abs(covDiffDiag)) % normWeights_, "inf");    // deprecated
  double cov_distance_to_goal = arma::norm(abs(covDiffDiag) % normWeights_, 2);

  // check for position and orientation thresholds
  if(cov_distance_to_goal > nEpsilon * reachDistCov_) {
    return false;
  }
  return true;
}

bool FlatQuadBeliefSpace::StateType::sampleBorderBeliefState(ompl::base::State* borderBelief) const {
  // get the mean of the center belief state
  arma::colvec center_mean = this->getArmaData();

  // compute the offset for from epsilon-relaxation parameters for isReached() condition
  // 1) consider mean offset
//     arma::colvec pos_rand(2, arma::fill::randu);  // range: [ 0.0, 1.0]
//     pos_rand -= 0.5;                              // range: [-0.5, 0.5]
//     double ori_rand_dir = (std::rand()%2>0) ? +1.0 : -1.0;
//     arma::colvec pos_offset = reachDistPos_ * pos_rand / arma::norm(pos_rand, 2);
//     arma::colvec ori_offset = reachDistOri_ * arma::colvec({ori_rand_dir});
  // 2) ignore mean offset
  arma::colvec pos_offset = {0, 0, 0};
  arma::colvec ori_offset = {0};

  // set the new state property
  borderBelief->as<StateType>()->setX(center_mean[0] + pos_offset[0]);
  borderBelief->as<StateType>()->setY(center_mean[1] + pos_offset[1]);
  borderBelief->as<StateType>()->setZ(center_mean[2] + pos_offset[2]);
  borderBelief->as<StateType>()->setYaw(center_mean[3] + ori_offset[0]);

  // get the covariance of the center belief state
  arma::mat center_cov = this->getCovariance();

  // compute the offset for from epsilon-relaxation parameters for isReached() condition
  // 1) random relaxation
  // NOTE this can lead to high variance in edge cost computation with not-too-tight reachDeisCov_ value and a limited number of particles
//     arma::colvec cov_rand(3, arma::fill::randu);  // range: [ 0.0, 1.0]    // NOTE consider border belief's covariance larger, but not less, than the center belief's covariance only
  // 2) uniform relaxation

  // same proportional relaxation for each coordinate; will be weighted according to normWeights_
  arma::colvec cov_rand = {1, 1, 1,
                          1, 1, 1,
                          1, 1, 1,
                          1, 1, 1};

  arma::colvec cov_offset = reachDistCov_ * cov_rand / arma::norm(cov_rand % normWeights_, 2);


  // set the new state property
  arma::mat border_cov = center_cov;
  for (int ii=0; ii<cov_offset.size(); ii++) {
    border_cov(ii,ii) += cov_offset[ii];
  }
  borderBelief->as<StateType>()->setCovariance(border_cov);

  return true;
}

bool FlatQuadBeliefSpace::StateType::sampleTrueStateFromBelief(ompl::base::State* sampState, const double nSigma) const {
  // Cholesky decomposition such that covariance_ = transform * transform.t()
  arma::mat transform;
  if(!arma::chol(transform, covariance_, "lower")) {
    OMPL_ERROR("Failed to decompose the covariance matrix for random sampling!");
    return false;
  }

  // draw a random sample from standard normal distribution
  arma::colvec randvec(4, arma::fill::randn);

  // transform this random sample for this Gaussian distribution
  arma::colvec mean = getArmaData();
  arma::colvec randvec_transformed = mean + nSigma * transform * randvec;

  // set the new state property
  sampState->as<StateType>()->setX(randvec_transformed[0]);
  sampState->as<StateType>()->setY(randvec_transformed[1]);
  sampState->as<StateType>()->setZ(randvec_transformed[2]);
  sampState->as<StateType>()->setYaw(randvec_transformed[3]);
  sampState->as<StateType>()->setCovariance(covariance_);    // REVIEW set the covariance as the same with this state

  return true;
}

double FlatQuadBeliefSpace::StateType::getStateDistanceTo(const ompl::base::State *state) const {
  // subtract the two beliefs and get the norm
  arma::colvec stateDiff = state->as<FlatQuadBeliefSpace::StateType>()->getArmaData() - this->getArmaData();

  if(stateDiff[3] > boost::math::constants::pi<double>()) {
    stateDiff[3] = (stateDiff[3] - 2*boost::math::constants::pi<double>());
  }
  if(stateDiff[3] < -boost::math::constants::pi<double>()) {
    stateDiff[3] = stateDiff[3] + 2*boost::math::constants::pi<double>();
  }

  // compute weighted sum of position and orientation errors
  double state_distance = arma::norm(stateDiff % normWeights_, 2);

  return state_distance;
}

double FlatQuadBeliefSpace::StateType::getPosDistanceTo(const ompl::base::State *state) const {
  // subtract the two beliefs and get the norm
  arma::colvec stateDiff = state->as<FlatQuadBeliefSpace::StateType>()->getArmaData() - this->getArmaData();

  // compute position error
  return arma::norm(stateDiff.subvec(0,2), 2);
}

double FlatQuadBeliefSpace::StateType::getOriDistanceTo(const ompl::base::State *state) const {
  // subtract the two beliefs and get the norm
  arma::colvec stateDiff = state->as<FlatQuadBeliefSpace::StateType>()->getArmaData() - this->getArmaData();

  if(stateDiff[3] > boost::math::constants::pi<double>()) {
    stateDiff[3] = (stateDiff[3] - 2*boost::math::constants::pi<double>());
  }
  if(stateDiff[3] < -boost::math::constants::pi<double>()) {
    stateDiff[3] = stateDiff[3] + 2*boost::math::constants::pi<double>();
  }

  // compute orientation error
  return std::abs(stateDiff[3]);
}

bool FlatQuadBeliefSpace::StateType::mergeBeliefIntoThis(const ompl::base::State *newBelief) {
  // NOTE this function should be called before incrementing N(h) by +1

  // retrieve means and covariances
  arma::mat covThis = this->getCovariance();
  arma::colvec meanThis = this->getArmaData();
  arma::mat covNew = newBelief->as<StateType>()->getCovariance();
  arma::colvec meanNew = newBelief->as<StateType>()->getArmaData();

  // compute natural parameters
  arma::mat invCovThis = arma::inv(covThis);
  arma::mat invCovNew = arma::inv(covNew);
  arma::colvec invMeanThis = invCovThis * meanThis;
  arma::colvec invMeanNew = invCovNew * meanNew;

  // compute merged parameter values
  double nVisitThis = this->getThisQVvisit();
  arma::mat invCovMerged = invCovThis + (invCovNew - invCovThis) / (nVisitThis + 1);
  arma::colvec invMeanMerged = invMeanThis + (invMeanNew - invMeanThis) / (nVisitThis + 1);

  // convert to cannocial parameters
  arma::mat covMerged = arma::inv(invCovMerged);
  arma::colvec meanMerged = covMerged * invMeanMerged;

  // apply the change of belief state
  this->setCovariance(covMerged);
  this->setArmaData(meanMerged);

  return true;
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
  double dz = state1->as<StateType>()->getZ() - state2->as<StateType>()->getZ();

  arma::colvec stateDiff = state1->as<FlatQuadBeliefSpace::StateType>()->getArmaData() - 
    state2->as<FlatQuadBeliefSpace::StateType>()->getArmaData();

  // TODO(acauligi): use full 12D state distance from 2PBVP?
  return arma::norm(stateDiff.subvec(0,2), 2);
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

int FlatQuadBeliefSpace::StateType::flatTransformMethod() const {
  arma::colvec acc(4), t(3), xc(3), yc(3), zb(3), yb(3), xb(3);
  // double psi = state->as<StateType>()->getYaw();
  double psi = this->getYaw(); 

  // acc = state->as<StateType>()->getAcceleration();
  acc = this->getAcceleration();
  t[0] = acc[0];
  t[1] = acc[1];
  t[2] = acc[2] + gravity_;
  zb = arma::normalise(t);

  xc[0] = cos(psi);   xc[1] = sin(psi); xc[2] = 0;
  yc[0] = -sin(psi);  yc[2] = cos(psi); yc[2] = 0;

  // check whether xc or yc is more normal with zb
  double xc_zb_dist = arma::norm(arma::cross(zb,xc));
  double yc_zb_dist = arma::norm(arma::cross(zb,yc));

  if (xc_zb_dist > yc_zb_dist) {
    yb = arma::cross(zb,xc);
    yb = arma::normalise(yb);
    xb = arma::cross(yb,zb);
    return ZYX;
  }
  return ZXY;
}

arma::mat FlatQuadBeliefSpace::StateType::flatToDCM() const {
  // returns rotation matrix from body to world frame
  arma::mat Rwb(3,3);
  arma::colvec acc(4), t(3), xb(3), yb(3), zb(3);
  // double psi = state->as<StateType>()->getYaw();
  double psi = this->getYaw(); 

  // acc = state->as<StateType>()->getAcceleration();
  acc = this->getAcceleration(); 
  t[0] = acc[0];
  t[1] = acc[1];
  t[2] = acc[2] + gravity_;
  zb = arma::normalise(t);

  // Check which version of differential flatness transform to use
  int flat_transform_method = this->flatTransformMethod();
  if (flat_transform_method == ZYX) {
    arma::colvec xc(3);
    xc[0] = cos(psi);   xc[1] = sin(psi); xc[2] = 0;
    yb = arma::cross(zb,xc);
    yb = arma::normalise(yb);
    xb = arma::cross(yb,zb);
  } else {
    arma::colvec yc(3);
    yc[0] = -sin(psi);  yc[2] = cos(psi); yc[2] = 0;
    xb = arma::cross(yc,zb);
    xb = arma::normalise(xb);
    yb = arma::cross(zb,xb);
  }

  Rwb.col(0) = xb;
  Rwb.col(1) = yb;
  Rwb.col(2) = zb;
  return Rwb;
}

arma::colvec FlatQuadBeliefSpace::StateType::flatToAxisAngle() const {
  // returns 4D vector with [ax, ay, az, angle] where (ax,ay,az) are axis of rotation
  arma::mat Rwb(3,3);
  Rwb = this->flatToDCM();
  arma::colvec ax_angle = arma::zeros<arma::colvec>(4);

  // Eq. 102 and 103a in "A Survey of Attitude Representations", M.D. Shuster (1993) 
  double th = acos(0.5*(arma::trace(Rwb)-1));
  if (th > boost::math::constants::pi<double>()) {
    th -= 2*boost::math::constants::pi<double>();
  }
  if(th < -boost::math::constants::pi<double>()) {
    th += 2*boost::math::constants::pi<double>();
  }

  // Ensure sin(th) != 0
  double sin_th = sin(th);
  if (abs(sin_th) > 1e-3) {
    ax_angle[0] = 1/(2*sin_th) * (Rwb[1,2]-Rwb[2,1]);
    ax_angle[1] = 1/(2*sin_th) * (Rwb[2,0]-Rwb[0,2]);
    ax_angle[2] = 1/(2*sin_th) * (Rwb[0,1]-Rwb[1,0]);
  }

  ax_angle = arma::normalise(ax_angle);
  ax_angle[3] = th;

  return ax_angle;
}

ompl::base::StateSpacePtr FlatQuadBeliefSpace::StateType::flatToOMPLSO3() const {
  // returns ompl SO3StateSpace object corresponding to orientation 
  arma::mat Rwb(3,3);
  Rwb = this->flatToDCM();

  // TODO(acauligi)
  arma::colvec ax(4);
  ax = this->flatToAxisAngle();

  ompl::base::StateSpacePtr rot(new ompl::base::SO3StateSpace());
  // rot->as<ompl::base::SO3StateSpace*>()->setAxisAngle(ax[0], ax[1], ax[2], ax[3]);
  return rot;
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
