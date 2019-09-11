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
#include "ObservationModels/LandmarkObservationModel.h"
#include "Visualization/Visualizer.h"
#include "Utils/FIRMUtils.h"

typename LandmarkObservationModel::ObservationType 
LandmarkObservationModel::getObservation(const ompl::base::State *state, bool isSimulation) {
	using namespace arma;

  ObservationType z(landmarkInfoDim*landmarks_.size());

  //generate observation from state, and corrupt with the given noise
  for(unsigned int ii = 0; ii < landmarks_.size(); ii++) {

    colvec noise = zeros<colvec>(obsNoiseDim);

    if(isSimulation) {
      // generate gaussian noise
      // get standard deviations of noise (sqrt of covariance matrix)
      // extract state from Cfg and normalize
      // generate noise scaling/shifting factor
      // generate raw noise
      colvec randNoiseVec = randn<colvec>(obsNoiseDim);

      // generate noise from a distribution scaled and shifted from
      // normal distribution N(0,1) to N(0,eta*range + sigma)
      // (shifting was done in GetNoiseCovariance)
      noise = sigma_%randNoiseVec;
    }

    colvec image_meas = zeros<colvec>(landmarkInfoDim);
    if (getPerspectiveProjection(state, landmarks_[ii], image_meas)) {
      z.subvec(landmarkInfoDim*ii,landmarkInfoDim*ii+1) = image_meas + noise;
    }
  }
 	return z;
}

typename LandmarkObservationModel::ObservationType 
LandmarkObservationModel::getObservationPrediction(const ompl::base::State *state, const ObservationType& Zg) {
	using namespace arma;

  ObservationType z(landmarkInfoDim*landmarks_.size());

  //generate observation from state
  for(unsigned int ii = 0; ii < landmarks_.size(); ii++) {
    colvec image_meas = zeros<colvec>(landmarkInfoDim);
    if (getPerspectiveProjection(state, landmarks_[ii], image_meas)) {
      z.subvec(landmarkInfoDim*ii,landmarkInfoDim*ii+1) = image_meas; 
    }
  }

 	return z;
}

typename LandmarkObservationModel::JacobianType LandmarkObservationModel::getObservationJacobian(const ompl::base::State *state, const ObsNoiseType& v, const ObservationType& z) {
	using namespace arma;

  unsigned int number_of_landmarks = landmarks_.size();

  colvec xVec = state->as<FlatQuadBeliefSpace::StateType>()->getArmaData();

  mat H(2*number_of_landmarks, 12); // Since we are passing the common id list

  for(unsigned int ii=0; ii<number_of_landmarks; ii++) {
    mat H_i(2,12);
    
    colvec lk_w = landmarks_[ii];
    mat Rbw_(3,3);    // world to body SO(3) rotation matrix
    colvec pbw_(3);   // world to body R3 position

    // lmk coordinate transformed into camera frame and multiplied by intrinsics matrix K
    colvec lk_c = this->K_ *(this->Rcb_ * Rbw_ * this->Rcb_ * pbw_ + this->pcb_);

    // generate list of matrices corresponding to derivatives of columns of Rbw_ w.r.t. x
    // i.e. each element of dRbw_dx is 3x12 matrix
    std::vector<mat> dRbw_dx;
    for (size_t jj=0; jj<3; jj++) {
      dRbw_dx.push_back(zeros(3,12));
      this->getObservationHx(state,dRbw_dx[jj],jj);
    }
    
    std::vector<colvec> dpbw_dx;
    for (size_t jj=0; jj<3; jj++) {
      colvec local_grad = zeros(12); 
      local_grad[jj] = 1.;
      dpbw_dx.push_back(local_grad);
    }

    mat dlk_c_dx = zeros(3,12);
    for (int jj=0; jj<3; jj++) {
      for (int kk=0; kk<3; kk++) {
        dlk_c_dx.row(jj) += (lk_c[kk]*this->K_.row(jj)*this->Rcb_ * dRbw_dx[kk] + dpbw_dx[kk] * K_.row(jj) * Rbw_.col(kk));
      }
    }

    H_i.row(0) = (dlk_c_dx.row(0)*lk_c[2] - lk_c[0]*dlk_c_dx.row(2)) / (lk_c[2]*lk_c[2]);
    H_i.row(1) = (dlk_c_dx.row(1)*lk_c[2] - lk_c[1]*dlk_c_dx.row(2)) / (lk_c[2]*lk_c[2]);
    H.submat(2*ii, 0, 2*ii+1, 11) = H_i;
  }
  return H;
}

void LandmarkObservationModel::getObservationHx(const ompl::base::State *state, arma::mat& Hx, size_t coord_idx) {
  // derivative of Rbw(:,coord_idx) w.r.t. x
  using namespace arma;

  Hx = zeros(3,12);
  colvec x = state->as<FlatQuadBeliefSpace::StateType>()->getArmaData();
  double z1 = x[0]; 
  double z2 = x[1]; 
  double z3 = x[2]; 
  double z4 = x[3]; 
  double z5 = x[4]; 
  double z6 = x[5]; 
  double z7 = x[6]; 
  double z8 = x[7]; 
  double z9 = x[8]; 
  double z10 = x[9]; 
  double z11 = x[10]; 
  double z12 = x[11]; 

  if (coord_idx==0) {
    Hx.row(0) = rowvec({  0., 0., -(pow(z10,2)*sin(z3) + pow(z11,2)*sin(z3) + z9*z10*cos(z3))/(pow(z9,2) + pow(z10,2) + pow(z11,2)), 0., 0., 0., 0., 0.,  -(pow(z10,3)*sin(z3) + 2*z9*pow(z10,2)*cos(z3) + 2*z9*pow(z11,2)*cos(z3) - pow(z9,2)*z10*sin(z3) + z10*pow(z11,2)*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),2),  -(z9*(pow(z9,2)*sin(z3) - pow(z10,2)*sin(z3) + pow(z11,2)*sin(z3) - 2*z9*z10*cos(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2), (2*z9*z11*(z9*cos(z3) + z10*sin(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),2), 0. } );
    Hx.row(1) = rowvec({  0., 0., -(z11*cos(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),0.5), 0., 0., 0., 0., 0., (z9*z11*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),1.5), (z10*z11*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), -(sin(z3)*(pow(z9,2) + pow(z10,2)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),1.5), 0. });
    Hx.row(2) = rowvec({  0., 0., 0., 0., 0., 0., 0., 0., (pow(z10,2) + pow(z11,2))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),1.5), -(z9*z10)/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),1.5), -(z9*z11)/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),1.5), 0.});

    // Jacobian for second way of recovering Rbw
    if (false) {
      Hx.row(0) = rowvec( { 0, 0, -(z11*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 0.5), 0, 0, 0, 0, 0, -(z9*z11*cos(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), -(z10*z11*cos(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), (cos(z3)*(pow(z9,2) + pow(z10,2)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), 0.});
      Hx.row(1) = rowvec( { 0, 0, -(pow(z10,2)*cos(z3) + pow(z11,2)*cos(z3) - z9*z10*sin(z3))/(pow(z9,2) + pow(z10,2) + pow(z11,2)), 0, 0, 0, 0, 0, (- pow(z10,3)*cos(z3) + pow(z9,2)*z10*cos(z3) - z10*pow(z11,2)*cos(z3) + 2*z9*pow(z10,2)*sin(z3) + 2*z9*pow(z11,2)*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2),  -(z9*(pow(z9,2)*cos(z3) - pow(z10,2)*cos(z3) + pow(z11,2)*cos(z3) + 2*z9*z10*sin(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2), (2*z9*z11*(z10*cos(z3) - z9*sin(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2), 0.}); 
      Hx.row(2) = rowvec( { 0, 0, 0, 0, 0, 0, 0, 0, (pow(z10,2) + pow(z11,2))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), -(z9*z10)/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), -(z9*z11)/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), 0.});
    }
  } else if (coord_idx==1) {
    Hx.row(0) = rowvec( { 0., 0., (pow(z9,2)*cos(z3) + pow(z11,2)*cos(z3) + z9*z10*sin(z3))/(pow(z9,2) + pow(z10,2) + pow(z11,2)), 0., 0., 0., 0., 0., -(pow(z10,3)*cos(z3) - pow(z9,2)*z10*cos(z3) + z10*pow(z11,2)*cos(z3) - 2*z9*pow(z10,2)*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2), -(pow(z9,3)*cos(z3) - z9*pow(z10,2)*cos(z3) + z9*pow(z11,2)*cos(z3) + 2*pow(z9,2)*z10*sin(z3) + 2*z10*pow(z11,2)*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2), (2*z10*z11*(z9*cos(z3) + z10*sin(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2), 0});
    Hx.row(1) = rowvec( { 0., 0., -(z11*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),0.5), 0., 0., 0., 0., 0., -(z9*z11*cos(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), -(z10*z11*cos(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),1.5), (cos(z3)*(pow(z9,2) + pow(z10,2)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),1.5), 0.});
    Hx.row(2) = rowvec( { 0., 0., 0., 0., 0., 0., 0., 0., -(z9*z10)/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),1.5), (pow(z9,2) + pow(z11,2))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), -(z10*z11)/pow(pow(z9,2) + pow(z10,2) + pow(z11,2),1.5), 0.});

    // Jacobian for second way of recovering Rbw
    if (false) {
      Hx.row(0) = rowvec( { 0., 0., (z11*cos(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 0.5), 0., 0., 0., 0., 0., -(z9*z11*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), -(z10*z11*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), (sin(z3)*(pow(z9,2) + pow(z10,2)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), 0.});
      Hx.row(1) = rowvec( { 0., 0., -(pow(z9,2)*sin(z3) + pow(z11,2)*sin(z3) - z9*z10*cos(z3))/(pow(z9,2) + pow(z10,2) + pow(z11,2)), 0., 0., 0., 0., 0.,  (z10*(- pow(z9,2)*sin(z3) + pow(z10,2)*sin(z3) + pow(z11,2)*sin(z3) + 2*z9*z10*cos(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2),  -(- pow(z9,3)*sin(z3) + 2*pow(z9,2)*z10*cos(z3) + 2*z10*pow(z11,2)*cos(z3) + z9*pow(z10,2)*sin(z3) - z9*pow(z11,2)*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2),  (2*z10*z11*(z10*cos(z3) - z9*sin(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2), 0});
      Hx.row(2) = rowvec( { 0., 0., 0., 0., 0., 0., 0., 0., -(z9*z10)/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), (pow(z9,2) + pow(z11,2))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), -(z10*z11)/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), });
    }
  } else {
    // coord_idx=2
    Hx.row(0) = rowvec({ 0., 0., -(z11*(z10*cos(z3) - z9*sin(z3)))/(pow(z9,2) + pow(z10,2) + pow(z11,2)), 0, 0, 0, 0, 0, (z11*(pow(z9,2)*cos(z3) - pow(z10,2)*cos(z3) - pow(z11,2)*cos(z3) + 2*z9*z10*sin(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2), -(z11*(pow(z9,2)*sin(z3) - pow(z10,2)*sin(z3) + pow(z11,2)*sin(z3) - 2*z9*z10*cos(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2),  -((z9*cos(z3) + z10*sin(z3))*(pow(z9,2) + pow(z10,2) - pow(z11,2)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2), 0});
    Hx.row(1) = rowvec({ 0, 0, (z9*cos(z3) + z10*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 0.5), 0, 0, 0, 0, 0, (pow(z10,2)*sin(z3) + pow(z11,2)*sin(z3) + z9*z10*cos(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), -(pow(z9,2)*cos(z3) + pow(z11,2)*cos(z3) + z9*z10*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), (z11*(z10*cos(z3) - z9*sin(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), 0});
    Hx.row(2) = rowvec({ 0., 0., 0., 0., 0., 0., 0., 0, -(z9*z11)/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), -(z10*z11)/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), (pow(z9,2) + pow(z10,2))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), 0.});

    // Jacobian for second way of recovering Rbw
    if (false) {
      Hx.row(0) = rowvec({ 0., 0., -(z10*cos(z3) - z9*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 0.5), 0., 0., 0., 0., 0., -(pow(z10,2)*cos(z3) + pow(z11,2)*cos(z3) - z9*z10*sin(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), -(pow(z9,2)*sin(z3) + pow(z11,2)*sin(z3) - z9*z10*cos(z3))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), (z11*(z9*cos(z3) + z10*sin(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), 0.});
      Hx.row(1) = rowvec({ 0., 0., (z11*(z9*cos(z3) + z10*sin(z3)))/(pow(z9,2) + pow(z10,2) + pow(z11,2)), 0., 0., 0., 0., 0., (z11*(- pow(z9,2)*sin(z3) + pow(z10,2)*sin(z3) + pow(z11,2)*sin(z3) + 2*z9*z10*cos(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2), -(z11*(pow(z9,2)*cos(z3) - pow(z10,2)*cos(z3) + pow(z11,2)*cos(z3) + 2*z9*z10*sin(z3)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2), -((z10*cos(z3) - z9*sin(z3))*(pow(z9,2) + pow(z10,2) - pow(z11,2)))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 2), 0.});
      Hx.row(2) = rowvec({ 0., 0., 0., 0., 0., 0., 0., 0., -(z9*z11)/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), -(z10*z11)/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), (pow(z9,2) + pow(z10,2))/pow(pow(z9,2) + pow(z10,2) + pow(z11,2), 1.5), 0.});
    }
  }
}

typename LandmarkObservationModel::JacobianType LandmarkObservationModel::getNoiseJacobian(const ompl::base::State *state, const ObsNoiseType& v, const ObservationType& z) {
	using namespace arma;

  unsigned int number_of_landmarks = landmarks_.size() ;

  mat M(2*number_of_landmarks,2*number_of_landmarks);
  M.eye();
  return M;
}

typename LandmarkObservationModel::ObservationType LandmarkObservationModel::computeInnovation(const ompl::base::State *predictedState, const ObservationType& Zg) {
	using namespace arma;

	colvec Z_pred = getObservationPrediction(predictedState, Zg);

	return Zg - Z_pred; 

}

arma::mat LandmarkObservationModel::getObservationNoiseCovariance(const ompl::base::State *state, const ObservationType& z) {
	using namespace arma;

  unsigned int number_of_landmarks = landmarks_.size();

  mat R(2*number_of_landmarks,2*number_of_landmarks);
  R.eye();
  R = R*pow(this->sigma_(0),2);

  return R;
}
    
bool LandmarkObservationModel::getPerspectiveProjection(const ompl::base::State *state, const arma::colvec& landmark, arma::colvec& image_meas) {
	using namespace arma;
	colvec x_vec = state->as<FlatQuadBeliefSpace::StateType>()->getArmaData();

  // position of quad pbw and orientation of world w.r.t. quad Rwb 
  colvec pbw = x_vec.subvec(0,2); 
	mat Rwb(3,3); 
	// mat Rwb = state->as<FlatQuadBeliefSpace::StateType>()->flatToDCM();
	
  colvec lmk_world = this->K_ * (this->Rcb_* Rwb.t() * landmark + this->Rcb_ * pbw + this->pcb_);
  image_meas[0] = lmk_world[0] / lmk_world[2];
  image_meas[1] = lmk_world[1] / lmk_world[2];

  // TODO(acauligi): add check to determine whether feature is within viewing frustum 
  return true;
}

bool LandmarkObservationModel::isStateObservable(const ompl::base::State *state) {
  using namespace arma;

  colvec obs = this->getObservation(state, false);

  if(obs.n_rows > 2) {
    return true;
  }

  return false;

}

void LandmarkObservationModel::loadLandmarks(const char *pathToSetupFile) {
	using namespace arma;
	// Load XML containing landmarks
	TiXmlDocument doc(pathToSetupFile);
	bool loadOkay = doc.LoadFile();

	if (!loadOkay) {
		printf( "Could not load Landmark list . Error='%s'. Exiting.\n", doc.ErrorDesc() );
		exit(1);
	}

	TiXmlNode* node = 0;
	TiXmlElement* landmarkElement = 0;
	TiXmlElement* itemElement = 0;

	// Get the landmarklist node
	node = doc.FirstChild( "LandmarkList" );
	assert( node );
	landmarkElement = node->ToElement(); //convert node to element
	assert( landmarkElement  );

	TiXmlNode* child = 0;

	//Iterate through all the landmarks and put them into the "landmarks_" list
	while( (child = landmarkElement ->IterateChildren(child))) {
		assert( child );
		itemElement = child->ToElement();
		assert( itemElement );

		ObservationType landmark(1+landmarkInfoDim);    // +1 for ID tag
		landmark.zeros();
		double attribute_val;
		itemElement->QueryDoubleAttribute("id", &attribute_val) ;
		landmark[0] = attribute_val;
		itemElement->QueryDoubleAttribute("x", &attribute_val) ;
		landmark[1] = attribute_val;
		itemElement->QueryDoubleAttribute("y", &attribute_val) ;
		landmark[2] = attribute_val;
		itemElement->QueryDoubleAttribute("z", &attribute_val) ;
		landmark[3] = attribute_val;

		this->landmarks_.push_back(landmark);
	}

	OMPL_INFORM("LandmarkObservationModel: Total number of landmarks loaded successfully : %u", landmarks_.size() );

	Visualizer::addLandmarks(landmarks_);
}

void LandmarkObservationModel::loadParameters(const char *pathToSetupFile) {
	using namespace arma;
  // Load XML containing landmarks
  TiXmlDocument doc(pathToSetupFile);
  bool loadOkay = doc.LoadFile();

  if ( !loadOkay ) {
    printf( "Could not load setup file . Error='%s'. Exiting.\n", doc.ErrorDesc() );
    exit( 1 );
  }

  TiXmlNode* node = 0;

  TiXmlElement* itemElement = 0;

  // Get the landmarklist node
  node = doc.FirstChild( "ObservationModels" );
  assert( node );

  TiXmlNode* child = 0;

  child = node->FirstChild("LandmarkObservationModel");

  //Iterate through all the landmarks and put them into the "landmarks_" list
  assert( child );
  itemElement = child->ToElement();
  assert( itemElement );
  double attribute_val;

  // Read camera intrinsics matrix parameters 
  this->K_ = arma::zeros(3,3);
  itemElement->QueryDoubleAttribute("K_alpha_u", &attribute_val);
  this->K_[0,0] = attribute_val;
  itemElement->QueryDoubleAttribute("K_gamma", &attribute_val);
  this->K_[0,1] = attribute_val;
  itemElement->QueryDoubleAttribute("K_u0", &attribute_val);
  this->K_[0,2] = attribute_val;
  itemElement->QueryDoubleAttribute("K_alpha_v", &attribute_val);
  this->K_[1,1] = attribute_val;
  itemElement->QueryDoubleAttribute("K_v0", &attribute_val);
  this->K_[1,2] = attribute_val;
  this->K_[2,2] = 1;

  // Read body to camera rotation matrix 
  this->Rcb_ = arma::zeros(3,3);
  itemElement->QueryDoubleAttribute("Rcb_00", &attribute_val);
  this->Rcb_[0,0] = attribute_val;
  itemElement->QueryDoubleAttribute("Rcb_01", &attribute_val);
  this->Rcb_[0,1] = attribute_val;
  itemElement->QueryDoubleAttribute("Rcb_02", &attribute_val);
  this->Rcb_[0,2] = attribute_val;
  itemElement->QueryDoubleAttribute("Rcb_10", &attribute_val);
  this->Rcb_[1,0] = attribute_val;
  itemElement->QueryDoubleAttribute("Rcb_11", &attribute_val);
  this->Rcb_[1,1] = attribute_val;
  itemElement->QueryDoubleAttribute("Rcb_12", &attribute_val);
  this->Rcb_[1,2] = attribute_val;
  itemElement->QueryDoubleAttribute("Rcb_20", &attribute_val);
  this->Rcb_[2,0] = attribute_val;
  itemElement->QueryDoubleAttribute("Rcb_21", &attribute_val);
  this->Rcb_[2,1] = attribute_val;
  itemElement->QueryDoubleAttribute("Rcb_22", &attribute_val);
  this->Rcb_[2,2] = attribute_val;

  // Read body to camera translation 
  this->pcb_ = arma::zeros(3);
  itemElement->QueryDoubleAttribute("pcb_x", &attribute_val);
  this->pcb_[0] = attribute_val;
  itemElement->QueryDoubleAttribute("pcb_y", &attribute_val);
  this->pcb_[1] = attribute_val;
  itemElement->QueryDoubleAttribute("pcb_z", &attribute_val);
  this->pcb_[2] = attribute_val;

  // Read in noise parameter for adding to image measurements
  itemElement->QueryDoubleAttribute("sigma_ss", &attribute_val) ;
  this->sigma_ << attribute_val << endr;

  OMPL_INFORM("LandmarkObservationModel: sigma_ss = ");
  std::cout<<this->sigma_<<std::endl;
}