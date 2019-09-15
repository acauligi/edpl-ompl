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
#include "ObservationModels/FlatQuadFullObservationModel.h"
#include "Visualization/Visualizer.h"
#include "Utils/FIRMUtils.h"

typename FlatQuadFullObservationModel::ObservationType 
FlatQuadFullObservationModel::getObservation(const ompl::base::State *state, bool isSimulation) {
	using namespace arma;

  ObservationType z(stateDim);

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

  colvec x_vec = state->as<FlatQuadBeliefSpace::StateType>()->getArmaData();
  z = x_vec + noise;
}

typename FlatQuadFullObservationModel::ObservationType 
FlatQuadFullObservationModel::getObservationPrediction(const ompl::base::State *state, const ObservationType& Zg) {
	using namespace arma;

  ObservationType z(stateDim);

  colvec x_vec = state->as<FlatQuadBeliefSpace::StateType>()->getArmaData();
  z = x_vec; 
 	return z;
}

typename FlatQuadFullObservationModel::JacobianType FlatQuadFullObservationModel::getObservationJacobian(const ompl::base::State *state, const ObsNoiseType& v, const ObservationType& z) {
	using namespace arma;
  mat H = eye(stateDim, stateDim);
  return H;
}

typename FlatQuadFullObservationModel::JacobianType FlatQuadFullObservationModel::getNoiseJacobian(const ompl::base::State *state, const ObsNoiseType& v, const ObservationType& z) {
	using namespace arma;
  mat H = eye(stateDim, stateDim);
  return H;
}

typename FlatQuadFullObservationModel::ObservationType FlatQuadFullObservationModel::computeInnovation(const ompl::base::State *predictedState, const ObservationType& Zg) {
	using namespace arma;

	colvec Z_pred = getObservationPrediction(predictedState, Zg);

	return Zg - Z_pred; 

}

arma::mat FlatQuadFullObservationModel::getObservationNoiseCovariance(const ompl::base::State *state, const ObservationType& z) {
	using namespace arma;

  mat R(stateDim, stateDim);
  R.eye();
  R = R*pow(this->sigma_(0),2);

  return R;
}
    

bool FlatQuadFullObservationModel::isStateObservable(const ompl::base::State *state) {
  using namespace arma;
  return true;
}

void FlatQuadFullObservationModel::loadParameters(const char *pathToSetupFile) {
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

  node = doc.FirstChild( "ObservationModels" );
  assert( node );

  TiXmlNode* child = 0;

  child = node->FirstChild("FlatQuadFullObservationModel");

  assert( child );
  itemElement = child->ToElement();
  assert( itemElement );
  double attribute_val;

  // Read in noise parameter for adding to image measurements
  itemElement->QueryDoubleAttribute("sigma_ss", &attribute_val) ;
  this->sigma_ << attribute_val << endr;

  OMPL_INFORM("FlatQuadFullObservationModel: sigma_ss = ");
  std::cout<<this->sigma_<<std::endl;
}
