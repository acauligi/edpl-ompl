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

/* Author: Abhishek Cauligi & Saurav Agarwal */

#ifndef FLAT_QUAD_SETUP_H
#define FLAT_QUAD_SETUP_H

#include <tinyxml.h>
#include <omplapp/geometry/RigidBodyGeometry.h>
#include <ompl/base/samplers/MaximizeClearanceValidStateSampler.h>
#include <ompl/base/samplers/MinimumClearanceValidStateSampler.h>
#include "omplapp/geometry/detail/FCLStateValidityChecker.h"

#include "Planner/FIRM.h"
#include "edplompl.h"
#include "Visualization/Window.h"
#include "Visualization/Visualizer.h"

/** \brief Wrapper for ompl::app::RigidBodyPlanning that plans for rigid bodies in FlatQuadBeliefSpace for a MAV using FIRM */
class FlatQuadSetup : public ompl::app::RigidBodyGeometry {
    typedef FlatQuadBeliefSpace::StateType StateType;

public:

    FlatQuadSetup() : ompl::app::RigidBodyGeometry(ompl::app::Motion_2D, ompl::app::FCL),
    ss_(ompl::base::StateSpacePtr(new SE2BeliefSpace())) {
      // set static variables
      // Not in use as of now
      Controller<FiniteTimeLQR, ExtendedKF>::setNodeReachedAngle(10); // degrees
      Controller<FiniteTimeLQR, ExtendedKF>::setNodeReachedDistance(0.1);// meters
      Controller<FiniteTimeLQR, ExtendedKF>::setMaxTries(30);
      Controller<FiniteTimeLQR, ExtendedKF>::setMaxTrajectoryDeviation(0.5); // meters

      // [2] NodeController (typedef in include/Planner/FIRM.h)
      Controller<StationaryLQR, LinearizedKF>::setNodeReachedAngle(10); // degrees
      Controller<StationaryLQR, LinearizedKF>::setNodeReachedDistance(0.1);// meters
//         Controller<StationaryLQR, LinearizedKF>::setMaxTries(30);
      Controller<StationaryLQR, LinearizedKF>::setMaxTries(300);
//         Controller<StationaryLQR, LinearizedKF>::setMaxTrajectoryDeviation(0.5); // meters
//         Controller<StationaryLQR, LinearizedKF>::setMaxTrajectoryDeviation(2.25); // meters
      Controller<StationaryLQR, LinearizedKF>::setMaxTrajectoryDeviation(2.50); // meters

      // [1] EdgeController (typedef in include/Planner/FIRM.h)
      RHCICreate::setControlQueueSize(5);
      RHCICreate::setTurnOnlyDistance(0.01);
      Controller<RHCICreate, ExtendedKF>::setNodeReachedAngle(10.0); // degrees
//         Controller<RHCICreate, ExtendedKF>::setNodeReachedDistance(0.2);// meters
      Controller<RHCICreate, ExtendedKF>::setNodeReachedDistance(0.1);// meters
      Controller<RHCICreate, ExtendedKF>::setMaxTries(30);
//         Controller<RHCICreate, ExtendedKF>::setMaxTrajectoryDeviation(0.5); // meters
      Controller<RHCICreate, ExtendedKF>::setMaxTrajectoryDeviation(1.0); // meters

      // setting the mean and norm weights (used in reachability check)
      // NOTE these values will be overwritten by loadParameters()
      StateType::covNormWeight_  =  1.0;
      StateType::meanNormWeight_ =  2.0;
      StateType::reachDist_ =  0.009;    // 0.015     // distance threshold for position, orientation, and covariance
      StateType::reachDistPos_ = 0.1;    // meters    // distance threshold for position
      StateType::reachDistOri_ = 10.0/180.0*boost::math::constants::pi<double>();    // radian    // distance threshold for orientation
      StateType::reachDistCov_ = 0.0004;    // distance threshold for covariance

      // set the state component norm weights
      arma::colvec normWeights(3);
      normWeights(0) = 2.0/std::sqrt(9);
      normWeights(1) = 2.0/std::sqrt(9);
      normWeights(2) = 1.0/std::sqrt(9);
      StateType::normWeights_ = normWeights;

      // The bounds should be inferred from the geometry files,
      // there is a function in Apputils to do this, so use that.        
      ompl::base::RealVectorBounds bounds(2);
      
      // set X & Y bound
      bounds.setLow(0.0);
      bounds.setHigh(20.0);

      ss_->as<SE2BeliefSpace>()->setBounds(bounds);

      //Construct the control space
      ompl::control::ControlSpacePtr controlspace( new ompl::control::RealVectorControlSpace(ss_,3) );

      cs_ = controlspace;

      // construct an instance of space information from this state space
      firm::SpaceInformation::SpaceInformationPtr si(new firm::SpaceInformation(ss_, cs_));
      siF_ = si;

      //siF_->setValidStateSamplerAllocator(FlatQuadSetup::allocMaxClearanceValidStateSampler);

      ompl::base::ProblemDefinitionPtr prblm(new ompl::base::ProblemDefinition(siF_));

      pdef_ = prblm;

      start_ = siF_->allocState();

      setup_ = false;

      dynamicObstacles_ = false;

      plannerMethod_ = 0; // by default we use FIRM
    }

    virtual ~FlatQuadSetup(void) {
    }

    const firm::SpaceInformation::SpaceInformationPtr& getSpaceInformation() const {
      return siF_;
    }

    void setPathToSetupFile(const std::string &path) {
      pathToSetupFile_  = path;       
    }

    void setStartState(const double X, const double Y) {
      ompl::base::State *temp = siF_->allocState();
      temp->as<StateType>()->setXY(X,Y);
      siF_->copyState(start_, temp);
      siF_->freeState(temp);
    }

    void addGoalState(const double X, const double Y) {
      ompl::base::State *temp = siF_->allocState();
      temp->as<StateType>()->setXY(X,Y);
      goalList_.push_back(temp);
    }


    virtual void setup() {
      if(!setup_) {
          this->loadParameters();

          if(pathToSetupFile_.length() == 0) {
            throw ompl::Exception("Path to setup file not set!");
          }

          if(!hasEnvironment() || !hasRobot()) {
            throw ompl::Exception("Robot/Environment mesh files not setup!");
          }

          ss_->as<SE2BeliefSpace>()->setBounds(inferEnvironmentBounds());

          // Create an FCL state validity checker and assign to space information
          const ompl::base::StateValidityCheckerPtr &fclSVC = this->allocStateValidityChecker(siF_, getGeometricStateExtractor(), false);
          siF_->setStateValidityChecker(fclSVC);

          // provide the observation model to the space
          ObservationModelMethod::ObservationModelPointer om(new HeadingBeaconObservationModel(siF_, pathToSetupFile_.c_str()));
          siF_->setObservationModel(om);

          // Provide the motion model to the space
          // We use the omnidirectional model because collision checking requires SE2
          MotionModelMethod::MotionModelPointer mm(new OmnidirectionalMotionModel(siF_, pathToSetupFile_.c_str()));            
          siF_->setMotionModel(mm);

          ompl::control::StatePropagatorPtr prop(ompl::control::StatePropagatorPtr(new OmnidirectionalStatePropagator(siF_)));
          statePropagator_ = prop;
          siF_->setStatePropagator(statePropagator_);
          siF_->setPropagationStepSize(0.1); // this is the duration that a control is applied
          siF_->setStateValidityCheckingResolution(0.005);
          siF_->setMinMaxControlDuration(1,100);

          if(!start_ || goalList_.size() == 0) {
            throw ompl::Exception("Start/Goal not set");
          }

          pdef_->setStartAndGoalStates(start_, goalList_[0], 1.0);

          ompl::base::PlannerPtr plnr(new FIRM(siF_, false));

          planner_ = plnr;

          planner_->setProblemDefinition(pdef_);

          planner_->as<FIRM>()->setMinFIRMNodes(minNodes_);

          planner_->as<FIRM>()->setMaxFIRMNodes(maxNodes_);

          planner_->as<FIRM>()->setKidnappedState(kidnappedState_);

          planner_->as<FIRM>()->loadParametersFromFile(pathToSetupFile_.c_str());
          
          planner_->setup();            

          // Setup visualizer because it is needed while loading roadmap and visualizing it
          Visualizer::updateSpaceInformation(this->getSpaceInformation());

          Visualizer::updateRenderer(*dynamic_cast<const ompl::app::RigidBodyGeometry*>(this), this->getGeometricStateExtractor());

          if (useSavedRoadMap_ == 1) planner_->as<FIRM>()->loadRoadMapFromFile(pathToRoadMapFile_.c_str());

          setup_ = true;
      }
    }

    ompl::base::PlannerStatus solve() {
      if(!setup_) {
        this->setup();
      }
      return planner_->solve(planningTime_);
    }

    virtual void executeSolution(int choice=0) {
      switch(choice) {
        case 0:
          planner_->as<FIRM>()->executeFeedback();
          break;

        case 1:
          planner_->as<FIRM>()->executeFeedbackWithRollout();
          break;

        case 2:
          planner_->as<FIRM>()->executeFeedbackWithKidnapping();
          break;

        default:
          OMPL_ERROR("PlanningMode method %d is not valid... Check the setup file!", choice);
          exit(0);
          break;
      }
    }

    void  Run() {
      executeSolution(plannerMethod_);

      // Need a function to get terminated state
      if(goalList_.size() > 1) {
        updateEnvironmentMesh();

        for(int i=0; i < goalList_.size()-1;i++) {
          pdef_->setStartAndGoalStates(goalList_[i], goalList_[i+1], 1.0);

          planner_->setProblemDefinition(pdef_);

          if(this->solve()) {
            executeSolution(plannerMethod_);
          }
        }
      }
    }

    virtual void updateEnvironmentMesh(int obindx = 0) {
      if(dynamicObstacles_) {
        // Set environment to new mesh with some dynamic / additional obstacles
        if(!this->setEnvironmentMesh(dynObstList_[obindx])) {
          OMPL_ERROR("Couldn't set mesh with path: %s",dynObstList_[obindx]);
        }
        
        const ompl::base::StateValidityCheckerPtr &svc = std::make_shared<ompl::app::FCLStateValidityChecker<ompl::app::Motion_2D>>(siF_,  getGeometrySpecification(), getGeometricStateExtractor(), false);

        siF_->setStateValidityChecker(svc);

        planner_->as<FIRM>()->updateCollisionChecker(svc);
      }
    }

    ompl::app::GeometricStateExtractor getGeometricStateExtractor(void) const {
      return boost::bind(&FlatQuadSetup::getGeometricComponentStateInternal, this, _1, _2);
    }

    virtual void saveRoadmap() {
      planner_->as<FIRM>()->savePlannerData();
    }

protected:

    static ompl::base::ValidStateSamplerPtr allocMaxClearanceValidStateSampler(const ompl::base::SpaceInformation *si) {
      // we can perform any additional setup / configuration of a sampler here,
      // but there is nothing to tweak in case of the ObstacleBasedValidStateSampler.
      std::shared_ptr<ompl::base::MinimumClearanceValidStateSampler> ss = std::make_shared<ompl::base::MinimumClearanceValidStateSampler>(si);
      
      ss->setMinimumObstacleClearance(0.365);

      return ss;
    }

    const ompl::base::State* getGeometricComponentStateInternal(const ompl::base::State *state, unsigned int /*index*/) const {
      return state;
    }

    void loadGoals() {
      using namespace arma;
      // Load XML containing landmarks
      TiXmlDocument doc(pathToSetupFile_);
      bool loadOkay = doc.LoadFile();

      if(!loadOkay) {
        printf( "Could not load setup file. Error='%s'. Exiting.\n", doc.ErrorDesc() );
        exit( 1 );
      }

      TiXmlNode* node = 0;
      TiXmlElement* landmarkElement = 0;
      TiXmlElement* itemElement = 0;

      // Get the landmarklist node
      node = doc.FirstChild( "GoalList" );
      assert( node );
      landmarkElement = node->ToElement(); //convert node to element
      assert( landmarkElement  );

      TiXmlNode* child = 0;

      //Iterate through all the landmarks and put them into the "landmarks_" list
      while( (child = landmarkElement ->IterateChildren(child))) {
        assert( child );
        itemElement = child->ToElement();
        assert( itemElement );

        double goalX = 0 , goalY = 0, goalTheta = 0;

        itemElement->QueryDoubleAttribute("x", &goalX);
        itemElement->QueryDoubleAttribute("y", &goalY);

        std::cout<<"Loaded Goal Pose X: "<<goalX<<" Y: "<<goalY<<std::endl;

        addGoalState(goalX, goalY);
      }
    }

    void loadDynamicObstaclesList() {
      // Load XML containing dynamic obstacle model paths
      TiXmlDocument doc(pathToSetupFile_);
      bool loadOkay = doc.LoadFile();

      if(!loadOkay) {
        printf( "Could not load setup file. Error='%s'. Exiting.\n", doc.ErrorDesc() );
        exit( 1 );
      }

      TiXmlNode* node = 0;
      TiXmlElement* landmarkElement = 0;
      TiXmlElement* itemElement = 0;

      // Get the dynamic obstacle list node
      node = doc.FirstChild( "DynamicObstacleList" );
      assert( node );

      landmarkElement = node->ToElement(); //convert node to element
      assert( landmarkElement  );

      TiXmlNode* child = 0;

      //Iterate through all the obstacles and put them into a list
      while( (child = landmarkElement ->IterateChildren(child))) {
        assert( child );
        itemElement = child->ToElement();
        assert( itemElement );

        double id = 0;

        string modelPath;

        itemElement->QueryDoubleAttribute("id", &id);
        itemElement->QueryStringAttribute("path", &modelPath);

        std::cout<<"Loaded Dynamic Obstacle: "<<modelPath<<std::endl;

        dynObstList_.push_back(modelPath);
      }
    }

  void loadParameters() {
    using namespace arma;

    TiXmlDocument doc(pathToSetupFile_);

    bool loadOkay = doc.LoadFile();

    if ( !loadOkay ) {
      printf( "Could not load setup file in planning problem. Error='%s'. Exiting.\n", doc.ErrorDesc() );
      exit( 1 );
    }

    TiXmlNode* node = 0;

    TiXmlElement* itemElement = 0;

    node = doc.FirstChild( "PlanningProblem" );
    assert( node );

    TiXmlNode* child = 0;

    // Planner Mode
    child = node->FirstChild("PlannerMode");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    // planner method
    int methodChoice = 0;
    itemElement->QueryIntAttribute("method", &methodChoice);
    plannerMethod_ = methodChoice;

    //Dynamic Obstacles
    int dynobst = 0;
    itemElement->QueryIntAttribute("dynobst", &dynobst);

    if(dynobst == 1) {
      dynamicObstacles_ = true;
      loadDynamicObstaclesList();
    } else {
      dynamicObstacles_ = false;
    }

    // Read the env mesh file
    child = node->FirstChild("Environment");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    std::string environmentFilePath;
    itemElement->QueryStringAttribute("environmentFile", &environmentFilePath);
    pathToEnvironmentMesh_ = environmentFilePath;
    this->addEnvironmentMesh(environmentFilePath);

    // Read the robot mesh file
    child  = node->FirstChild("Robot");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    std::string robotFilePath;
    itemElement->QueryStringAttribute("robotFile", &robotFilePath);

    this->setRobotMesh(robotFilePath);
   
    // Read the roadmap filename
    child  = node->FirstChild("RoadMap");
    assert( child );
    itemElement = child->ToElement();
    assert( itemElement );

    std::string tempPathStr;
    itemElement->QueryStringAttribute("roadmapFile", &tempPathStr);
    pathToRoadMapFile_ = tempPathStr;

    int usermap = 0;
    itemElement->QueryIntAttribute("useRoadMap", &usermap);
    useSavedRoadMap_ = usermap;

    // Read the start Pose
    child  = node->FirstChild("StartPose");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    double startX = 0,startY = 0;

    itemElement->QueryDoubleAttribute("x", &startX);
    itemElement->QueryDoubleAttribute("y", &startY);

    setStartState(startX, startY);

    // Read the Goal Pose
    /*
    child  = node->FirstChild("GoalPose");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    double goalX = 0 , goalY = 0, goalTheta = 0;

    itemElement->QueryDoubleAttribute("x", &goalX);
    itemElement->QueryDoubleAttribute("y", &goalY);
    itemElement->QueryDoubleAttribute("theta", &goalTheta);

    setGoalState(goalX, goalY, goalTheta);
    */

    // read planning time
    child  = node->FirstChild("PlanningTime");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    double time = 0;

    itemElement->QueryDoubleAttribute("maxTime", &time) ;

    planningTime_ = time;

    // read planning time
    child  = node->FirstChild("FIRMNodes");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    int minNodeNum = 0;
    itemElement->QueryIntAttribute("minNodes", &minNodeNum) ;
    minNodes_ = minNodeNum;

    int maxNodeNum = 0;
    itemElement->QueryIntAttribute("maxNodes", &maxNodeNum) ;
    maxNodes_ = maxNodeNum;

    // Read Kidnapped State
    // Read the Goal Pose
    child  = node->FirstChild("KidnappedState");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    double kX = 0 , kY = 0;

    itemElement->QueryDoubleAttribute("x", &kX);
    itemElement->QueryDoubleAttribute("y", &kY);

    kidnappedState_ = siF_->allocState();

    kidnappedState_->as<SE2BeliefSpace::StateType>()->setXY(kX, kY);

    loadGoals();

    OMPL_INFORM("Problem configuration is");
    std::cout<<"Path to environment mesh: "<<environmentFilePath<<std::endl;
    std::cout<<"Path to robot mesh: "<<robotFilePath<<std::endl;
    std::cout<<"Path to Roadmap File: "<<pathToRoadMapFile_<<std::endl;
    std::cout<<"Start Pose X: "<<startX<<" Y: "<<startY<<std::endl;
    std::cout<<"Planning Time: "<<planningTime_<<" seconds"<<std::endl;
    std::cout<<"Min Nodes: "<<minNodes_<<std::endl;
    std::cout<<"Max Nodes: "<<maxNodes_<<std::endl;
    std::cout<<"Kidnapped Pose x:"<<kX<<" y:"<<kY<<std::endl;

    // Goal Constraints
    node = doc.FirstChild( "GoalConstraints" );
    assert( node );

    // Weights
    child = 0;
    child = node->FirstChild("Weights");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    double meanNormWeight = 0.0;
    itemElement->QueryDoubleAttribute("meanNormWeight", &meanNormWeight);
    StateType::meanNormWeight_ = meanNormWeight;

    double covNormWeight = 0.0;
    itemElement->QueryDoubleAttribute("covNormWeight", &covNormWeight);
    StateType::covNormWeight_ = covNormWeight;

    // Threshold
    child = 0;
    child = node->FirstChild("Threshold");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    double reachDist = 0.0;
    itemElement->QueryDoubleAttribute("reachDist", &reachDist);
    StateType::reachDist_ = reachDist;

    OMPL_INFORM("Goal Constraints are");
    std::cout<<"meanNormWeight: "<<StateType::meanNormWeight_<<std::endl;
    std::cout<<"covNormWeight: "<<StateType::covNormWeight_<<std::endl;
    std::cout<<"reachDist: "<<StateType::reachDist_<<std::endl;
  }

protected:
  ompl::base::State *start_;

  std::vector<ompl::base::State*> goalList_;

  ompl::base::State *kidnappedState_;

  firm::SpaceInformation::SpaceInformationPtr siF_;

  ompl::control::StatePropagatorPtr statePropagator_;

  ompl::control::ControlSpacePtr cs_;

  ompl::base::StateSpacePtr ss_;

  ompl::base::PlannerPtr planner_;

  ompl::base::ProblemDefinitionPtr pdef_;

  ompl::base::StateValidityCheckerPtr vc_;

  std::string pathToSetupFile_;

  std::string pathToRoadMapFile_;

  std::string pathToEnvironmentMesh_;

  int useSavedRoadMap_;

  double planningTime_;

  unsigned int minNodes_;
  unsigned int maxNodes_;

  bool setup_;

  /** \brief Activate dynamic obstacles */
  bool dynamicObstacles_;

  std::vector<string> dynObstList_;

  int plannerMethod_;
};
#endif
