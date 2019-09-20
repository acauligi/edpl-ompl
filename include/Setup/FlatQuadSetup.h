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
    
// class myMotionValidator : public base::MotionValidator {
//   public:
//     bool checkMotion(const State *s1, const State *s2) const {
//       // TODO(acauligi)
//       return true;
//     }
// };

/** \brief Wrapper for ompl::app::RigidBodyPlanning that plans for rigid bodies in FlatQuadBeliefSpace for a MAV using FIRM */
class FlatQuadSetup : public ompl::app::RigidBodyGeometry {
    typedef FlatQuadBeliefSpace::StateType StateType;

public:
    FlatQuadSetup() : ompl::app::RigidBodyGeometry(ompl::app::Motion_3D, ompl::app::FCL),
    ss_(ompl::base::StateSpacePtr(new FlatQuadBeliefSpace())) {
      // set static variables
      // Not in use as of now
      Controller<FiniteTimeLQR, ExtendedKF>::setNodeReachedAngle(10); // degrees
      Controller<FiniteTimeLQR, ExtendedKF>::setNodeReachedDistance(0.1);// meters
      Controller<FiniteTimeLQR, ExtendedKF>::setMaxTries(30);
      Controller<FiniteTimeLQR, ExtendedKF>::setMaxTrajectoryDeviation(0.5); // meters

      // [2] NodeController (typedef in include/Planner/FIRM.h)
      Controller<StationaryLQR, LinearizedKF>::setNodeReachedAngle(10); // degrees
      Controller<StationaryLQR, LinearizedKF>::setNodeReachedDistance(0.1);// meters
      Controller<StationaryLQR, LinearizedKF>::setMaxTries(300);
      Controller<StationaryLQR, LinearizedKF>::setMaxTrajectoryDeviation(2.50); // meters

      // [1] EdgeController (typedef in include/Planner/FIRM.h)
      RHCICreate::setControlQueueSize(5);
      RHCICreate::setTurnOnlyDistance(0.01);
      Controller<RHCICreate, ExtendedKF>::setNodeReachedAngle(10.0); // degrees
      Controller<RHCICreate, ExtendedKF>::setNodeReachedDistance(0.1);// meters
      Controller<RHCICreate, ExtendedKF>::setMaxTries(30);
      Controller<RHCICreate, ExtendedKF>::setMaxTrajectoryDeviation(1.0); // meters

      // setting the mean and norm weights (used in reachability check)
      // NOTE these values will be overwritten by loadParameters()
      StateType::covNormWeight_ =  1.0;
      StateType::meanNormWeight_=  2.0;
      StateType::reachDist_ =  0.009;    // distance threshold for position, orientation, and covariance
      StateType::reachDistPos_ = 0.1;    // distance threshold for position [m]
      StateType::reachDistOri_ = 10.0/180.0*boost::math::constants::pi<double>(); // distance threshold for orientation [rad]
      StateType::reachDistCov_ = 0.0004;    // distance threshold for covariance

      // set the state component norm weights
      arma::colvec normWeights(12);
      normWeights(0)  = 2.0/std::sqrt(9);
      normWeights(1)  = 2.0/std::sqrt(9);
      normWeights(2)  = 2.0/std::sqrt(9);
      normWeights(3)  = 1.0/std::sqrt(9);
      normWeights(4)  = 1.0/std::sqrt(9);
      normWeights(5)  = 1.0/std::sqrt(9);
      normWeights(6)  = 1.0/std::sqrt(9);
      normWeights(7)  = 1.0/std::sqrt(9);
      normWeights(8)  = 1.0/std::sqrt(9);
      normWeights(9)  = 1.0/std::sqrt(9);
      normWeights(10) = 1.0/std::sqrt(9);
      normWeights(11) = 1.0/std::sqrt(9);
      StateType::normWeights_ = normWeights;

      // The bounds should be inferred from the geometry files,
      // there is a function in Apputils to do this, so use that.        
      ompl::base::RealVectorBounds bounds(12);
      bounds.setLow(0, -20.);     bounds.setHigh(0, 20.);
      bounds.setLow(1, -20.);     bounds.setHigh(1, 20.);
      bounds.setLow(2, -20.);     bounds.setHigh(2, 20.);
      bounds.setLow(3, -20.);     bounds.setHigh(3, 20.);    // yaw angle value not used

      bounds.setLow(4, -0.0001);     bounds.setHigh(4, 0.0);
      bounds.setLow(5, -0.0001);     bounds.setHigh(5, 0.0);
      bounds.setLow(6, -0.0001);     bounds.setHigh(6, 0.0);
      bounds.setLow(7, -0.0001);     bounds.setHigh(7, 0.0);

      bounds.setLow(8, -0.0001);     bounds.setHigh(8, 0.);
      bounds.setLow(9, -0.0001);     bounds.setHigh(9, 0.);
      bounds.setLow(10, -0.0001);     bounds.setHigh(10, 0.);
      bounds.setLow(11, -0.0001);     bounds.setHigh(11, 0.);

      ss_->as<FlatQuadBeliefSpace>()->setBounds(bounds);

      //Construct the control space
      ompl::control::ControlSpacePtr controlspace( new ompl::control::RealVectorControlSpace(ss_,4) );

      cs_ = controlspace;

      // construct an instance of space information from this state space
      firm::SpaceInformation::SpaceInformationPtr si(new firm::SpaceInformation(ss_, cs_));
      siF_ = si;

      //siF_->setValidStateSamplerAllocator(FlatQuadSetup::allocMaxClearanceValidStateSampler);

      ompl::base::ProblemDefinitionPtr prblm(new ompl::base::ProblemDefinition(siF_));

      pdef_ = prblm;

      start_ = siF_->allocState();

      setup_ = false;

      dynamic_obstacles_ = false;

      planner_method_ = 0; // by default we use FIRM
    }

    virtual ~FlatQuadSetup(void) {
    }

    const firm::SpaceInformation::SpaceInformationPtr& getSpaceInformation() const {
      return siF_;
    }

    void setPathToSetupFile(const std::string &path) {
      path_to_setup_file_  = path;       
    }

    void setStartState(const double X, const double Y, const double Z, const double Yaw) {
      ompl::base::State *temp = siF_->allocState();
      temp->as<StateType>()->setXYZYaw(X,Y,Z,Yaw);

      // Set velocity and acceleration to zero
      arma::colvec rest_state = arma::zeros<arma::colvec>(4);
      temp->as<StateType>()->setVelocity(rest_state);
      temp->as<StateType>()->setAcceleration(rest_state);

      siF_->copyState(start_, temp);
      siF_->freeState(temp);
    }

    void addGoalState(const double X, const double Y, const double Z, const double Yaw) {
      ompl::base::State *temp = siF_->allocState();
      temp->as<StateType>()->setXYZYaw(X,Y,Z,Yaw);

      // Set velocity and acceleration to zero
      arma::colvec rest_state = arma::zeros<arma::colvec>(4);
      temp->as<StateType>()->setVelocity(rest_state);
      temp->as<StateType>()->setAcceleration(rest_state);

      goal_list_.push_back(temp);
    }

    virtual void setup() {
      if(!setup_) {
        this->loadParameters();

        if(path_to_setup_file_.length() == 0) {
          throw ompl::Exception("Path to setup file not set!");
        }

        if(!hasEnvironment() || !hasRobot()) {
          throw ompl::Exception("Robot/Environment mesh files not setup!");
        }

        // TODO(acauligi): how to set bounds correctly given new environment mesh?
        // ss_->as<FlatQuadBeliefSpace>()->setBounds(inferEnvironmentBounds());

        // Create an FCL state validity checker and assign to space information
        const ompl::base::StateValidityCheckerPtr &fclSVC = this->allocStateValidityChecker(siF_, getGeometricStateExtractor(), false);
        siF_->setStateValidityChecker(fclSVC);

        // siF_->setMotionValidator(std::make_shared<myMotionValidator>(siF_));

        // provide the observation model to the space
        ObservationModelMethod::ObservationModelPointer om(new LandmarkObservationModel(siF_, path_to_setup_file_.c_str()));
        siF_->setObservationModel(om);

        // Provide the motion model to the space
        MotionModelMethod::MotionModelPointer mm(new FlatQuadMotionModel(siF_, path_to_setup_file_.c_str()));            
        siF_->setMotionModel(mm);

        ompl::control::StatePropagatorPtr prop(ompl::control::StatePropagatorPtr(new FlatQuadStatePropagator(siF_)));
        state_propagator_ = prop;
        siF_->setStatePropagator(state_propagator_);    // func that performs state propagation
        siF_->setPropagationStepSize(0.1); // this is the duration that a control is applied
        siF_->setStateValidityCheckingResolution(0.005);  // controls are applied a time duration that is an integer multiple of argument here
        siF_->setMinMaxControlDuration(1,100);  // min+max number of steps a control is propagated for

        if(!start_ || goal_list_.size() == 0) {
          throw ompl::Exception("Start/Goal not set");
        }

        pdef_->setStartAndGoalStates(start_, goal_list_[0], 1.0);

        // Instantiate new planner object
        ompl::base::PlannerPtr plnr(new FIRM(siF_, false));
        planner_ = plnr;
        planner_->setProblemDefinition(pdef_);
        planner_->as<FIRM>()->setMinFIRMNodes(min_nodes_);
        planner_->as<FIRM>()->setMaxFIRMNodes(max_nodes_);
        planner_->as<FIRM>()->setKidnappedState(kidnapped_state_);
        planner_->as<FIRM>()->loadParametersFromFile(path_to_setup_file_.c_str());
        planner_->setup();            

        // Setup visualizer because it is needed while loading roadmap and visualizing it
        Visualizer::updateSpaceInformation(this->getSpaceInformation());
        Visualizer::updateRenderer(*dynamic_cast<const ompl::app::RigidBodyGeometry*>(this), this->getGeometricStateExtractor());

        if (use_saved_road_map_ == 1) {
          planner_->as<FIRM>()->loadRoadMapFromFile(path_to_road_map_file_.c_str());
        }

        setup_ = true;
      }
    }

    ompl::base::PlannerStatus solve() {
      if(!setup_) {
        this->setup();
      }
      return planner_->solve(planning_time_);
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
      executeSolution(planner_method_);

      // Need a function to get terminated state
      if(goal_list_.size() > 1) {
        updateEnvironmentMesh();

        for(int ii=0; ii < goal_list_.size()-1;ii++) {
          pdef_->setStartAndGoalStates(goal_list_[ii], goal_list_[ii+1], 1.0);

          planner_->setProblemDefinition(pdef_);

          if(this->solve()) {
            executeSolution(planner_method_);
          }
        }
      }
    }

    virtual void updateEnvironmentMesh(int obindx = 0) {
      if(dynamic_obstacles_) {
        // Set environment to new mesh with some dynamic / additional obstacles
        if(!this->setEnvironmentMesh(dyn_obst_list_[obindx])) {
          OMPL_ERROR("Couldn't set mesh with path: %s",dyn_obst_list_[obindx]);
        }
        
        const ompl::base::StateValidityCheckerPtr &svc = std::make_shared<ompl::app::FCLStateValidityChecker<ompl::app::Motion_3D>>(siF_,  getGeometrySpecification(), getGeometricStateExtractor(), false);

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
      // Convert FlatQuadBeliefSpace representation to SE3StateSpace pose used for collision checking
      ompl::base::StateSpacePtr pose_space(new ompl::base::SE3StateSpace());
      ompl::base::State * pose_state = pose_space->allocState();

      pose_state->as<ompl::base::SE3StateSpace::StateType>()->setXYZ(
        state->as<FlatQuadBeliefSpace::StateType>()->getX(),
        state->as<FlatQuadBeliefSpace::StateType>()->getY(),
        state->as<FlatQuadBeliefSpace::StateType>()->getZ()
      );

      arma::colvec axis_angle = arma::zeros<arma::colvec>(4);
      axis_angle = state->as<FlatQuadBeliefSpace::StateType>()->flatToAxisAngle(); 
      pose_state->as<ompl::base::SE3StateSpace::StateType>()->rotation().setAxisAngle(axis_angle[0], axis_angle[1], axis_angle[2], axis_angle[3]);

      return pose_state;
    }

    void loadGoals() {
      using namespace arma;
      // Load XML containing landmarks
      TiXmlDocument doc(path_to_setup_file_);
      bool loadOkay = doc.LoadFile();

      if(!loadOkay) {
        printf( "Could not load setup file. Error='%s'. Exiting.\n", doc.ErrorDesc() );
        exit( 1 );
      }

      TiXmlNode* node = 0;
      TiXmlElement* goalElement = 0;
      TiXmlElement* itemElement = 0;

      // Get the goallist node
      node = doc.FirstChild( "GoalList" );
      assert( node );
      goalElement = node->ToElement(); //convert node to element
      assert( goalElement  );

      TiXmlNode* child = 0;

      //Iterate through all the goals and assign last goal as goal state 
      while( (child = goalElement ->IterateChildren(child))) {
        assert( child );
        itemElement = child->ToElement();
        assert( itemElement );

        double goalX = 0 , goalY = 0, goalZ = 0, goalYaw = 0;

        itemElement->QueryDoubleAttribute("x", &goalX);
        itemElement->QueryDoubleAttribute("y", &goalY);
        itemElement->QueryDoubleAttribute("z", &goalZ);
        itemElement->QueryDoubleAttribute("yaw", &goalYaw);

        std::cout<<"Loaded Goal Pose X: "<< goalX << 
          " Y: "<< goalY  << 
          " Z: "<< goalZ <<
          " Yaw: "<< goalYaw<<
          std::endl;
        addGoalState(goalX, goalY, goalZ, goalYaw);
      }
    }

    void loadDynamicObstaclesList() {
      // Load XML containing dynamic obstacle model paths
      TiXmlDocument doc(path_to_setup_file_);
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

        dyn_obst_list_.push_back(modelPath);
      }
    }

  void loadParameters() {
    using namespace arma;

    TiXmlDocument doc(path_to_setup_file_);

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
    planner_method_ = methodChoice;

    //Dynamic Obstacles
    int dynobst = 0;
    itemElement->QueryIntAttribute("dynobst", &dynobst);

    if(dynobst == 1) {
      dynamic_obstacles_ = true;
      loadDynamicObstaclesList();
    } else {
      dynamic_obstacles_ = false;
    }

    // Read the env mesh file
    child = node->FirstChild("Environment");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    std::string environment_file_path;
    itemElement->QueryStringAttribute("environmentFile", &environment_file_path);
    path_to_environment_mesh_ = environment_file_path;
    this->addEnvironmentMesh(environment_file_path);

    // Read the robot mesh file
    child  = node->FirstChild("Robot");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    std::string robot_file_path;
    itemElement->QueryStringAttribute("robotFile", &robot_file_path);
    this->setRobotMesh(robot_file_path);
   
    // Read the roadmap filename
    child  = node->FirstChild("RoadMap");
    assert( child );
    itemElement = child->ToElement();
    assert( itemElement );

    std::string temp_path_str;
    itemElement->QueryStringAttribute("roadmapFile", &temp_path_str);
    path_to_road_map_file_ = temp_path_str;

    int usermap = 0;
    itemElement->QueryIntAttribute("useRoadMap", &usermap);
    use_saved_road_map_ = usermap;

    // Read the start Pose
    child  = node->FirstChild("StartPose");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    double startX = 0.,startY = 0., startZ = 0., startYaw = 0.;

    itemElement->QueryDoubleAttribute("x", &startX);
    itemElement->QueryDoubleAttribute("y", &startY);
    itemElement->QueryDoubleAttribute("z", &startZ);
    itemElement->QueryDoubleAttribute("yaw", &startYaw);
    setStartState(startX, startY, startZ, startYaw);
    std::cout<<"Loaded Start Pose X: "<< startX << 
      " Y: "<< startY << 
      " Z: "<< startZ <<
      " Yaw: "<< startYaw <<
      std::endl;

    // read planning time
    child  = node->FirstChild("PlanningTime");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    double time = 0;
    itemElement->QueryDoubleAttribute("maxTime", &time) ;
    planning_time_ = time;

    // read planning time
    child  = node->FirstChild("FIRMNodes");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    int min_node_num = 0;
    itemElement->QueryIntAttribute("minNodes", &min_node_num) ;
    min_nodes_ = min_node_num;

    int max_node_num = 0;
    itemElement->QueryIntAttribute("maxNodes", &max_node_num) ;
    max_nodes_ = max_node_num;

    // Read Kidnapped State
    // Read the Goal Pose
    child  = node->FirstChild("KidnappedState");
    assert( child );

    itemElement = child->ToElement();
    assert( itemElement );

    double kX = 0 , kY = 0, kZ = 0, kYaw = 0;
    itemElement->QueryDoubleAttribute("x", &kX);
    itemElement->QueryDoubleAttribute("y", &kY);
    itemElement->QueryDoubleAttribute("z", &kZ);
    itemElement->QueryDoubleAttribute("yaw", &kYaw);

    kidnapped_state_ = siF_->allocState();
    kidnapped_state_->as<FlatQuadBeliefSpace::StateType>()->setXYZYaw(kX, kY, kZ, kYaw);

    loadGoals();

    OMPL_INFORM("Problem configuration is");
    std::cout<<"Path to environment mesh: "<<environment_file_path<<std::endl;
    std::cout<<"Path to robot mesh: "<<robot_file_path<<std::endl;
    std::cout<<"Path to Roadmap File: "<<path_to_road_map_file_<<std::endl;
    std::cout<<"Start Pose X: " << startX << " Y: " << startY << " Z: " << startZ << std::endl;
    std::cout<<"Planning Time: "<<planning_time_<<" seconds"<<std::endl;
    std::cout<<"Min Nodes: "<<min_nodes_<<std::endl;
    std::cout<<"Max Nodes: "<<max_nodes_<<std::endl;
    std::cout<<"Kidnapped Pose x:" << kX << " y:" << kY << " z: " << kZ <<std::endl;

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

  std::vector<ompl::base::State*> goal_list_;

  ompl::base::State *kidnapped_state_;

  firm::SpaceInformation::SpaceInformationPtr siF_;

  ompl::control::StatePropagatorPtr state_propagator_;

  ompl::control::ControlSpacePtr cs_;

  ompl::base::StateSpacePtr ss_;

  ompl::base::PlannerPtr planner_;

  ompl::base::ProblemDefinitionPtr pdef_;

  ompl::base::StateValidityCheckerPtr vc_;

  std::string path_to_setup_file_;

  std::string path_to_road_map_file_;

  std::string path_to_environment_mesh_;

  int use_saved_road_map_;

  double planning_time_;

  unsigned int min_nodes_;
  unsigned int max_nodes_;

  bool setup_;

  /** \brief Activate dynamic obstacles */
  bool dynamic_obstacles_;

  std::vector<string> dyn_obst_list_;

  int planner_method_;
};
#endif
