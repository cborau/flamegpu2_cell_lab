FLAMEGPU_AGENT_FUNCTION(ecm_boundary_concentration_conditions, flamegpu::MessageNone, flamegpu::MessageNone) {
  // Agent properties in local register
  int id = FLAMEGPU->getVariable<int>("id");

  // Agent position
  float agent_x = FLAMEGPU->getVariable<float>("x");
  float agent_y = FLAMEGPU->getVariable<float>("y");
  float agent_z = FLAMEGPU->getVariable<float>("z");
    
  // Agent concentration
  uint8_t N_SPECIES = FLAMEGPU->getVariable<uint8_t>("N_SPECIES");
  float agent_conc = FLAMEGPU->getVariable<float>("concentration"); 
 
  float separation_x_pos = 0.0;
  float separation_x_neg = 0.0;
  float separation_y_pos = 0.0;
  float separation_y_neg = 0.0;
  float separation_z_pos = 0.0;
  float separation_z_neg = 0.0;
  const float ECM_BOUNDARY_INTERACTION_RADIUS = FLAMEGPU->environment.getProperty<float>("ECM_BOUNDARY_INTERACTION_RADIUS");
  const float ECM_BOUNDARY_EQUILIBRIUM_DISTANCE = FLAMEGPU->environment.getProperty<float>("ECM_BOUNDARY_EQUILIBRIUM_DISTANCE");
  float EPSILON = FLAMEGPU->environment.getProperty<float>("EPSILON");

  // Get position of the boundaries
  const float COORD_BOUNDARY_X_POS = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",0);
  const float COORD_BOUNDARY_X_NEG = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",1);
  const float COORD_BOUNDARY_Y_POS = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",2);
  const float COORD_BOUNDARY_Y_NEG = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",3);
  const float COORD_BOUNDARY_Z_POS = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",4);
  const float COORD_BOUNDARY_Z_NEG = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",5);
  
  const float BOUNDARY_CONC_INIT_X_POS = FLAMEGPU->environment.getProperty<float>("BOUNDARY_CONC_INIT", 0);
  const float BOUNDARY_CONC_INIT_X_NEG = FLAMEGPU->environment.getProperty<float>("BOUNDARY_CONC_INIT", 1);
  const float BOUNDARY_CONC_INIT_Y_POS = FLAMEGPU->environment.getProperty<float>("BOUNDARY_CONC_INIT", 2);
  const float BOUNDARY_CONC_INIT_Y_NEG = FLAMEGPU->environment.getProperty<float>("BOUNDARY_CONC_INIT", 3);
  const float BOUNDARY_CONC_INIT_Z_POS = FLAMEGPU->environment.getProperty<float>("BOUNDARY_CONC_INIT", 4);
  const float BOUNDARY_CONC_INIT_Z_NEG = FLAMEGPU->environment.getProperty<float>("BOUNDARY_CONC_INIT", 5);

  const float BOUNDARY_CONC_FIXED_X_POS = FLAMEGPU->environment.getProperty<float>("BOUNDARY_CONC_FIXED", 0);
  const float BOUNDARY_CONC_FIXED_X_NEG = FLAMEGPU->environment.getProperty<float>("BOUNDARY_CONC_FIXED", 1);
  const float BOUNDARY_CONC_FIXED_Y_POS = FLAMEGPU->environment.getProperty<float>("BOUNDARY_CONC_FIXED", 2);
  const float BOUNDARY_CONC_FIXED_Y_NEG = FLAMEGPU->environment.getProperty<float>("BOUNDARY_CONC_FIXED", 3);
  const float BOUNDARY_CONC_FIXED_Z_POS = FLAMEGPU->environment.getProperty<float>("BOUNDARY_CONC_FIXED", 4);
  const float BOUNDARY_CONC_FIXED_Z_NEG = FLAMEGPU->environment.getProperty<float>("BOUNDARY_CONC_FIXED", 5);

  // Check for ecm-boundary separations
  //this takes into account the distance with respect to the actual boundary position, while forces are calculated with respect to boundary initial position
  separation_x_pos = (agent_x - COORD_BOUNDARY_X_POS); 
  separation_x_neg = (agent_x - COORD_BOUNDARY_X_NEG);
  separation_y_pos = (agent_y - COORD_BOUNDARY_Y_POS);
  separation_y_neg = (agent_y - COORD_BOUNDARY_Y_NEG);
  separation_z_pos = (agent_z - COORD_BOUNDARY_Z_POS);
  separation_z_neg = (agent_z - COORD_BOUNDARY_Z_NEG);

  // TODO: CHECK ELEMENTS TOUCHING MULTIPLE BOUNDARIES
  // Check concentration conditions
  if (fabsf(separation_x_pos) < (ECM_BOUNDARY_INTERACTION_RADIUS)){
	  if (BOUNDARY_CONC_FIXED_X_POS >= 0.0){
		  agent_conc = BOUNDARY_CONC_FIXED_X_POS;
	  }
	  if (BOUNDARY_CONC_INIT_X_POS >= 0.0){
		  agent_conc = BOUNDARY_CONC_INIT_X_POS;
	  }
  }
  if (fabsf(separation_x_neg) < (ECM_BOUNDARY_INTERACTION_RADIUS)){
	  if (BOUNDARY_CONC_FIXED_X_NEG >= 0.0){
		  agent_conc = BOUNDARY_CONC_FIXED_X_NEG;
	  }
	  if (BOUNDARY_CONC_INIT_X_NEG >= 0.0){
		  agent_conc = BOUNDARY_CONC_INIT_X_NEG;
	  }
  }
  if (fabsf(separation_y_pos) < (ECM_BOUNDARY_INTERACTION_RADIUS)){
	  if (BOUNDARY_CONC_FIXED_Y_POS >= 0.0){
		  agent_conc = BOUNDARY_CONC_FIXED_Y_POS;
	  }
	  if (BOUNDARY_CONC_INIT_Y_POS >= 0.0){
		  agent_conc = BOUNDARY_CONC_INIT_Y_POS;
	  }
  }
  if (fabsf(separation_y_neg) < (ECM_BOUNDARY_INTERACTION_RADIUS)){
	  if (BOUNDARY_CONC_FIXED_Y_NEG >= 0.0){
		  agent_conc = BOUNDARY_CONC_FIXED_Y_NEG;
	  }
	  if (BOUNDARY_CONC_INIT_Y_NEG >= 0.0){
		  agent_conc = BOUNDARY_CONC_INIT_Y_NEG;
	  }
  }
  if (fabsf(separation_z_pos) < (ECM_BOUNDARY_INTERACTION_RADIUS)){
	  if (BOUNDARY_CONC_FIXED_Z_POS >= 0.0){
		  agent_conc = BOUNDARY_CONC_FIXED_Z_POS;
	  }
	  if (BOUNDARY_CONC_INIT_Z_POS >= 0.0){
		  agent_conc = BOUNDARY_CONC_INIT_Z_POS;
	  }
  }
  if (fabsf(separation_z_neg) < (ECM_BOUNDARY_INTERACTION_RADIUS)){
	  if (BOUNDARY_CONC_FIXED_Z_NEG >= 0.0){
		  agent_conc = BOUNDARY_CONC_FIXED_Z_NEG;
	  }
	  if (BOUNDARY_CONC_INIT_Z_NEG >= 0.0){
		  agent_conc = BOUNDARY_CONC_INIT_Z_NEG;
	  }
  }
  

  FLAMEGPU->setVariable<float>("concentration", agent_conc);

  return flamegpu::ALIVE;
}