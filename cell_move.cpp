FLAMEGPU_DEVICE_FUNCTION float vec3Length(const float x, const float y, const float z) {
  return sqrtf(x * x + y * y + z * z);
}

FLAMEGPU_AGENT_FUNCTION(cell_move, flamegpu::MessageArray3D, flamegpu::MessageNone) {
  //Agent position vector
  int id = FLAMEGPU->getVariable<int>("id");
  float agent_x = FLAMEGPU->getVariable<float>("x");
  float agent_y = FLAMEGPU->getVariable<float>("y");
  float agent_z = FLAMEGPU->getVariable<float>("z");
  float agent_vx = FLAMEGPU->getVariable<float>("vx");
  float agent_vy = FLAMEGPU->getVariable<float>("vy");
  float agent_vz = FLAMEGPU->getVariable<float>("vz");
  float agent_orx = FLAMEGPU->getVariable<float>("orx");
  float agent_ory = FLAMEGPU->getVariable<float>("ory");
  float agent_orz = FLAMEGPU->getVariable<float>("orz");
  float agent_k_elast = FLAMEGPU->getVariable<float>("k_elast");
  
  // Get number of agents per direction
  const int Nx = FLAMEGPU->environment.getProperty<int>("ECM_AGENTS_PER_DIR",0);
  const int Ny = FLAMEGPU->environment.getProperty<int>("ECM_AGENTS_PER_DIR",1);
  const int Nz = FLAMEGPU->environment.getProperty<int>("ECM_AGENTS_PER_DIR",2);
  // Get position of the boundaries
  const float COORD_BOUNDARY_X_POS = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",0);
  const float COORD_BOUNDARY_X_NEG = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",1);
  const float COORD_BOUNDARY_Y_POS = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",2);
  const float COORD_BOUNDARY_Y_NEG = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",3);
  const float COORD_BOUNDARY_Z_POS = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",4);
  const float COORD_BOUNDARY_Z_NEG = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",5);
  
  // transform x,y,z positions to i,j,k grid positions
  int agent_grid_i = roundf(((agent_x - COORD_BOUNDARY_X_NEG) / (COORD_BOUNDARY_X_POS - COORD_BOUNDARY_X_NEG)) * (Nx - 1));
  int agent_grid_j = roundf(((agent_y - COORD_BOUNDARY_Y_NEG) / (COORD_BOUNDARY_Y_POS - COORD_BOUNDARY_Y_NEG)) * (Ny - 1));
  int agent_grid_k = roundf(((agent_z - COORD_BOUNDARY_Z_NEG) / (COORD_BOUNDARY_Z_POS - COORD_BOUNDARY_Z_NEG)) * (Nz - 1));
  
  // direction: the vector joining interacting agents
  float dir_x = 0.0; 
  float dir_y = 0.0; 
  float dir_z = 0.0; 
  float distance = 0.0;
  const float ECM_ECM_EQUILIBRIUM_DISTANCE = FLAMEGPU->environment.getProperty<float>("ECM_ECM_EQUILIBRIUM_DISTANCE");
  const float DELTA_TIME = FLAMEGPU->environment.getProperty<float>("DELTA_TIME");
  float min_distance = ECM_ECM_EQUILIBRIUM_DISTANCE; //initialize to maximum possible value between ECM agents
  
  int message_id = 0;
  float message_x = 0.0;
  float message_y = 0.0;
  float message_z = 0.0;
  float message_vx = 0.0;
  float message_vy = 0.0;
  float message_vz = 0.0;
  
  for (const auto &message : FLAMEGPU->message_in(agent_grid_i, agent_grid_j, agent_grid_k)) { // find the closest ECM agent and move with it
	message_id = message.getVariable<int>("id");
	message_x = message.getVariable<float>("x");
    message_y = message.getVariable<float>("y");
    message_z = message.getVariable<float>("z");  
	message_vx = message.getVariable<float>("vx");
    message_vy = message.getVariable<float>("vy");
    message_vz = message.getVariable<float>("vz");
	
	dir_x = agent_x - message_x; 
    dir_y = agent_y - message_y; 
    dir_z = agent_z - message_z; 
    distance = vec3Length(dir_x, dir_y, dir_z); 

	//printf("MESSAGE ID: %d, dist: %g, pos -> (%g, %g, %g), vel -> (%g, %g, %g)\n", message_id, distance, message_x, message_y, message_z, message_vx, message_vy, message_vz);
	
               
    if (distance < min_distance) {	
		min_distance = distance;
		//agent_vx = message_vx;
		//agent_vy = message_vy;
		//agent_vz = message_vz;
		agent_vx = 0.0;
		agent_vy = 0.0;
		agent_vz = 0.0; // TODO: cell migration and orientation aligment with fibers
    }  
  }
 
  agent_x += agent_vx * DELTA_TIME;
  agent_y += agent_vy * DELTA_TIME;
  agent_z += agent_vz * DELTA_TIME;
  
  //printf("VASC AFTER move ID: %d, min_dist: %g, pos -> (%g, %g, %g), vel -> (%g, %g, %g)\n", id, min_distance, agent_x, agent_y, agent_z, agent_vx, agent_vy, agent_vz);

  
  FLAMEGPU->setVariable<float>("x",agent_x);
  FLAMEGPU->setVariable<float>("y",agent_y);
  FLAMEGPU->setVariable<float>("z",agent_z);
  FLAMEGPU->setVariable<float>("vx",agent_vx);
  FLAMEGPU->setVariable<float>("vy",agent_vy);
  FLAMEGPU->setVariable<float>("vz",agent_vz);

  return flamegpu::ALIVE;
}