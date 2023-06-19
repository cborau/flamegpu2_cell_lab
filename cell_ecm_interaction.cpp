FLAMEGPU_DEVICE_FUNCTION void vec3CrossProd(float &x, float &y, float &z, float x1, float y1, float z1, float x2, float y2, float z2) {
  x = (y1 * z2 - z1 * y2);
  y = (z1 * x2 - x1 * z2);
  z = (x1 * y2 - y1 * x2);
}
FLAMEGPU_DEVICE_FUNCTION void vec3Div(float &x, float &y, float &z, const float divisor) {
  x /= divisor;
  y /= divisor;
  z /= divisor;
}
FLAMEGPU_DEVICE_FUNCTION float vec3Length(const float x, const float y, const float z) {
  return sqrtf(x * x + y * y + z * z);
}
FLAMEGPU_DEVICE_FUNCTION void vec3Normalize(float &x, float &y, float &z) {
  float length = vec3Length(x, y, z);
  vec3Div(x, y, z, length);
}
FLAMEGPU_DEVICE_FUNCTION float getAngleBetweenVec(const float x1, const float y1, const float z1, const float x2, const float y2, const float z2) {
  float dot_dir = x1 * x2 + y1 * y2 + z1 * z2;
  float cross_x_dir = 0.0;
  float cross_y_dir = 0.0;
  float cross_z_dir = 0.0;
  float angle = 0.0;
  float EPSILON = 0.0000000001;
  vec3CrossProd(cross_x_dir, cross_y_dir, cross_z_dir, x1, y1, z1, x2, y2, z2);
  float det_dir = vec3Length(cross_x_dir, cross_y_dir, cross_z_dir);
  if (fabsf(dot_dir) > EPSILON) {
    angle = atan2f(det_dir, dot_dir);
  }
  else {
    angle = 0.0;
  }
  
  return angle; //in radians
}
// This function includes cell reorientation after deformations (Cell agent is the caller, ECM agents are the messages). WARNING: not to be confused with ecm_cell_interaction, which computes the ECM deformation
FLAMEGPU_AGENT_FUNCTION(cell_ecm_interaction, flamegpu::MessageArray3D, flamegpu::MessageNone) {
  // Agent properties in local register
  int id = FLAMEGPU->getVariable<int>("id");
  
  // Agent position
  float agent_x = FLAMEGPU->getVariable<float>("x");
  float agent_y = FLAMEGPU->getVariable<float>("y");
  float agent_z = FLAMEGPU->getVariable<float>("z");
  
   // Agent velocity
  float agent_vx = FLAMEGPU->getVariable<float>("vx");
  float agent_vy = FLAMEGPU->getVariable<float>("vy");
  float agent_vz = FLAMEGPU->getVariable<float>("vz");
   
  // Agen orientation
  float agent_orx = FLAMEGPU->getVariable<float>("orx");
  float agent_ory = FLAMEGPU->getVariable<float>("ory");
  float agent_orz = FLAMEGPU->getVariable<float>("orz");
  
  const float MAX_SEARCH_RADIUS_CELLS = FLAMEGPU->environment.getProperty<float>("MAX_SEARCH_RADIUS_CELLS");
  const float DELTA_TIME = FLAMEGPU->environment.getProperty<float>("DELTA_TIME");
  float EPSILON = FLAMEGPU->environment.getProperty<float>("EPSILON");
  int INCLUDE_CELL_ORIENTATION = FLAMEGPU->environment.getProperty<int>("INCLUDE_CELL_ORIENTATION");
  
  
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
  
  int message_id = 0;
  float message_x = 0.0;
  float message_y = 0.0;
  float message_z = 0.0;
  float message_orx = 0.0;
  float message_ory = 0.0;
  float message_orz = 0.0;
  
  // direction: the vector joining interacting agents
  float dir_x = 0.0; 
  float dir_y = 0.0; 
  float dir_z = 0.0; 
  float distance = 0.0;
   
  float dir_max_strain_x = 0.0;
  float dir_max_strain_y = 0.0;
  float dir_max_strain_z = 0.0;
  float max_strain = 0.0;
  
  // The closest ECM agent
  const auto message = FLAMEGPU->message_in.at(agent_grid_i, agent_grid_j, agent_grid_k);

  message_id = message.getVariable<int>("id");
  message_x = message.getVariable<float>("x");
  message_y = message.getVariable<float>("y");
  message_z = message.getVariable<float>("z");
  message_orx = message.getVariable<float>("orx");
  message_ory = message.getVariable<float>("ory");
  message_orz = message.getVariable<float>("orz");

	
	//printf("agent %d -> message xyz (%d) = %2.6f, %2.6f, %2.6f \n", id, message_id, message_x, message_y, message_z);
	    
  dir_x = agent_x - message_x; 
  dir_y = agent_y - message_y; 
  dir_z = agent_z - message_z; 
  distance = vec3Length(dir_x, dir_y, dir_z); 
  float strain = (MAX_SEARCH_RADIUS_CELLS - distance) / MAX_SEARCH_RADIUS_CELLS;
		
  if (strain > max_strain){
		max_strain = strain;
		dir_max_strain_x = dir_x / distance;
		dir_max_strain_y = dir_y / distance;
		dir_max_strain_z = dir_z / distance;
  }
 
    
  if ((max_strain > EPSILON) && (INCLUDE_CELL_ORIENTATION == 1)){
	const float CELL_ORIENTATION_RATE = FLAMEGPU->environment.getProperty<float>("CELL_ORIENTATION_RATE");
	float inc_dir_x = 0.0;
	float inc_dir_y = 0.0;
	float inc_dir_z = 0.0;
	float dir_fx = 0.0;
	float dir_fy = 0.0;
	float dir_fz = 0.0;
	
	// reorient towards ECM orientation direction 
	dir_fx = message_orx;
	dir_fy = message_ory;
	dir_fz = message_orz;
	float cos_force_ori = cosf(getAngleBetweenVec(agent_orx,agent_ory,agent_orz,dir_fx,dir_fy,dir_fz));
	if (cos_force_ori < 0.0){ // invert direction to find the closest angle between force and orientation directions
	  dir_fx = -1 * dir_fx;
	  dir_fy = -1 * dir_fy;
	  dir_fz = -1 * dir_fz;
	}
	  
	// Multiply again by strain magnitude to make orientation rate strain-dependent
	dir_fx *= max_strain;
	dir_fy *= max_strain;
	dir_fz *= max_strain;
	  
	float tmpx = 0.0;
	float tmpy = 0.0;
	float tmpz = 0.0;
	  	  
	vec3CrossProd(tmpx, tmpy, tmpz, dir_fx, dir_fy, dir_fz, agent_orx, agent_ory, agent_orz);
	vec3CrossProd(inc_dir_x, inc_dir_y, inc_dir_z, agent_orx, agent_ory, agent_orz, tmpx, tmpy, tmpz); 
	
	printf("cell %d -> ecm %d, cell_orientation_rate = %2.6f, delta_time = %2.6f, inc_dir_x = %2.6f, inc_dir_y = %2.6f, inc_dir_z = %2.6f, \n", id, message_id, CELL_ORIENTATION_RATE, DELTA_TIME, inc_dir_x, inc_dir_y, inc_dir_z);

	  
	agent_orx += inc_dir_x * CELL_ORIENTATION_RATE * DELTA_TIME;
	agent_ory += inc_dir_y * CELL_ORIENTATION_RATE * DELTA_TIME;
	agent_orz += inc_dir_z * CELL_ORIENTATION_RATE * DELTA_TIME;
	  
	float ori_length = vec3Length(agent_orx,agent_ory,agent_orz);	  
	vec3Div(agent_orx, agent_ory, agent_orz, ori_length);
	  	  
	FLAMEGPU->setVariable<float>("orx", agent_orx);
	FLAMEGPU->setVariable<float>("ory", agent_ory);
    FLAMEGPU->setVariable<float>("orz", agent_orz);
	
	float new_cos_force_ori = cosf(getAngleBetweenVec(agent_orx,agent_ory,agent_orz,dir_fx/max_strain,dir_fy/max_strain,dir_fz/max_strain));
	
	FLAMEGPU->setVariable<float>("alignment", fabsf(new_cos_force_ori));
  }


  return flamegpu::ALIVE;
}