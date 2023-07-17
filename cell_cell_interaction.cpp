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
// This function computes the interaction between cells
FLAMEGPU_AGENT_FUNCTION(cell_cell_interaction, flamegpu::MessageSpatial3D, flamegpu::MessageNone) {
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
  
  // Agent orientation
  float agent_orx = FLAMEGPU->getVariable<float>("orx");
  float agent_ory = FLAMEGPU->getVariable<float>("ory");
  float agent_orz = FLAMEGPU->getVariable<float>("orz");
  
  // Agent force
  float agent_fx = 0.0;
  float agent_fy = 0.0;
  float agent_fz = 0.0;
  float agent_f_extension = 0.0;
  float agent_f_compression = 0.0;
  float agent_elastic_energy = 0.0;
  float agent_fx_abs = 0.0; // if there are opposing forces (F) in the same direction, agent_fx = 0, but agent_fx_abs = 2*F
  float agent_fy_abs = 0.0;
  float agent_fz_abs = 0.0; 
  
  // Agent radius
  float agent_radius = FLAMEGPU->getVariable<float>("radius");
  
  const float MAX_SEARCH_RADIUS_CELL_CELL_INTERACTION = FLAMEGPU->environment.getProperty<float>("MAX_SEARCH_RADIUS_CELL_CELL_INTERACTION");
  const float DELTA_TIME = FLAMEGPU->environment.getProperty<float>("DELTA_TIME");
  float EPSILON = FLAMEGPU->environment.getProperty<float>("EPSILON");
  int INCLUDE_CELL_ORIENTATION = FLAMEGPU->environment.getProperty<int>("INCLUDE_CELL_ORIENTATION");
  
  int message_id = 0;
  float message_x = 0.0;
  float message_y = 0.0;
  float message_z = 0.0;
  float message_vx = 0.0;
  float message_vy = 0.0;
  float message_vz = 0.0;
  float message_orx = 0.0;
  float message_ory = 0.0;
  float message_orz = 0.0;
  float message_k_elast = 0.0;
  float message_d_dumping = 0.0;
  float message_radius = 0.0;
  
  // direction: the vector joining interacting agents
  float dir_x = 0.0; 
  float dir_y = 0.0; 
  float dir_z = 0.0; 
  float distance = 0.0;
  
  // director cosines (with respect to global axis) of the direction vector
  float cos_x = 0.0;
  float cos_y = 0.0;
  float cos_z = 0.0;
   // angle (in radians) between agent orientation vector and direction vector
  float angle_agent_ori_dir = 0.0;
  float angle_message_ori_dir = 0.0;
  float cos_ori_agent = 0.0;
  float cos_ori_message = 0.0;
  // angle (in radians) between agent velocity vector and direction vector
  float angle_agent_v_dir = 0.0;
  float angle_message_v_dir = 0.0;
  // relative speed between agents
  float relative_speed = 0.0;
  // total force between agents
  float total_f = 0.0;
  
  float agent_k_elast = FLAMEGPU->getVariable<float>("k_elast");
  float agent_d_dumping = FLAMEGPU->getVariable<float>("d_dumping");
  
  for (const auto &message : FLAMEGPU->message_in(agent_x, agent_y, agent_z)) { // find cell agents within radius
    message_id = message.getVariable<int>("id");
	message_x = message.getVariable<float>("x");
    message_y = message.getVariable<float>("y");
    message_z = message.getVariable<float>("z");
	message_vx = message.getVariable<float>("vx");
	message_vy = message.getVariable<float>("vy");
	message_vz = message.getVariable<float>("vz");
	message_orx = message.getVariable<float>("orx");
	message_ory = message.getVariable<float>("ory");
	message_orz = message.getVariable<float>("orz");
	message_k_elast = message.getVariable<float>("k_elast");
	message_d_dumping = message.getVariable<float>("d_dumping");
	message_radius = message.getVariable<float>("radius");
	
	//printf("agent %d -> message xyz (%d) = %2.6f, %2.6f, %2.6f \n", id, message_id, message_x, message_y, message_z);
	    
    dir_x = agent_x - message_x; 
    dir_y = agent_y - message_y; 
    dir_z = agent_z - message_z; 
    distance = vec3Length(dir_x, dir_y, dir_z); 

	if ((distance < MAX_SEARCH_RADIUS_CELL_CELL_INTERACTION) && (distance > 0.0)) {
		// angles between agent orientation and the direction joining agents.
		angle_agent_ori_dir = getAngleBetweenVec(agent_orx,agent_ory,agent_orz,dir_x,dir_y,dir_z);
		angle_message_ori_dir = getAngleBetweenVec(message_orx,message_ory,message_orz,dir_x,dir_y,dir_z);
		cos_ori_agent = fabsf(cosf(angle_agent_ori_dir));
		cos_ori_message = fabsf(cosf(angle_message_ori_dir));
		
		if (INCLUDE_CELL_ORIENTATION != 1){
			cos_ori_agent = 1.0;
			cos_ori_message = 1.0;
		}
		
		if (cos_ori_agent < EPSILON){
			cos_ori_agent = EPSILON;
		}			

		if (cos_ori_message < EPSILON){
			cos_ori_message = EPSILON;
		}	
		
		cos_x = (1.0 * dir_x + 0.0 * dir_y + 0.0 * dir_z) / distance;
		cos_y = (0.0 * dir_x + 1.0 * dir_y + 0.0 * dir_z) / distance;
		cos_z = (0.0 * dir_x + 0.0 * dir_y + 1.0 * dir_z) / distance;
		
		// angles between agent & message velocity vector and the direction joining them		
		angle_agent_v_dir = getAngleBetweenVec(agent_vx,agent_vy,agent_vz,dir_x,dir_y,dir_z);
		angle_message_v_dir = getAngleBetweenVec(message_vx,message_vy,message_vz,dir_x,dir_y,dir_z);
		
		float k_elast = (cos_ori_agent * agent_k_elast * cos_ori_message * message_k_elast) / ((cos_ori_agent * agent_k_elast) + (cos_ori_message * message_k_elast));


		// relative speed <0 means particles are getting closer
		relative_speed = vec3Length(agent_vx, agent_vy, agent_vz) * cosf(angle_agent_v_dir) - vec3Length(message_vx, message_vy, message_vz) * cosf(angle_message_v_dir);
		// if total_f > 0, agents are attracted, if <0 agents are repelled. Since cells are always contracting, this should be always positive unless the ECM overlaps cell radius
		if (distance < (agent_radius + message_radius)){
			float offset = (MAX_SEARCH_RADIUS_CELL_CELL_INTERACTION - (agent_radius + message_radius)) * (k_elast); 
			total_f = ((offset / (agent_radius + message_radius)) * distance -1 * (offset)) + message_d_dumping * relative_speed;
		} 
		else {
			total_f = +1 * (MAX_SEARCH_RADIUS_CELL_CELL_INTERACTION - distance) * (k_elast) + message_d_dumping * relative_speed;
		}
		
				
		//printf("CELL %d - CELL %d -> distance = %2.6f, k_elast = %2.6f, total_f = %2.6f, relative_speed = %2.6f \n", id, message_id, distance, message_k_elast, total_f, relative_speed);
		
		
		if (total_f < 0) {
			agent_f_compression += total_f;
		}
		else {
			agent_f_extension += total_f;
			// store the absolute extensions in each direction
			agent_fx_abs += fabsf(total_f * cos_x);
			agent_fy_abs += fabsf(total_f * cos_y);
			agent_fz_abs += fabsf(total_f * cos_z);
		}

		agent_elastic_energy += 0.5 * (total_f * total_f) / k_elast;
		//printf("F antes CELL %d -> fx = %2.6f, fy = %2.6f, fz = %2.6f,\n", id, agent_fx, agent_fy, agent_fz);


		agent_fx += -1 * total_f * cos_x; // minus comes from the direction definition (agent-message)
		agent_fy += -1 * total_f * cos_y;
		agent_fz += -1 * total_f * cos_z;
		
		//printf("F despues CELL %d -> fx = %2.6f, fy = %2.6f, fz = %2.6f,\n", id, agent_fx, agent_fy, agent_fz);
  
    }	
	

  }
  
  FLAMEGPU->setVariable<float>("fx", agent_fx);
  FLAMEGPU->setVariable<float>("fy", agent_fy);
  FLAMEGPU->setVariable<float>("fz", agent_fz);
  FLAMEGPU->setVariable<float>("f_extension", agent_f_extension);
  FLAMEGPU->setVariable<float>("f_compression", agent_f_compression);
  FLAMEGPU->setVariable<float>("elastic_energy", agent_elastic_energy);


  return flamegpu::ALIVE;
}