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
// This function computes the ECM deformation (ECM agent is the caller, Cells are the messages) due to the action of cell agents. WARNING: not to be confused with cell_ecm_interaction, which includes cell reorientation after deformations
FLAMEGPU_AGENT_FUNCTION(ecm_cell_interaction, flamegpu::MessageSpatial3D, flamegpu::MessageNone) {
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
  
  // Agent force
  float agent_fx = FLAMEGPU->getVariable<float>("fx");
  float agent_fy = FLAMEGPU->getVariable<float>("fy");
  float agent_fz = FLAMEGPU->getVariable<float>("fz");
  float agent_f_extension = FLAMEGPU->getVariable<float>("f_extension");
  float agent_f_compression = FLAMEGPU->getVariable<float>("f_compression");
  float agent_elastic_energy = FLAMEGPU->getVariable<float>("elastic_energy");
  float agent_fx_abs = 0.0; // if there are opposing forces (F) in the same direction, agent_fx = 0, but agent_fx_abs = 2*F
  float agent_fy_abs = 0.0;
  float agent_fz_abs = 0.0; 
  float agent_k_elast = FLAMEGPU->getVariable<float>("k_elast");
  
  const float MAX_SEARCH_RADIUS_CELLS = FLAMEGPU->environment.getProperty<float>("MAX_SEARCH_RADIUS_CELLS");
  const float CELL_RADIUS = FLAMEGPU->environment.getProperty<float>("CELL_RADIUS");
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
    
    //printf("agent %d -> message xyz (%d) = %2.6f, %2.6f, %2.6f \n", id, message_id, message_x, message_y, message_z);
        
    dir_x = agent_x - message_x; 
    dir_y = agent_y - message_y; 
    dir_z = agent_z - message_z; 
    distance = vec3Length(dir_x, dir_y, dir_z); 

    if (distance < MAX_SEARCH_RADIUS_CELLS) {
        // angles between agent orientation and the direction joining agents.
      angle_message_ori_dir = getAngleBetweenVec(message_orx,message_ory,message_orz,dir_x,dir_y,dir_z);
      cos_ori_agent = 1.0; // ECM orientation is neglected here. Only cell orientation matters. 
      cos_ori_message = fabsf(cosf(angle_message_ori_dir));
      
      if (INCLUDE_CELL_ORIENTATION != 1){
        cos_ori_message = 1.0;
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
	  
	  float k_elast = (agent_k_elast * message_k_elast) / (agent_k_elast + message_k_elast); //series springs

      // relative speed <0 means particles are getting closer
      relative_speed = vec3Length(agent_vx, agent_vy, agent_vz) * cosf(angle_agent_v_dir) - vec3Length(message_vx, message_vy, message_vz) * cosf(angle_message_v_dir);
      // if total_f > 0, agents are attracted, if <0 agents are repelled. Since cells are always contracting, this should be always positive unless the ECM overlaps cell radius
      if (distance < CELL_RADIUS){
        float offset = (MAX_SEARCH_RADIUS_CELLS - CELL_RADIUS) * (k_elast); 
        total_f = ((offset / CELL_RADIUS) * distance -1 * (offset)) + message_d_dumping * relative_speed;
      } 
	  else {
         total_f = +1 * (MAX_SEARCH_RADIUS_CELLS - distance) * (k_elast) + message_d_dumping * relative_speed;
      }
      
      total_f *= cos_ori_message;
      // TODO: rewrite the force equation so that the ECM agent does not cross the cell radius
      
      //printf("ECM agent %d - cell %d -> distance = %2.6f, k_elast = %2.6f, total_f = %2.6f, relative_speed = %2.6f \n", id, message_id, distance, message_k_elast, total_f, relative_speed);
      
      
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

      agent_elastic_energy += 0.5 * (total_f * total_f) / message_k_elast;

      agent_fx += -1 * total_f * cos_x; // minus comes from the direction definition (agent-message)
      agent_fy += -1 * total_f * cos_y;
      agent_fz += -1 * total_f * cos_z;
  
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