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
  
  return angle;
}


FLAMEGPU_AGENT_FUNCTION(ecm_ecm_interaction, flamegpu::MessageArray3D, flamegpu::MessageNone) {
  // Agent properties in local register
  int id = FLAMEGPU->getVariable<int>("id");

  // Agent position
  float agent_x = FLAMEGPU->getVariable<float>("x");
  float agent_y = FLAMEGPU->getVariable<float>("y");
  float agent_z = FLAMEGPU->getVariable<float>("z");

  // Agent grid position
  uint8_t agent_grid_i = FLAMEGPU->getVariable<uint8_t>("grid_i");
  uint8_t agent_grid_j = FLAMEGPU->getVariable<uint8_t>("grid_j");
  uint8_t agent_grid_k = FLAMEGPU->getVariable<uint8_t>("grid_k");
  
  // Agent velocity
  float agent_vx = FLAMEGPU->getVariable<float>("vx");
  float agent_vy = FLAMEGPU->getVariable<float>("vy");
  float agent_vz = FLAMEGPU->getVariable<float>("vz");
  
  // Agent concentration
  const uint8_t N_SPECIES = 2; // WARNING: this variable must be hard coded to have the same value as the one defined in the main python function. TODO: declare it somehow at compile time
  float agent_conc_multi[N_SPECIES] = {}; 
  for (int i = 0; i < N_SPECIES; i++) {
	 agent_conc_multi[i] = FLAMEGPU->getVariable<float, N_SPECIES>("concentration_multi", i);
  }
  
  // Elastic constant of the ecm 
  //const float k_elast = FLAMEGPU->getVariable<float>("k_elast");
  float k_elast = 0.0; //Equivalent elastic constant of two springs in series (agent and message)
  // Elastic constant and orientation of the fibers
  float agent_k_elast = FLAMEGPU->getVariable<float>("k_elast");
  float agent_orx = FLAMEGPU->getVariable<float>("orx");
  float agent_ory = FLAMEGPU->getVariable<float>("ory");
  float agent_orz = FLAMEGPU->getVariable<float>("orz");
  
  // Dumping constant of the ecm 
  const float d_dumping = FLAMEGPU->getVariable<float>("d_dumping");
  const float ECM_ECM_EQUILIBRIUM_DISTANCE = FLAMEGPU->environment.getProperty<float>("ECM_ECM_EQUILIBRIUM_DISTANCE");
  // Equilibrium distance must be adapted depending on the message grid position respect to the agent. Messages in the Neuman neighbourhood will use the original
  // whereas messages in the diagonals will use increased values.
  float grid_equilibrium_distance = 0.0; 
  
  float agent_fx = 0.0;
  float agent_fy = 0.0;
  float agent_fz = 0.0; 
  float agent_f_extension = 0.0;
  float agent_f_compression = 0.0;
  float agent_elastic_energy = 0.0;

  
  float message_x = 0.0;
  float message_y = 0.0;
  float message_z = 0.0;
  float message_conc_multi[N_SPECIES] = {}; //initialize values to 0.0
  int message_id = 0;
  float message_vx = 0.0;
  float message_vy = 0.0;
  float message_vz = 0.0;
  uint8_t message_grid_i = 0;
  uint8_t message_grid_j = 0;
  uint8_t message_grid_k = 0;
  // Elastic constant and orientation of the message agent
  float message_k_elast = 0.0;
  float message_orx = 0.0;
  float message_ory = 0.0;
  float message_orz = 0.0;

  // Initialize other variables
  float EPSILON = FLAMEGPU->environment.getProperty<float>("EPSILON");
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

  int conn = 0;
  int i_diff = 0;
  int j_diff = 0;
  int k_diff = 0;
  int ct = 0;
  int DEBUG_PRINTING = FLAMEGPU->environment.getProperty<int>("DEBUG_PRINTING");
  
  // Distance data of Neuman neighbourhood. Needed to solve diffusion equation
  float n_up_dist = 0.0; 
  float n_down_dist = 0.0;
  float n_right_dist = 0.0; 
  float n_left_dist = 0.0; 
  float n_front_dist = 0.0; 
  float n_back_dist = 0.0; 

  // Concentration and distance data of Neuman neighbourhood. Needed to solve diffusion equation for multiple species
  float n_up_conc_multi[N_SPECIES] = {};  // concentration of agent on top of current one
  float n_down_conc_multi[N_SPECIES] = {}; 
  float n_right_conc_multi[N_SPECIES] = {};  
  float n_left_conc_multi[N_SPECIES] = {};  
  float n_front_conc_multi[N_SPECIES] = {};  
  float n_back_conc_multi[N_SPECIES] = {};  
  
  const float DELTA_TIME = FLAMEGPU->environment.getProperty<float>("DELTA_TIME");
  

  //printf("Interaction agent %d [%d %d %d]\n", id, agent_grid_i, agent_grid_j, agent_grid_k);

  // Iterate location messages, accumulating relevant data and counts.
  for (const auto &message : FLAMEGPU->message_in(agent_grid_i, agent_grid_j, agent_grid_k)) {
    message_id = message.getVariable<int>("id");
    message_x = message.getVariable<float>("x");
    message_y = message.getVariable<float>("y");
    message_z = message.getVariable<float>("z");
    message_grid_i = message.getVariable<uint8_t>("grid_i");
    message_grid_j = message.getVariable<uint8_t>("grid_j");
    message_grid_k = message.getVariable<uint8_t>("grid_k");
	for (int i = 0; i < N_SPECIES; i++) {
	  message_conc_multi[i] = message.getVariable<float, N_SPECIES>("concentration_multi", i);
	}

    i_diff = abs(agent_grid_i - message_grid_i);
    j_diff = abs(agent_grid_j - message_grid_j);
    k_diff = abs(agent_grid_k - message_grid_k);
    conn = i_diff + j_diff + k_diff;

    /*
    if (id == 9 || id == 10 || id == 13 || id == 14 || id == 25 || id == 26 || id == 29 || id == 30) {
        printf("agent id %d, agent grid [%d %d %d] -> (message %d): message grid [%d %d %d], conn = %d \n", id, agent_grid_i, agent_grid_j, agent_grid_k, message_id, message_grid_i, message_grid_j, message_grid_k, conn);
    }
    */
    // If conn < 2 only the Neuman neighbourhood is checked. conn < 4 checks the 26 surrounding agents
    // ��BEWARE!!: grid domain wraps itself, meaning that agents at the grid boundaries, read messages from opposite boundaries. A grid distance condition must be added to avoid that. 
    
    if ((id != message_id) && (conn < 4) && (i_diff < 2) && (j_diff < 2) && (k_diff < 2)){

        if (conn < 2) {
            grid_equilibrium_distance = ECM_ECM_EQUILIBRIUM_DISTANCE; //Neuman neighbourhood
        }
        else if (conn < 3) {
            grid_equilibrium_distance = sqrtf(2.0) * ECM_ECM_EQUILIBRIUM_DISTANCE; //in-plane diagonals
        }
        else {
            grid_equilibrium_distance = sqrtf(3.0) * ECM_ECM_EQUILIBRIUM_DISTANCE; //corners of the 3x3x3 adjacent agents
        }

        ct++;
        dir_x = agent_x - message_x; 
        dir_y = agent_y - message_y; 
        dir_z = agent_z - message_z; 
        distance = vec3Length(dir_x, dir_y, dir_z);    

		if (conn < 2) {			
			if (message_grid_i < agent_grid_i)				
				n_left_dist = distance;
			if (message_grid_i > agent_grid_i) 
				n_right_dist = distance;
			if (message_grid_j < agent_grid_j) 
				n_back_dist = distance;
			if (message_grid_j > agent_grid_j) 
				n_front_dist = distance;
			if (message_grid_k < agent_grid_k) 
				n_down_dist = distance;
			if (message_grid_k > agent_grid_k) 
				n_up_dist = distance;
				
			// For multiple species diffusion
			for (int i = 0; i < N_SPECIES; i++) {
			    if (message_grid_i < agent_grid_i) 
					n_left_conc_multi[i] = message_conc_multi[i];
				if (message_grid_i > agent_grid_i) 
					n_right_conc_multi[i] = message_conc_multi[i];
				if (message_grid_j < agent_grid_j) 
					n_back_conc_multi[i] = message_conc_multi[i];
				if (message_grid_j > agent_grid_j) 
					n_front_conc_multi[i] = message_conc_multi[i];
				if (message_grid_k < agent_grid_k) 
					n_down_conc_multi[i] = message_conc_multi[i];
				if (message_grid_k > agent_grid_k) 
					n_up_conc_multi[i] = message_conc_multi[i];
			}
		}		
		
		message_vx = message.getVariable<float>("vx");
		message_vy = message.getVariable<float>("vy");
		message_vz = message.getVariable<float>("vz");
		message_k_elast = message.getVariable<float>("k_elast");
		message_orx = message.getVariable<float>("orx");
		message_ory = message.getVariable<float>("ory");
		message_orz = message.getVariable<float>("orz");
		
		// angles between agent orientation and the direction joining agents.
		angle_agent_ori_dir = getAngleBetweenVec(agent_orx,agent_ory,agent_orz,dir_x,dir_y,dir_z);
        angle_message_ori_dir = getAngleBetweenVec(message_orx,message_ory,message_orz,dir_x,dir_y,dir_z);
		cos_ori_agent = fabsf(cosf(angle_agent_ori_dir));
		cos_ori_message = fabsf(cosf(angle_message_ori_dir));
		
		if (cos_ori_agent < EPSILON){
			cos_ori_agent = EPSILON;
		}
		if (cos_ori_message < EPSILON){
			cos_ori_message = EPSILON;
		}
		
		if (id == 9) {
            printf("id1: %d id2: %d cos ori [agent: %g - message: %g] positions [%g %g %g]-[%g %g %g] \n", id, message_id, cos_ori_agent,cos_ori_message,agent_x,agent_y,agent_z,message_x,message_y,message_z);            
        }
		cos_ori_agent = 1.0;
		cos_ori_message = 1.0;
		
		k_elast = (cos_ori_agent * agent_k_elast * cos_ori_message * message_k_elast) / ((cos_ori_agent * agent_k_elast) + (cos_ori_message * message_k_elast));

        cos_x = (1.0 * dir_x + 0.0 * dir_y + 0.0 * dir_z) / distance;
        cos_y = (0.0 * dir_x + 1.0 * dir_y + 0.0 * dir_z) / distance;
        cos_z = (0.0 * dir_x + 0.0 * dir_y + 1.0 * dir_z) / distance;
		
		// angles between agent & message velocity vector and the direction joining them		
		angle_agent_v_dir = getAngleBetweenVec(agent_vx,agent_vy,agent_vz,dir_x,dir_y,dir_z);
		angle_message_v_dir = getAngleBetweenVec(message_vx,message_vy,message_vz,dir_x,dir_y,dir_z);

        // relative speed <0 means particles are getting closer
        relative_speed = vec3Length(agent_vx, agent_vy, agent_vz) * cosf(angle_agent_v_dir) - vec3Length(message_vx, message_vy, message_vz) * cosf(angle_message_v_dir);
        // if total_f > 0, agents are attracted, if <0 agents are repelled
        total_f = (distance - grid_equilibrium_distance) * (k_elast) + d_dumping * relative_speed;

        if (total_f < 0) {
            agent_f_compression += total_f;
        }
        else {
            agent_f_extension += total_f;
        }

        agent_elastic_energy += 0.5 * (total_f * total_f) / k_elast;

        agent_fx += -1 * total_f * cos_x;
        agent_fy += -1 * total_f * cos_y;
        agent_fz += -1 * total_f * cos_z;

        if (DEBUG_PRINTING == 1 && (id == 9 || id == 10 || id == 13 || id == 14 || id == 25 || id == 26 || id == 29 || id == 30)) {
            printf("ECM interaction [id1: %d - id2: %d] agent_pos (%2.6f, %2.6f, %2.6f), message_pos (%2.6f, %2.6f, %2.6f)\n", id, message_id, agent_x, agent_y, agent_z, message_x, message_y, message_z);
            printf("ECM interaction id1: %d - id2: %d distance -> (%2.6f)\n", id, message_id, distance);
            printf("ECM interaction id1: %d - id2: %d total_f -> %2.6f (%2.6f , %2.6f, %2.6f)\n", id, message_id, total_f, -1 * total_f * cos_x, -1 * total_f * cos_y, -1 * total_f * cos_z);
        }
    }
  }
  if (threadIdx.x%2 == 0) {
    //printf("Array3D for agent %d read %d messages!\n", id, ct);
  }
  //printf("Array3D for agent %d read %d messages! grid [%d %d %d], pos (%2.6f , %2.6f, %2.6f) \n", id, ct, agent_grid_i, agent_grid_j, agent_grid_k, agent_x, agent_y, agent_z);

  //Apply diffusion equation
  
  float R = 0.0; // reactive term. Unused for now
  float dx = ((n_left_dist > 0.0) & (n_right_dist > 0.0)) ? (n_left_dist + n_right_dist) / 2.0 : fmaxf(n_left_dist,n_right_dist);
  float dy = ((n_front_dist > 0.0) & (n_back_dist > 0.0)) ? (n_front_dist + n_back_dist) / 2.0 : fmaxf(n_front_dist,n_back_dist);
  float dz = ((n_up_dist > 0.0) & (n_down_dist > 0.0)) ? (n_up_dist + n_down_dist) / 2.0 : fmaxf(n_up_dist,n_down_dist);
  
  //Apply diffusion equation for multiple species
  float agent_conc_prev_multi[N_SPECIES] = {}; 
  for (int i = 0; i < N_SPECIES; i++) {
	float DIFFUSION_COEFF = FLAMEGPU->environment.getProperty<float>("DIFFUSION_COEFF_MULTI",i);
	float Fx = DIFFUSION_COEFF * DELTA_TIME / powf(dx, 2.0);
    float Fy = DIFFUSION_COEFF * DELTA_TIME / powf(dy, 2.0);
    float Fz = DIFFUSION_COEFF * DELTA_TIME / powf(dz, 2.0);
	agent_conc_prev_multi[i] = agent_conc_multi[i];
    agent_conc_multi[i] = agent_conc_prev_multi[i] + Fx * (n_left_conc_multi[i] - (2 * agent_conc_prev_multi[i]) + n_right_conc_multi[i]) + Fy * (n_front_conc_multi[i] - (2 * agent_conc_prev_multi[i]) + n_back_conc_multi[i]) + Fz * (n_up_conc_multi[i] - (2 * agent_conc_prev_multi[i]) + n_down_conc_multi[i]) + R * DELTA_TIME;
    FLAMEGPU->setVariable<float, N_SPECIES>("concentration_multi", i, agent_conc_multi[i]);
	if ((id > 8) && (DEBUG_PRINTING == 1)){  
		printf("DIFFUSION for agent %d, species %d, [dx,dy,dz] = [%2.6f , %2.6f, %2.6f], [Fx,Fy,Fz] = [%2.6f , %2.6f, %2.6f] \n", id, i+1, dx, dy, dz, Fx, Fy, Fz);
		printf("agent %d: MULTI left conc = %2.6f, right conc = %2.6f \n", id, n_left_conc_multi[i], n_right_conc_multi[i]);
		printf("agent %d: MULTI front conc = %2.6f, back conc = %2.6f \n", id, n_front_conc_multi[i], n_back_conc_multi[i]);
		printf("agent %d: MULTI up conc = %2.6f, down conc = %2.6f \n", id, n_up_conc_multi[i], n_down_conc_multi[i]);
		printf("agent %d: MULTI conc prev = %2.6f, current conc = %2.6f \n", id, agent_conc_prev_multi[i], agent_conc_multi[i]); 
	}	
  }
  
  // Re-orientation of fibers towards the direction of the resultant force
  // inc_dir = ECM_ORIENTATION_RATE * DELTA_TIME * cross(agent_ori,cross(force_dir,agent_ori))
  
  float force_magnitude = vec3Length(agent_fx,agent_fy,agent_fz);
  
  if (force_magnitude > 0){
	  const float ECM_ORIENTATION_RATE = FLAMEGPU->environment.getProperty<float>("ECM_ORIENTATION_RATE");
	  float inc_dir_x = 0.0;
	  float inc_dir_y = 0.0;
	  float inc_dir_z = 0.0;
	  float dir_fx = agent_fx / force_magnitude;
	  float dir_fy = agent_fy / force_magnitude;
	  float dir_fz = agent_fz / force_magnitude;
	  float cos_force_ori = cosf(getAngleBetweenVec(agent_orx,agent_ory,agent_orz,dir_fx,dir_fy,dir_fz));
	  if (cos_force_ori < 0.0){ // invert direction to find the closest angle between force and orientation directions
		  dir_fx = -1 * dir_fx;
		  dir_fy = -1 * dir_fy;
		  dir_fz = -1 * dir_fz;
	  }
	  float tmpx = 0.0;
	  float tmpy = 0.0;
	  float tmpz = 0.0;
	  
	  if (id == 9 && (DEBUG_PRINTING == 1)) {
		printf("ORI ANTES id1: %d  ori  [%g %g %g] \n", id, agent_orx,agent_ory,agent_orz); 
		printf("inc_dir ANTES id1: %d  inc_dir  [%g %g %g] \n", id, inc_dir_x,inc_dir_y,inc_dir_z); 
		printf("dir_f ANTES id1: %d  dirf  [%g %g %g] \n", id, dir_fx,dir_fy,dir_fz); 		
	  }
	  
	  vec3CrossProd(tmpx, tmpy, tmpz, dir_fx, dir_fy, dir_fz, agent_orx, agent_ory, agent_orz);
	  vec3CrossProd(inc_dir_x, inc_dir_y, inc_dir_z, agent_orx, agent_ory, agent_orz, tmpx, tmpy, tmpz); 
	  
	  agent_orx += inc_dir_x * ECM_ORIENTATION_RATE * DELTA_TIME;
	  agent_ory += inc_dir_y * ECM_ORIENTATION_RATE * DELTA_TIME;
	  agent_orz += inc_dir_z * ECM_ORIENTATION_RATE * DELTA_TIME;
	  
	  FLAMEGPU->setVariable<float>("alignment", fabsf(cos_force_ori));
  }
  
  
    
  FLAMEGPU->setVariable<float>("orx", agent_orx);
  FLAMEGPU->setVariable<float>("ory", agent_ory);
  FLAMEGPU->setVariable<float>("orz", agent_orz);
  FLAMEGPU->setVariable<float>("fx", agent_fx);
  FLAMEGPU->setVariable<float>("fy", agent_fy);
  FLAMEGPU->setVariable<float>("fz", agent_fz);
  FLAMEGPU->setVariable<float>("f_extension", agent_f_extension);
  FLAMEGPU->setVariable<float>("f_compression", agent_f_compression);
  FLAMEGPU->setVariable<float>("elastic_energy", agent_elastic_energy);


  return flamegpu::ALIVE;
}