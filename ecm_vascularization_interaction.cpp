FLAMEGPU_DEVICE_FUNCTION float vec3Length(const float x, const float y, const float z) {
  return sqrtf(x * x + y * y + z * z);
}

FLAMEGPU_AGENT_FUNCTION(ecm_vascularization_interaction, flamegpu::MessageSpatial3D, flamegpu::MessageNone) {
  // Agent properties in local register
  int id = FLAMEGPU->getVariable<int>("id");
  

  // Agent position
  float agent_x = FLAMEGPU->getVariable<float>("x");
  float agent_y = FLAMEGPU->getVariable<float>("y");
  float agent_z = FLAMEGPU->getVariable<float>("z");
  
 
  // Agent concentration
  const uint8_t N_SPECIES = 2; // WARNING: this variable must be hard coded to have the same value as the one defined in the main python function. TODO: declare it somehow at compile time
  float agent_conc_multi[N_SPECIES];
  for (int i = 0; i < N_SPECIES; i++) {
	 agent_conc_multi[i] = FLAMEGPU->getVariable<float, N_SPECIES>("concentration_multi", i);
  }
  
  const float ECM_ECM_EQUILIBRIUM_DISTANCE = FLAMEGPU->environment.getProperty<float>("ECM_ECM_EQUILIBRIUM_DISTANCE");
  const float DELTA_TIME = FLAMEGPU->environment.getProperty<float>("DELTA_TIME");
  float EPSILON = FLAMEGPU->environment.getProperty<float>("EPSILON");
  
  int message_id = 0;
  float message_x = 0.0;
  float message_y = 0.0;
  float message_z = 0.0;
  float message_conc_multi[N_SPECIES]; //initialize values to 0.0
  
  // direction: the vector joining interacting agents
  float dir_x = 0.0; 
  float dir_y = 0.0; 
  float dir_z = 0.0; 
  float distance = 0.0;
  
  float min_distance = ECM_ECM_EQUILIBRIUM_DISTANCE; //initialize to maximum possible value between ECM agents
    
  
  for (const auto &message : FLAMEGPU->message_in(agent_x, agent_y, agent_x)) { // find the closest vascularization agent (if any)
    message_id = message.getVariable<int>("id");
	message_x = message.getVariable<float>("x");
    message_y = message.getVariable<float>("y");
    message_z = message.getVariable<float>("z");
    
    dir_x = agent_x - message_x; 
    dir_y = agent_y - message_y; 
    dir_z = agent_z - message_z; 
    distance = vec3Length(dir_x, dir_y, dir_z);      
               
    if (distance < ECM_ECM_EQUILIBRIUM_DISTANCE) {		
		if (distance < min_distance){
			min_distance = distance;
			for (int i = 0; i < N_SPECIES; i++) {
				message_conc_multi[i] = message.getVariable<float, N_SPECIES>("concentration_multi", i);
			}	
		}   
    }
  }
  
  if (min_distance < ECM_ECM_EQUILIBRIUM_DISTANCE){ // at least one vascularization agent is close
	    for (int i = 0; i < N_SPECIES; i++) {
			if (agent_conc_multi[i] < message_conc_multi[i]){
				// printf("agent %d -> min_distance = %2.6f, own conc before = %2.6f \n", id, min_distance, agent_conc_multi[i]);
				agent_conc_multi[i] += (1.0 - (min_distance / ECM_ECM_EQUILIBRIUM_DISTANCE)) * (message_conc_multi[i] - agent_conc_multi[i]) * DELTA_TIME;
				// printf("agent %d -> vasc_conc = %2.6f, own conc after = %2.6f \n", id, message_conc_multi[i] , agent_conc_multi[i]);
				FLAMEGPU->setVariable<float, N_SPECIES>("concentration_multi", i, agent_conc_multi[i]);
			}			
		}
  }



  return flamegpu::ALIVE;
}