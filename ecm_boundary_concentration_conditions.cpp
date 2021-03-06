FLAMEGPU_AGENT_FUNCTION(ecm_boundary_concentration_conditions, flamegpu::MessageNone, flamegpu::MessageNone) {
  // Agent properties in local register
  int id = FLAMEGPU->getVariable<int>("id");

  // Agent position
  float agent_x = FLAMEGPU->getVariable<float>("x");
  float agent_y = FLAMEGPU->getVariable<float>("y");
  float agent_z = FLAMEGPU->getVariable<float>("z");
    
  // Agent concentration
  const uint8_t N_SPECIES = 2; // WARNING: this variable must be hard coded to have the same value as the one defined in the main python function. TODO: declare it somehow at compile time
  float agent_conc_multi[N_SPECIES] = {}; 
  for (int i = 0; i < N_SPECIES; i++) {
    agent_conc_multi[i] = FLAMEGPU->getVariable<float, N_SPECIES>("concentration_multi", i);
  }
 
  float separation_x_pos = 0.0;
  float separation_x_neg = 0.0;
  float separation_y_pos = 0.0;
  float separation_y_neg = 0.0;
  float separation_z_pos = 0.0;
  float separation_z_neg = 0.0;
  const float ECM_BOUNDARY_INTERACTION_RADIUS = FLAMEGPU->environment.getProperty<float>("ECM_BOUNDARY_INTERACTION_RADIUS");
  const float ECM_BOUNDARY_EQUILIBRIUM_DISTANCE = FLAMEGPU->environment.getProperty<float>("ECM_BOUNDARY_EQUILIBRIUM_DISTANCE");
  float EPSILON = FLAMEGPU->environment.getProperty<float>("EPSILON");
  int DEBUG_PRINTING = FLAMEGPU->environment.getProperty<int>("DEBUG_PRINTING");

  // Get position of the boundaries
  const float COORD_BOUNDARY_X_POS = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",0);
  const float COORD_BOUNDARY_X_NEG = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",1);
  const float COORD_BOUNDARY_Y_POS = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",2);
  const float COORD_BOUNDARY_Y_NEG = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",3);
  const float COORD_BOUNDARY_Z_POS = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",4);
  const float COORD_BOUNDARY_Z_NEG = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",5);
  
  // Get concentration conditions from macroscopic variables
  auto BOUNDARY_CONC_INIT_MULTI = FLAMEGPU->environment.getMacroProperty<float, N_SPECIES, 6>("BOUNDARY_CONC_INIT_MULTI");
  auto BOUNDARY_CONC_FIXED_MULTI = FLAMEGPU->environment.getMacroProperty<float, N_SPECIES, 6>("BOUNDARY_CONC_FIXED_MULTI");

  // Check for ecm-boundary separations
  //this takes into account the distance with respect to the actual boundary position, while forces are calculated with respect to boundary initial position
  separation_x_pos = (agent_x - COORD_BOUNDARY_X_POS); 
  separation_x_neg = (agent_x - COORD_BOUNDARY_X_NEG);
  separation_y_pos = (agent_y - COORD_BOUNDARY_Y_POS);
  separation_y_neg = (agent_y - COORD_BOUNDARY_Y_NEG);
  separation_z_pos = (agent_z - COORD_BOUNDARY_Z_POS);
  separation_z_neg = (agent_z - COORD_BOUNDARY_Z_NEG);
  float separations[6];
  separations[0] = separation_x_pos;
  separations[1] = separation_x_neg;
  separations[2] = separation_y_pos;
  separations[3] = separation_y_neg;
  separations[4] = separation_z_pos;
  separations[5] = separation_z_neg;
  
  //
  float max_conc = 0.0;
  int touching_boundary = 0;
  for (int i = 0; i < N_SPECIES; i++) { // loop through the species
    max_conc = 0.0;             // if an agent is touching several boundaries, the maximum concentration is considered
    for (int j = 0; j < 6; j++) {     // loop through the 6 boundaries
      if ((id == 9) && (DEBUG_PRINTING == 1)){             // print first agent for debugging
        printf("species id: %d, boundary: [%d] , initial conc -> %g  \n", i+1, j+1, (float)BOUNDARY_CONC_INIT_MULTI[i][j]);
        printf("species id: %d, boundary: [%d] , fixed conc -> %g  \n", i+1, j+1, (float)BOUNDARY_CONC_FIXED_MULTI[i][j]);
      }    
      if (fabsf(separations[j]) < (ECM_BOUNDARY_INTERACTION_RADIUS)){
        touching_boundary = 1;
        if (BOUNDARY_CONC_FIXED_MULTI[i][j] > max_conc){
          max_conc = BOUNDARY_CONC_FIXED_MULTI[i][j];
          agent_conc_multi[i] = max_conc; 
        }
        if (BOUNDARY_CONC_INIT_MULTI[i][j] > max_conc){
          max_conc = BOUNDARY_CONC_INIT_MULTI[i][j];
          agent_conc_multi[i] = max_conc; 
        }   
      }    
    }     
  }
  
  if (touching_boundary == 1){
    for (int i = 0; i < N_SPECIES; i++) {
      //printf("agent id: %d, species id: %d, max_conc -> %2.6f, conc -> %2.6f  \n", id, i+1, max_conc, agent_conc_multi[i]);
      FLAMEGPU->setVariable<float, N_SPECIES>("concentration_multi", i, agent_conc_multi[i]);    
    } 
  }

  
  return flamegpu::ALIVE;
}