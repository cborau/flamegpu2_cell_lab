FLAMEGPU_DEVICE_FUNCTION void vec3Div(float &x, float &y, float &z, const float divisor) {
  x /= divisor;
  y /= divisor;
  z /= divisor;
}
FLAMEGPU_DEVICE_FUNCTION float vec3Length(const float x, const float y, const float z) {
  return sqrtf(x * x + y * y + z * z);
}
FLAMEGPU_AGENT_FUNCTION(cell_move, flamegpu::MessageNone, flamegpu::MessageNone) {
  
  int id = FLAMEGPU->getVariable<int>("id");
  //Agent position vector
  float agent_x = FLAMEGPU->getVariable<float>("x");
  float agent_y = FLAMEGPU->getVariable<float>("y");
  float agent_z = FLAMEGPU->getVariable<float>("z");

  int DEBUG_PRINTING = FLAMEGPU->environment.getProperty<int>("DEBUG_PRINTING");

  // Agent velocity is reset every step as cells are not tied to each other
  float agent_vx = 0.0;
  float agent_vy = 0.0;
  float agent_vz = 0.0;
  const float CELL_SPEED_REF = FLAMEGPU->environment.getProperty<float>("CELL_SPEED_REF");
  
  // Agent orientation
  float agent_orx = FLAMEGPU->getVariable<float>("orx");
  float agent_ory = FLAMEGPU->getVariable<float>("ory");
  float agent_orz = FLAMEGPU->getVariable<float>("orz");
     
  // Mass of the cell agents (current value shared with ECM. Modify rest of paramenters instead)
  float mass = 1.0;

  //Forces acting on the agent
  float agent_fx = FLAMEGPU->getVariable<float>("fx");
  float agent_fy = FLAMEGPU->getVariable<float>("fy");
  float agent_fz = FLAMEGPU->getVariable<float>("fz");
  
  //printf("ANTES cell %d: x: %2.3f, y: %2.3f, z: %2.3f vx: %2.3f, vy: %2.3f, vz: %2.3f fx: %2.3f, fy: %2.3f, fz: %2.3f\n",id, agent_x,agent_y,agent_z, agent_vx,agent_vy,agent_vz, agent_fx,agent_fy,agent_fz);

  //Get the new position and velocity: 
  // a(t) = f(t) / m;
  // v(t) = v(t-1) + a(t) * dt; 
  // x(t) = x(t-1) + v(t) * dt
  const float DELTA_TIME = FLAMEGPU->environment.getProperty<float>("DELTA_TIME");
  // Random movement
  float rand_dir_x = FLAMEGPU->random.uniform<float>(-1.0,1.0);
  float rand_dir_y = FLAMEGPU->random.uniform<float>(-1.0,1.0);
  float rand_dir_z = FLAMEGPU->random.uniform<float>(-1.0,1.0); 
  float rand_dir_length = vec3Length(rand_dir_x,rand_dir_y,rand_dir_z);      
  vec3Div(rand_dir_x, rand_dir_y, rand_dir_z, rand_dir_length);
  float rand_dir_factor = FLAMEGPU->random.uniform<float>(0.0,2.0); // from 0 to 2 times the ref speed in a random direction 
  // Migration in the direction of cell orientation
  float rand_prot_factor = FLAMEGPU->random.uniform<float>(0.0,2.0); // from 0 to 2 times the ref speed in the direction of cell orientation
  //printf("cell %d: rand_prot_factor: %2.3f",id,rand_prot_factor);
  
  agent_vx += (agent_fx / mass) * DELTA_TIME;
  agent_vx += CELL_SPEED_REF * rand_dir_factor * rand_dir_x;
  agent_vx += CELL_SPEED_REF * rand_prot_factor * agent_orx;
  agent_x += agent_vx * DELTA_TIME;
  
  agent_vy += (agent_fy / mass) * DELTA_TIME;
  agent_vy += CELL_SPEED_REF * rand_dir_factor * rand_dir_y;
  agent_vy += CELL_SPEED_REF * rand_prot_factor * agent_ory;
  agent_y += agent_vy * DELTA_TIME;
  
  agent_vz += (agent_fz / mass) * DELTA_TIME;
  agent_vz += CELL_SPEED_REF * rand_dir_factor * rand_dir_z;
  agent_vz += CELL_SPEED_REF * rand_prot_factor * agent_orz;
  agent_z += agent_vz * DELTA_TIME;
  
  //printf("DESPUES cell %d: x: %2.3f, y: %2.3f, z: %2.3f vx: %2.3f, vy: %2.3f, vz: %2.3f fx: %2.3f, fy: %2.3f, fz: %2.3f\n",id, agent_x,agent_y,agent_z, agent_vx,agent_vy,agent_vz, agent_fx,agent_fy,agent_fz);

  // Check boundaries
  int PERIODIC_BOUNDARIES_FOR_CELLS = FLAMEGPU->environment.getProperty<int>("PERIODIC_BOUNDARIES_FOR_CELLS");
  const float COORD_BOUNDARY_X_POS = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",0);
  const float COORD_BOUNDARY_X_NEG = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",1);
  const float COORD_BOUNDARY_Y_POS = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",2);
  const float COORD_BOUNDARY_Y_NEG = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",3);
  const float COORD_BOUNDARY_Z_POS = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",4);
  const float COORD_BOUNDARY_Z_NEG = FLAMEGPU->environment.getProperty<float>("COORDS_BOUNDARIES",5);
  
  if (PERIODIC_BOUNDARIES_FOR_CELLS == 1){
    if (agent_x > COORD_BOUNDARY_X_POS){
      agent_x -= (COORD_BOUNDARY_X_POS - COORD_BOUNDARY_X_NEG);
    }
    if (agent_x < COORD_BOUNDARY_X_NEG){
      agent_x += (COORD_BOUNDARY_X_POS - COORD_BOUNDARY_X_NEG);
    }
    if (agent_y > COORD_BOUNDARY_Y_POS){
      agent_y -= (COORD_BOUNDARY_Y_POS - COORD_BOUNDARY_Y_NEG);
    }
    if (agent_y < COORD_BOUNDARY_Y_NEG){
      agent_y += (COORD_BOUNDARY_Y_POS - COORD_BOUNDARY_Y_NEG);
    }
    if (agent_z > COORD_BOUNDARY_Z_POS){
      agent_z -= (COORD_BOUNDARY_Z_POS - COORD_BOUNDARY_Z_NEG);
    }
    if (agent_z < COORD_BOUNDARY_Z_NEG){
      agent_z += (COORD_BOUNDARY_Z_POS - COORD_BOUNDARY_Z_NEG);
    } 
  }
  else { // don't let them go out of boundaries
    if (agent_x > COORD_BOUNDARY_X_POS){
      agent_x = COORD_BOUNDARY_X_POS;
      agent_vx = 0.0;
    }
    if (agent_x < COORD_BOUNDARY_X_NEG){
      agent_x = COORD_BOUNDARY_X_NEG;
      agent_vx = 0.0;
    }
    if (agent_y > COORD_BOUNDARY_Y_POS){
      agent_y = COORD_BOUNDARY_Y_POS;
      agent_vy = 0.0;
    }
    if (agent_y < COORD_BOUNDARY_Y_NEG){
      agent_y = COORD_BOUNDARY_Y_NEG;
      agent_vy = 0.0;
    }
    if (agent_z > COORD_BOUNDARY_Z_POS){
      agent_z = COORD_BOUNDARY_Z_POS;
      agent_vz = 0.0;
    }
    if (agent_z < COORD_BOUNDARY_Z_NEG){
      agent_z = COORD_BOUNDARY_Z_NEG;
      agent_vz = 0.0;
    }  
  }

  //Update the agents position and velocity
  FLAMEGPU->setVariable<float>("x",agent_x);
  FLAMEGPU->setVariable<float>("y",agent_y);
  FLAMEGPU->setVariable<float>("z",agent_z);
  FLAMEGPU->setVariable<float>("vx",agent_vx);
  FLAMEGPU->setVariable<float>("vy",agent_vy);
  FLAMEGPU->setVariable<float>("vz",agent_vz);
  

  return flamegpu::ALIVE;
}