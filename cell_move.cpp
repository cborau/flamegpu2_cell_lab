FLAMEGPU_AGENT_FUNCTION(cell_move, flamegpu::MessageNone, flamegpu::MessageNone) {
  
  int id = FLAMEGPU->getVariable<int>("id");
  //Agent position vector
  float agent_x = FLAMEGPU->getVariable<float>("x");
  float agent_y = FLAMEGPU->getVariable<float>("y");
  float agent_z = FLAMEGPU->getVariable<float>("z");

  int DEBUG_PRINTING = FLAMEGPU->environment.getProperty<int>("DEBUG_PRINTING");

  // Agent velocity
  float agent_vx = FLAMEGPU->getVariable<float>("vx");
  float agent_vy = FLAMEGPU->getVariable<float>("vy");
  float agent_vz = FLAMEGPU->getVariable<float>("vz");
  const float CELL_SPEED_REF = FLAMEGPU->environment.getProperty<float>("CELL_SPEED_REF");
  
     
  // Mass of the cell agents (current value shared with ECM. Modify rest of paramenters instead)
  float mass = 1.0;

  //Forces acting on the agent
  float agent_fx = FLAMEGPU->getVariable<float>("fx");
  float agent_fy = FLAMEGPU->getVariable<float>("fy");
  float agent_fz = FLAMEGPU->getVariable<float>("fz");
  
  printf("ANTES cell %d: x: %2.3f, y: %2.3f, z: %2.3f vx: %2.3f, vy: %2.3f, vz: %2.3f fx: %2.3f, fy: %2.3f, fz: %2.3f\n",id, agent_x,agent_y,agent_z, agent_vx,agent_vy,agent_vz, agent_fx,agent_fy,agent_fz);

  //Get the new position and velocity: 
  // a(t) = f(t) / m;
  // v(t) = v(t-1) + a(t) * dt; 
  // x(t) = x(t-1) + v(t) * dt
  const float DELTA_TIME = FLAMEGPU->environment.getProperty<float>("DELTA_TIME");
  float rand_dir_x = 1.0;
  float rand_dir_y = 1.0;
  float rand_dir_z = 1.0;
  
  
  agent_vx += (agent_fx / mass) * DELTA_TIME;
  agent_vx += CELL_SPEED_REF * DELTA_TIME * rand_dir_x;
  agent_x += agent_vx * DELTA_TIME;
  agent_vy += (agent_fy / mass) * DELTA_TIME;
  agent_vy += CELL_SPEED_REF * DELTA_TIME * rand_dir_y;
  agent_y += agent_vy * DELTA_TIME;
  agent_vz += (agent_fz / mass) * DELTA_TIME;
  agent_vz += CELL_SPEED_REF * DELTA_TIME * rand_dir_z;
  agent_z += agent_vz * DELTA_TIME;
  
  printf("DESPUES cell %d: x: %2.3f, y: %2.3f, z: %2.3f vx: %2.3f, vy: %2.3f, vz: %2.3f fx: %2.3f, fy: %2.3f, fz: %2.3f\n",id, agent_x,agent_y,agent_z, agent_vx,agent_vy,agent_vz, agent_fx,agent_fy,agent_fz);


    //Update the agents position and velocity
  FLAMEGPU->setVariable<float>("x",agent_x);
  FLAMEGPU->setVariable<float>("y",agent_y);
  FLAMEGPU->setVariable<float>("z",agent_z);
  FLAMEGPU->setVariable<float>("vx",agent_vx);
  FLAMEGPU->setVariable<float>("vy",agent_vy);
  FLAMEGPU->setVariable<float>("vz",agent_vz);
  

  return flamegpu::ALIVE;
}