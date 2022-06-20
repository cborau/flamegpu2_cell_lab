FLAMEGPU_AGENT_FUNCTION(bcorner_move, flamegpu::MessageNone, flamegpu::MessageNone) {
  //Agent position vector
  float agent_x = FLAMEGPU->getVariable<float>("x");
  float agent_y = FLAMEGPU->getVariable<float>("y");
  float agent_z = FLAMEGPU->getVariable<float>("z");
  int agent_id = FLAMEGPU->getVariable<int>("id");
  
  //TODO: implement clapping to boundaries and elastic behaviour?

  return flamegpu::ALIVE;
}