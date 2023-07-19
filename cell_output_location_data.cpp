FLAMEGPU_AGENT_FUNCTION(cell_output_location_data, flamegpu::MessageNone, flamegpu::MessageSpatial3D) {
  FLAMEGPU->message_out.setVariable<int>("id", FLAMEGPU->getVariable<int>("id"));
  FLAMEGPU->message_out.setVariable<float>("x", FLAMEGPU->getVariable<float>("x"));
  FLAMEGPU->message_out.setVariable<float>("y", FLAMEGPU->getVariable<float>("y"));
  FLAMEGPU->message_out.setVariable<float>("z", FLAMEGPU->getVariable<float>("z"));
  FLAMEGPU->message_out.setVariable<float>("vx", FLAMEGPU->getVariable<float>("vx"));
  FLAMEGPU->message_out.setVariable<float>("vy", FLAMEGPU->getVariable<float>("vy"));
  FLAMEGPU->message_out.setVariable<float>("vz", FLAMEGPU->getVariable<float>("vz"));
  FLAMEGPU->message_out.setVariable<float>("fmag", FLAMEGPU->getVariable<float>("fmag"));
  FLAMEGPU->message_out.setVariable<float>("orx", FLAMEGPU->getVariable<float>("orx"));
  FLAMEGPU->message_out.setVariable<float>("ory", FLAMEGPU->getVariable<float>("ory"));
  FLAMEGPU->message_out.setVariable<float>("orz", FLAMEGPU->getVariable<float>("orz"));
  FLAMEGPU->message_out.setVariable<float>("k_elast", FLAMEGPU->getVariable<float>("k_elast"));
  FLAMEGPU->message_out.setVariable<float>("d_dumping", FLAMEGPU->getVariable<float>("d_dumping"));
  FLAMEGPU->message_out.setVariable<float>("alignment", FLAMEGPU->getVariable<float>("alignment"));
  FLAMEGPU->message_out.setVariable<float>("radius", FLAMEGPU->getVariable<float>("radius"));
  FLAMEGPU->message_out.setVariable<int>("cycle_phase", FLAMEGPU->getVariable<int>("cycle_phase"));
  FLAMEGPU->message_out.setVariable<float>("clock", FLAMEGPU->getVariable<float>("clock"));
  FLAMEGPU->message_out.setVariable<int>("completed_cycles", FLAMEGPU->getVariable<int>("completed_cycles"));

  return flamegpu::ALIVE;
}