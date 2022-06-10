FLAMEGPU_AGENT_FUNCTION(ecm_output_grid_location_data, flamegpu::MsgNone, flamegpu::MsgArray3D) {
    FLAMEGPU->message_out.setIndex(FLAMEGPU->getVariable<uint8_t>("grid_i"), FLAMEGPU->getVariable<uint8_t>("grid_j"), FLAMEGPU->getVariable<uint8_t>("grid_k"));
    FLAMEGPU->message_out.setVariable<int>("id", FLAMEGPU->getVariable<int>("id"));
    FLAMEGPU->message_out.setVariable<float>("x", FLAMEGPU->getVariable<float>("x"));
    FLAMEGPU->message_out.setVariable<float>("y", FLAMEGPU->getVariable<float>("y"));
    FLAMEGPU->message_out.setVariable<float>("z", FLAMEGPU->getVariable<float>("z"));
    FLAMEGPU->message_out.setVariable<float>("vx", FLAMEGPU->getVariable<float>("vx"));
    FLAMEGPU->message_out.setVariable<float>("vy", FLAMEGPU->getVariable<float>("vy"));
    FLAMEGPU->message_out.setVariable<float>("vz", FLAMEGPU->getVariable<float>("vz"));
    FLAMEGPU->message_out.setVariable<uint8_t>("grid_i", FLAMEGPU->getVariable<uint8_t>("grid_i"));
    FLAMEGPU->message_out.setVariable<uint8_t>("grid_j", FLAMEGPU->getVariable<uint8_t>("grid_j"));
    FLAMEGPU->message_out.setVariable<uint8_t>("grid_k", FLAMEGPU->getVariable<uint8_t>("grid_k"));
	FLAMEGPU->message_out.setVariable<float>("concentration", FLAMEGPU->getVariable<float>("concentration"));
    return flamegpu::ALIVE;
}