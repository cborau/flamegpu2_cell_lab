FLAMEGPU_AGENT_FUNCTION(ecm_output_grid_location_data, flamegpu::MessageNone, flamegpu::MessageArray3D) {
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
	const uint8_t N_SPECIES = 2; // WARNING: this variable must be hard coded to have the same value as the one defined in the main python function. TODO: declare it somehow at compile time
	for (int i = 0; i < N_SPECIES; i++) {
		float c = FLAMEGPU->getVariable<float, N_SPECIES>("concentration_multi", i);
		FLAMEGPU->message_out.setVariable<float, N_SPECIES>("concentration_multi", i, c);
	}
	
    return flamegpu::ALIVE;
}